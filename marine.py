import argparse
from glob import glob
from multiprocessing import Pool
import multiprocessing
import os
import pandas as pd
import polars as pl
import psutil
import pysam
from scipy.special import betainc
import sys
from sys import getsizeof
import time
from tqdm import tqdm
import tracemalloc

sys.path.append('src/')

from read_process import incorporate_replaced_pos_info,incorporate_insertions_and_deletions,\
get_positions_from_md_tag,reverse_complement,get_edit_information,get_edit_information_wrapper,\
has_edits,get_total_coverage_for_contig_at_position,\
print_read_info, get_read_information, get_hamming_distance, \
remove_softclipped_bases,find

from utils import get_intervals, index_bam, write_rows_to_info_file, write_header_to_edit_info, \
write_read_to_bam_file, remove_file_if_exists, make_folder, concat_and_write_bams_wrapper, \
pretty_print, read_barcode_whitelist_file, get_contigs_that_need_bams_written

from core import run_edit_identifier, run_bam_reconfiguration, \
gather_edit_information_across_subcontigs, run_coverage_calculator, generate_site_level_information

from annotate import annotate_sites, get_strand_specific_conversion 


def edit_finder(bam_filepath, output_folder, reverse_stranded, barcode_tag="CB", barcode_whitelist=None, contigs=[], num_intervals_per_contig=16, 
                verbose=False, cores=64, min_read_quality = 0):
    
    pretty_print("Each contig is being split into {} subsets...".format(num_intervals_per_contig))
    
    overall_label_to_list_of_contents, results, overall_time, overall_total_reads, \
    total_seconds_for_reads = run_edit_identifier(
        bam_filepath, 
        output_folder, 
        reverse_stranded=reverse_stranded,
        barcode_tag=barcode_tag,
        barcode_whitelist=barcode_whitelist,
        contigs=contigs,
        num_intervals_per_contig=num_intervals_per_contig, 
        verbose=verbose,
        cores=cores,
        min_read_quality=min_read_quality
    )
    
    #print(overall_label_to_list_of_contents.keys())
    #print(overall_label_to_list_of_contents.get(list(overall_label_to_list_of_contents.keys())[0]))
    
    pretty_print(
        [
            "Reads processed:\t{}".format(overall_total_reads), 
            "Time to process reads in min:\t{}".format(round(overall_time/60, 5))
        ],
        style="-"
    )
    
    
    total_seconds_for_reads_df = pd.DataFrame.from_dict(total_seconds_for_reads, orient='index')
    total_seconds_for_reads_df.columns = ['seconds']
    total_seconds_for_reads_df['reads'] = total_seconds_for_reads_df.index
    total_seconds_for_reads_df.index = range(len(total_seconds_for_reads_df))
    
    
    return overall_label_to_list_of_contents, results, total_seconds_for_reads_df


def bam_processing(overall_label_to_list_of_contents, output_folder, barcode_tag='CB', cores=1, number_of_expected_bams=4,
                   verbose=False):
    split_bams_folder = '{}/split_bams'.format(output_folder)
    make_folder(split_bams_folder)
    contigs_to_generate_bams_for = get_contigs_that_need_bams_written(list(overall_label_to_list_of_contents.keys()),
                                                                      split_bams_folder, 
                                                                      barcode_tag=barcode_tag,
                                                                      number_of_expected_bams=number_of_expected_bams
                                                                     )
    pretty_print("Will split and reconfigure the following contigs: {}".format(",".join(contigs_to_generate_bams_for)))
    
    
    # BAM Generation
    total_bam_generation_time, total_seconds_for_bams = run_bam_reconfiguration(split_bams_folder, bam_filepath, overall_label_to_list_of_contents, contigs_to_generate_bams_for, barcode_tag=barcode_tag, cores=cores, 
                                                                                number_of_expected_bams=number_of_expected_bams,
                                                                                verbose=verbose)
    
    total_seconds_for_bams_df = pd.DataFrame.from_dict(total_seconds_for_bams, orient='index')
    total_seconds_for_bams_df.columns = ['seconds']
    total_seconds_for_bams_df['contigs'] = total_seconds_for_bams_df.index
    total_seconds_for_bams_df.index = range(len(total_seconds_for_bams_df))
    
    return total_bam_generation_time, total_seconds_for_bams_df
    
    
import subprocess 
def coverage_processing(output_folder, barcode_tag='CB', paired_end=False, verbose=False, cores=1, number_of_expected_bams=4,
                       min_read_quality=0, bam_filepath=''):

    # Single-cell or long read version:
    edit_info_grouped_per_contig_combined = gather_edit_information_across_subcontigs(output_folder, 
                                                                                      barcode_tag=barcode_tag,
                                                                                      number_of_expected_bams=number_of_expected_bams
                                                                                     )
    
    if verbose:
        print('edit_info_grouped_per_contig_combined', edit_info_grouped_per_contig_combined.keys())

        
    results, total_time, total_seconds_for_contig = run_coverage_calculator(edit_info_grouped_per_contig_combined, output_folder,
                                                                            barcode_tag=barcode_tag,
                                                                            paired_end=paired_end,
                                                                            verbose=verbose,
                                                                            processes=cores
                                                                           )

    concat_command = 'for f in {}/coverage/*.tsv; do cat $f; done > {}/all_edit_info.tsv'.format(output_folder, output_folder)

    command_bash = '{}/concat_command.sh'.format(output_folder)
    with open(command_bash, 'w') as f:
        f.write(concat_command)
        
    print("Concatenating results...")
    subprocess.run(['bash', command_bash])
    print("Done concatenating.")
    
    total_seconds_for_contig_df = pd.DataFrame.from_dict(total_seconds_for_contig, orient='index')
    total_seconds_for_contig_df.columns = ['seconds']
    total_seconds_for_contig_df['contig sections'] = total_seconds_for_contig_df.index
    total_seconds_for_contig_df.index = range(len(total_seconds_for_contig_df))
    
    return results, total_time, total_seconds_for_contig_df
    

def print_marine_logo():
    logo_lines = [
    "::::    ::::      :::     :::::::::  ::::::::::: ::::    ::: :::::::::: ",
    "+:+:+: :+:+:+   :+: :+:   :+:    :+:     :+:     :+:+:   :+: :+:        ",
    "+:+ +:+:+ +:+  +:+   +:+  +:+    +:+     +:+     :+:+:+  +:+ +:+        ",
    "+#+  +:+  +#+ +#++:++#++: +#++:++#:      +#+     +#+ +:+ +#+ +#++:++#   ",
    "+#+       +#+ +#+     +#+ +#+    +#+     +#+     +#+  +#+#+# +#+        ",
    "#+#       #+# #+#     #+# #+#    #+#     #+#     #+#   #+#+# #+#        ",
    "###       ### ###     ### ###    ### ########### ###    #### ########## "
    ]
    for l in logo_lines:
        pretty_print(l)
        
    pretty_print("Multi-core Algorithm for Rapid Identification of Nucleotide Edits", style="=")

def calculate_sailor_score(sailor_row):
    edits = sailor_row['count']
    num_reads = sailor_row['coverage']
    original_base_counts = num_reads - edits
    
    cov_margin = 0.01
    alpha, beta = 0, 0

    num_failures = 0
    failure_messages = []
    try:
        edit_frac = float(edits) / float(num_reads)
    except Exception as e:
        num_failures += 1
        print("Bad Row: {}".format(sailor_row))
        return None
        
    # calc smoothed counts and confidence
    destination_base_smoothed = edits + alpha
    origin_base_smoothed = original_base_counts + beta

    ########  MOST IMPORTANT LINE  ########
    # calculates the confidence of theta as
    # P( theta < cov_margin | A, G) ~ Beta_theta(G, A)
    confidence = 1 - betainc(destination_base_smoothed, origin_base_smoothed, cov_margin)
    return confidence
    

def get_sailor_sites(final_site_level_information_df, conversion="C>T", skip_coverage=False):
    final_site_level_information_df = final_site_level_information_df[final_site_level_information_df['strand_conversion'] == conversion]

    if skip_coverage:
        # Case where we want to skip coverage counting, we will just set it to -1 in the sailor-style output
        coverage_col = "-1"
        weird_sites = pd.DataFrame()
    else:
        # Normally, assume that there is a coverage column
        coverage_col = final_site_level_information_df['coverage'].astype(str)

        weird_sites = final_site_level_information_df[
            (final_site_level_information_df.coverage == 0) |\
            (final_site_level_information_df.coverage < final_site_level_information_df['count'])]
    
        print("{} rows had coverage of 0 or more edits than coverage... filtering these out, but look into them...".format(
            len(weird_sites)))
              
        final_site_level_information_df = final_site_level_information_df[
        (final_site_level_information_df.coverage > 0) & \
        (final_site_level_information_df.coverage >= final_site_level_information_df['count'])
        ]

    
    final_site_level_information_df['combo'] = final_site_level_information_df['count'].astype(str) + ',' + coverage_col

    if skip_coverage:
        final_site_level_information_df['score'] = -1
    else:
        final_site_level_information_df['score'] = final_site_level_information_df.apply(calculate_sailor_score, axis=1)

    final_site_level_information_df['start'] = final_site_level_information_df['position']
    final_site_level_information_df['end'] = final_site_level_information_df['position'] + 1
    
    final_site_level_information_df = final_site_level_information_df[['contig', 'start', 'end', 'score', 'combo', 'strand']]
    return final_site_level_information_df, weird_sites
    

def collate_edit_info_shards(output_folder):
    edit_info_folder = "{}/edit_info".format(output_folder)
    edit_info_files = glob("{}/*".format(edit_info_folder))
    print("\tFound {} files in {}...".format(len(edit_info_files), edit_info_folder))
    all_dfs = []
    for f in edit_info_files:
        edit_info_df = pd.read_csv(f, sep='\t')
        all_dfs.append(edit_info_df)
        
    collated_df = pd.concat(all_dfs)
    print("\tShape of collated edit info df: {}".format(collated_df.shape))
    print("\tColumns of collated edit info df: {}".format(collated_df.columns))
    return collated_df
    
def run(bam_filepath, annotation_bedfile_path, output_folder, contigs=[], num_intervals_per_contig=16, reverse_stranded=True, barcode_tag="CB", paired_end=False, barcode_whitelist_file=None, verbose=False, coverage_only=False, filtering_only=False, annotation_only=False, sailor=False, min_base_quality = 15, min_read_quality = 0, min_dist_from_end = 10, max_edits_per_read = None, cores = 64, number_of_expected_bams=4):
    
    print_marine_logo()
    
    
    logging_folder = "{}/metadata".format(output_folder)
    make_folder(logging_folder)

    with open('{}/manifest.txt'.format(logging_folder), 'w') as f:
        f.write('bam_filepath\t{}\n'.format(bam_filepath)) 
        f.write('annotation_bedfile_path\t{}\n'.format(annotation_bedfile_path))
        f.write('output_folder\t{}\n'.format(output_folder))  
        f.write('reverse_stranded\t{}\n'.format(reverse_stranded))  
        f.write('barcode_tag\t{}\n'.format(barcode_tag))  
        f.write('barcode_whitelist_file\t{}\n'.format(barcode_whitelist_file))  
        f.write('contigs\t{}\n'.format(contigs))  
        f.write('num_intervals_per_contig\t{}\n'.format(num_intervals_per_contig))  
        f.write('verbose\t{}\n'.format(verbose))
        f.write('cores\t{}\n'.format(cores))
        f.write('number_of_expected_bams\t{}\n'.format(number_of_expected_bams))
        f.write('paired_end\t{}\n'.format(paired_end))
        f.write('min_base_quality\t{}\n'.format(min_base_quality))
        f.write('min_read_quality\t{}\n'.format(min_read_quality))
        f.write('min_dist_from_end\t{}\n'.format(min_base_quality))

    # Flag to turn of coverage calculation and column addition -- only needs to be turned off for debugging
    skip_coverage = False
        
    if not (coverage_only or filtering_only):
        if barcode_whitelist_file:
            barcode_whitelist = read_barcode_whitelist_file(barcode_whitelist_file)
        else:
            barcode_whitelist = None

        # Edit identification
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pretty_print("Identifying edits", style="~")
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        overall_label_to_list_of_contents, results, total_seconds_for_reads_df = edit_finder(
            bam_filepath, 
            output_folder, 
            reverse_stranded,
            barcode_tag,
            barcode_whitelist,
            contigs,
            num_intervals_per_contig,
            verbose,
            cores=cores,
            min_read_quality=min_read_quality
        )

        total_seconds_for_reads_df.to_csv("{}/edit_finder_timing.tsv".format(logging_folder), sep='\t')
        # REMOVE
        #sys.exit()

        if barcode_tag:
            # Make a subfolder into which the split bams will be placed
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            pretty_print("Contigs processed:\n\n\t{}".format(sorted(list(overall_label_to_list_of_contents.keys()))))
            pretty_print("Splitting and reconfiguring BAMs to optimize coverage calculations", style="~")
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
            total_bam_generation_time, total_seconds_for_bams_df = bam_processing(overall_label_to_list_of_contents, output_folder, barcode_tag=barcode_tag, cores=cores, number_of_expected_bams=number_of_expected_bams, verbose=verbose)
            total_seconds_for_bams_df.to_csv("{}/bam_reconfiguration_timing.tsv".format(logging_folder), sep='\t')
            pretty_print("Total time to concat and write bams: {} minutes".format(round(total_bam_generation_time/60, 3)))

        
    if not filtering_only:
        # Coverage calculation
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pretty_print("Calculating coverage at edited sites, minimum read quality is {}...".format(min_read_quality), style='~')

        coverage_subfolder = '{}/coverage'.format(output_folder)
        make_folder(coverage_subfolder)
        
        results, total_time, total_seconds_for_contig_df = coverage_processing(output_folder, 
                                                                               barcode_tag=barcode_tag, 
                                                                               paired_end=paired_end,
                                                                               verbose=verbose,
                                                                               cores=cores,
                                                                               number_of_expected_bams=number_of_expected_bams,
                                                                               min_read_quality=min_read_quality,
                                                                               bam_filepath=bam_filepath
                                                                              )
        total_seconds_for_contig_df.to_csv("{}/coverage_calculation_timing.tsv".format(logging_folder), sep='\t')


        pretty_print("Total time to calculate coverage: {} minutes".format(round(total_time/60, 3)))
        #all_edit_info_pd = pd.concat(results)

        #if verbose:
        #    all_edit_info_pd.to_csv('{}/all_edit_info.tsv'.format(output_folder), sep='\t')

        # Filtering and site-level information
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pretty_print("Filtering and calculating site-level statistics", style="~")
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        print("Loading edit information...")
        all_edit_info_pd = pd.read_csv('{}/all_edit_info.tsv'.format(output_folder),
                                       sep='\t', 
                                       names=['barcode', 'contig', 'position', 'ref', 'alt', 'read_id', 'strand',
       'dist_from_end', 'base_quality', 'mapping_quality', 'barcode_position',
       'coverage', 'source', 'position_barcode'], dtype={'base_quality': int, 'dist_from_end': int, 'contig': str})
        
        #all_edit_info_pd['contig'] = all_edit_info_pd['contig'].astype(str)

        # Ensure that we are taking cleaning for only unique edits
        # TODO fix thisssssss
        #all_edit_info_pd['position_barcode'] = all_edit_info_pd['position'].astype(str) + '_' + all_edit_info_pd['barcode'].astype(str)
        
        coverage_per_unique_position_df = pd.DataFrame(all_edit_info_pd.groupby(
        [
            "position_barcode"
        ]).coverage.max())

        distinguishing_columns = [
            "barcode",
            "contig",
            "position",
            "ref",
            "alt",
            "read_id",
            "strand",
            "dist_from_end",
            "base_quality",
            "mapping_quality"
        ]
        all_edit_info_unique_position_df = all_edit_info_pd.drop_duplicates(distinguishing_columns)[distinguishing_columns]

        all_edit_info_unique_position_df.index = all_edit_info_unique_position_df['position'].astype(str)\
+ '_' + all_edit_info_unique_position_df['barcode']

        all_edit_info_unique_position_with_coverage_df = all_edit_info_unique_position_df.join(coverage_per_unique_position_df)
        
        
        all_edit_info_unique_position_with_coverage_df.to_csv('{}/final_edit_info.tsv'.format(output_folder), sep='\t')

    if skip_coverage:
        # If we skip coverage steps, we have to just collate edit info results from the edit finding steps into one giant 
        # final_edit_info.tsv file (where coverage is None)
        pretty_print("Skipped coverage, so combining edit info alone...")
        all_edit_info_unique_position_with_coverage_df = collate_edit_info_shards(output_folder)

    
    # Check if filtering step finished
    final_filtered_sites_path = '{}/final_filtered_site_info.tsv'.format(output_folder)
    final_path_already_exists = False
    final_annotated_path_already_exists = False
    
    if os.path.exists(final_filtered_sites_path):
        print("{} exists...".format(final_filtered_sites_path))
        final_path_already_exists = True

    # Filtering steps
    if not final_path_already_exists:
        print("Filtering..")
        if filtering_only:
            all_edit_info_unique_position_with_coverage_df = pd.read_csv('{}/final_edit_info.tsv'.format(output_folder), sep='\t', dtype={"contig": str})

        pretty_print("Filtering edited sites", style='~')
        pretty_print("Minimum distance from end = {}, Minimum base-calling quality = {}".format(min_dist_from_end, min_base_quality))

        #all_edit_info_unique_position_with_coverage_df['min_base_quality'] = all_edit_info_unique_position_with_coverage_df['min_base_quality'].astype(int)
        
        all_edit_info_filtered = all_edit_info_unique_position_with_coverage_df[
            (all_edit_info_unique_position_with_coverage_df["base_quality"] >= min_base_quality) & 
            (all_edit_info_unique_position_with_coverage_df["dist_from_end"] >= min_dist_from_end)]
        
        
        print("Deduplicating....")
        distinguishing_columns = [
            "barcode",
            "contig",
            "position",
            "ref",
            "alt",
            "read_id",
            "strand",
            "mapping_quality",
        ]
        if not skip_coverage:
            distinguishing_columns.append("coverage")
            
        all_edit_info_filtered_deduped = all_edit_info_filtered.drop_duplicates(
            distinguishing_columns)[distinguishing_columns]
        
        pretty_print("\tNumber of edits before filtering:\n\t{}".format(len(all_edit_info_unique_position_with_coverage_df)))
        pretty_print("\tNumber of edits after filtering:\n\t{}".format(len(all_edit_info_filtered)))
        pretty_print("\tNumber of edits after deduplicating:\n\t{}".format(len(all_edit_info_filtered_deduped)))

        if max_edits_per_read:
            pretty_print("\t\tFiltering out reads with more than {} edits...".format(max_edits_per_read))
            read_edits = all_edit_info_filtered_deduped.groupby('read_id').count().sort_values('barcode')
            all_edit_info_filtered_deduped = all_edit_info_filtered_deduped[all_edit_info_filtered_deduped.read_id.isin(read_edits[read_edits['barcode'] <= max_edits_per_read].index)]
        
        all_edit_info_filtered_deduped.to_csv('{}/final_filtered_edit_info.tsv'.format(output_folder), 
                                         sep='\t', index=False)
        
            
        all_edit_info_unique_position_with_coverage_deduped_df = pd.read_csv('{}/final_filtered_edit_info.tsv'.format(output_folder), sep='\t', dtype={"contig": str})
            
    
        all_edit_info_filtered_pl = pl.from_pandas(all_edit_info_unique_position_with_coverage_deduped_df)

        final_site_level_information_df = generate_site_level_information(all_edit_info_filtered_pl, skip_coverage=skip_coverage)
        pretty_print("\tNumber of unique edit sites:\n\t{}".format(len(final_site_level_information_df)))
        pretty_print("Writing sites...\n")
        final_site_level_information_df.write_csv('{}/final_filtered_site_info.tsv'.format(output_folder), 
                                                  separator='\t')

        
        pretty_print("Adding strand-specific conversion...\n")
        final_site_level_information_df = pd.read_csv('{}/final_filtered_site_info.tsv'.format(output_folder), 
                                                  sep='\t')
        final_site_level_information_df['strand_conversion'] = final_site_level_information_df.apply(get_strand_specific_conversion, args=(reverse_stranded,), axis=1)
        final_site_level_information_df.to_csv('{}/final_filtered_site_info.tsv'.format(output_folder), 
                                                  sep='\t', index=False)
        final_path_already_exists = True
        
    if not annotation_bedfile_path:
        print("annotation_bedfile_path argument not provided ...\
        not annotating with feature information and strand-specific conversions.")
        
    if final_path_already_exists and annotation_bedfile_path:
        final_site_level_information_df = pd.read_csv('{}/final_filtered_site_info.tsv'.format(output_folder), 
                                                  sep='\t')
        final_site_level_information_annotated_df = annotate_sites(final_site_level_information_df,
                                                                   annotation_bedfile_path)
        final_site_level_information_annotated_df.to_csv('{}/final_filtered_site_info_annotated.tsv'.format(output_folder), 
                                                  sep='\t', index=False)
        final_annotated_path_already_exists = True

    
    if final_annotated_path_already_exists:
        final_annotated_site_level_information_df = pd.read_csv('{}/final_filtered_site_info_annotated.tsv'.format(output_folder), 
                                                  sep='\t')
        
        if sailor:
            print("{} sites being converted to SAILOR format...".format(len(final_site_level_information_df)))

            # Output SAILOR-formatted file for use in FLARE downstream
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # 1       629275  629276  0.966040688     2,30    +
            # 1       629309  629310  2.8306e-05      1,1043  +
    
            conversion = 'C>T'
            sailor_sites,weird_sites = get_sailor_sites(final_annotated_site_level_information_df, conversion, skip_coverage=skip_coverage)
            sailor_sites = sailor_sites.drop_duplicates()

            print("{} final deduplicated SAILOR-formatted sites".format(len(sailor_sites)))
            sailor_sites.to_csv('{}/sailor_style_sites_{}.bed'.format(
                output_folder, 
                conversion.replace(">", "-")), 
                header=False,
                index=False,       
                sep='\t')
    
            weird_sites.to_csv('{}/problematic_sites.tsv'.format(output_folder), sep='\t')

    # Check memory usage
    current, peak = tracemalloc.get_traced_memory()

    with open('{}/manifest.txt'.format(logging_folder), 'a+') as f:
        f.write(f'peak_memory_mb: {peak/1e6}\n') 
        f.write(f'time_elapsed_seconds\t{time.time()-start_time:.2f}s\n') 

    print(f"Current memory usage {current/1e6}MB; Peak: {peak/1e6}MB")
    print(f'Time elapsed: {time.time()-start_time:.2f}s')

    
    pretty_print("Done!", style="+")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MARINE')
            
    parser.add_argument('--bam_filepath', type=str, default=None)
    parser.add_argument('--annotation_bedfile_path', type=str, default=None)

    parser.add_argument('--output_folder', type=str, default=None)
    
    parser.add_argument('--barcode_whitelist_file', type=str, default=None)
    
    parser.add_argument('--cores', type=int, default=multiprocessing.cpu_count())
    
    parser.add_argument('--reverse_stranded', dest='reverse_stranded', action='store_true')

    parser.add_argument('--coverage', dest='coverage_only', action='store_true')
    parser.add_argument('--filtering', dest='filtering_only', action='store_true')
    parser.add_argument('--annotation', dest='annotation_only', action='store_true')

    parser.add_argument('--barcode_tag', type=str, default=None)
    
    parser.add_argument('--min_dist_from_end', type=int, default=0)

    parser.add_argument('--min_base_quality', type=int, default=15)
    parser.add_argument('--contigs', type=str, default='all')
    parser.add_argument('--min_read_quality', type=int, default=0)
    
    parser.add_argument('--sailor', dest='sailor', action='store_true')
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    parser.add_argument('--paired_end', dest='paired_end', action='store_true')
    parser.add_argument('--max_edits_per_read', type=int, default=None)
    parser.add_argument('--num_intervals_per_contig', type=int, default=200)
    
    args = parser.parse_args()
    bam_filepath = args.bam_filepath
    annotation_bedfile_path = args.annotation_bedfile_path
    output_folder = args.output_folder
    barcode_whitelist_file = args.barcode_whitelist_file
    cores = args.cores
    reverse_stranded = args.reverse_stranded
    contigs = args.contigs
    annotation_bedfile_path = args.annotation_bedfile_path
    
    coverage_only = args.coverage_only
    filtering_only = args.filtering_only
    annotation_only= args.annotation_only
    
    sailor = args.sailor
    verbose = args.verbose
    paired_end = args.paired_end
    
    barcode_tag = args.barcode_tag
    min_base_quality = args.min_base_quality
    min_read_quality = args.min_read_quality
    min_dist_from_end = args.min_dist_from_end
    max_edits_per_read = args.max_edits_per_read
    
    num_intervals_per_contig = args.num_intervals_per_contig
    
    if cores is None:
        cores = 16
    pretty_print("Assuming {} cores available for multiprocessing".format(cores))
   
    
    assert(not(coverage_only and filtering_only))
    
    pretty_print(["Arguments:",
                  "\tBAM filepath:\t{}".format(bam_filepath), 
                  "\tAnnotation bedfile filepath:\t{}".format(annotation_bedfile_path),
                  "\tOutput folder:\t{}".format(output_folder),
                  "\tBarcode whitelist:\t{}".format(barcode_whitelist_file),
                  "\tReverse Stranded:\t{}".format(reverse_stranded),
                  "\tBarcode Tag:\t{}".format(barcode_tag),
                  "\tPaired End:\t{}".format(paired_end),
                  "\tCoverage only:\t{}".format(coverage_only),
                  "\tFiltering only:\t{}".format(filtering_only),
                  "\tAnnotation only:\t{}".format(annotation_only),
                  "\tSailor outputs:\t{}".format(sailor),
                  "\tMinimum base quality:\t{}".format(min_base_quality),
                  "\tMinimum read quality:\t{}".format(min_read_quality),
                  "\tMinimum distance from end:\t{}".format(min_dist_from_end),
                  "\tmax_edits_per_read:\t{}".format(max_edits_per_read),
                  "\tContigs:\t{}".format(contigs),
                  "\tNumber of intervals:\t{}".format(num_intervals_per_contig),
                  "\tCores:\t{}".format(cores),
                  "\tVerbose:\t{}".format(verbose)
                 ])

    # Whether to only run for certain contigs 
    if contigs == 'all':
        contigs = []
    else:
        contigs = contigs.split(",")

    start_time = time.time()
    tracemalloc.start()
    
    run(bam_filepath, 
        annotation_bedfile_path,
        output_folder, 
        contigs=contigs,
        reverse_stranded=reverse_stranded,
        barcode_tag=barcode_tag,
        paired_end=paired_end,
        barcode_whitelist_file=barcode_whitelist_file,
        num_intervals_per_contig=num_intervals_per_contig,
        coverage_only=coverage_only,
        filtering_only=filtering_only,
        annotation_only=annotation_only,
        sailor=sailor,
        min_base_quality = min_base_quality, 
        min_read_quality = min_read_quality,
        min_dist_from_end = min_dist_from_end,
        max_edits_per_read = max_edits_per_read,
        cores = cores,
        verbose = verbose,
        number_of_expected_bams=num_intervals_per_contig
       )
    
    