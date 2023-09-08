#!/usr/bin/env python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import multiprocessing
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pybedtools
from tqdm import trange
import sys
from collections import defaultdict
from pandas.errors import EmptyDataError
import pysam
from functools import partial
import gc
from pathlib import Path
import logging
import shutil
pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('display.max_rows', 500)

os.environ["TMPDIR"] = "/home/bay001/projects/eric_amanda_g3bp1_20230724/permanent_data/01_mpileup/temp"

ETYPES = {
    'ct': {'refp': 'C', 'altp': 'T', 'refm': 'G', 'altm': 'A'},
    'ag': {'refp': 'A', 'altp': 'G', 'refm': 'T', 'altm': 'C'},
    'ai': {'refp': 'A', 'altp': 'G', 'refm': 'T', 'altm': 'C'},
}

def build_fasta_dict(fasta):
    """returns dictionary of fasta id:sequence"""
    fasta_dict = {}
    with open(fasta, 'rU') as FastaFile:
        for rec in SeqIO.parse(FastaFile, 'fasta'):
            name = rec.id
            seq = str(rec.seq)
            fasta_dict[name] = seq
    return fasta_dict


def get_ct_edit_fraction(row, strand, etype):
    """
    Calculates the C>T (or G>A) edit fraction at a given position.
    """
    
    if strand == '+':
        ref = ETYPES[etype]['refp']
        alt = ETYPES[etype]['altp']
    elif strand == '-':
        ref = ETYPES[etype]['refm']
        alt = ETYPES[etype]['altm']
    else:
        ref, alt = 0
    if row['ref'] == ref:
        try:
            return row[alt] / (row[ref] + row[alt])
        except ZeroDivisionError:
            return -1
    return -1


def count_bases(st):
    return {
        'A': st.count('A'),
        'T': st.count('T'),
        'C': st.count('C'),
        'G': st.count('G'),
        'del': st.count('d'),
        'ins': st.count('i')
    }


def get_position_matrix(bam_file, region, reffile, etype='ct', stepper='nofilter', debug_read="", min_base_quality=0):
    # print('region: {}'.format(region))
    bam = pysam.AlignmentFile(bam_file, "rb", reference_filename=reffile)
    chrom = region.chrom
    start = region.start
    end = region.end
    strand = region.strand
    
    reference = pybedtools.BedTool.seq([chrom, 0, end], reffile)
    
    count = start
    alphabet = {}
    edits = defaultdict(dict)
    positions = []
    offset = 0
    max_offset = 0
    MAX_DEPTH = 10000000
    # progress = trange(end-start, leave=False, position=0, desc='position')
    for position in bam.pileup(chrom, start, end, stepper=stepper, max_depth=MAX_DEPTH, min_base_quality=min_base_quality):
        if position.pos >= start:
            st = ""
            if count >= end or position.pos >= end:
                break
            
            while count < position.pos:
                alphabet = {
                    'A': 0,
                    'T': 0,
                    'C': 0,
                    'G': 0,
                    'del': 0,
                    'ins': 0,
                    'ref': reference[position.reference_pos].upper()
                }
                positions.append(alphabet)
                count += 1
            
            for read in position.pileups:
                if not read.alignment.is_secondary and not read.alignment.is_supplementary:
                    if not read.is_del and not read.is_refskip:
                        try:
                            if strand == '+' and reference[count].upper() == ETYPES[etype]['refp'] and read.alignment.query_sequence[read.query_position] == ETYPES[etype]['altp']:
                                edits[read.alignment.query_name][count] = 1
                            elif strand == '-' and reference[count].upper() == ETYPES[etype]['refm'] and read.alignment.query_sequence[read.query_position] == ETYPES[etype]['altm']:
                                edits[read.alignment.query_name][count] = 1
                        except Exception as e:
                            print(e, strand, reference[count].upper(), count)
                        st += read.alignment.query_sequence[read.query_position]

                        if read.alignment.query_name == debug_read:
                            print(f"st at refpos: {count} and querypos: {read.query_position}: {st}")
                    elif read.is_del:
                        st += 'd'
                        if read.alignment.query_name == debug_read:
                            print(f"st at refpos: {count} and querypos: {read.query_position}: {st}")
                    elif read.is_refskip:
                        st += 'i'
                        if read.alignment.query_name == debug_read:
                            print(f"st at refpos: {count} and querypos: {read.query_position}: {st}")
                    else:
                        print("THIS READ", read.alignment.query_name, count, position.reference_pos, read.query_position)
                        print(read.is_del, read.is_refskip, read.indel)
                        st += '-'
                        if read.alignment.query_name == debug_read:
                            print(f"st at refpos: {count} and querypos: {read.query_position}: {st}")
                
            alphabet = count_bases(st)
            try:
                alphabet['ref'] = reference[count].upper()
            except IndexError:
                print('count: {} and alphabet: {} and position: {}'.format(count, alphabet, position.pos))
                raise
            count += 1
            positions.append(alphabet)
        # progress.update(1)
    while count < end:
        alphabet = {
            'A': 0,
            'T': 0,
            'C': 0,
            'G': 0,
            'del': 0,
            'ins': 0,
            'ref': reference[count].upper()
        }
        count += 1
        positions.append(alphabet)
    
    df = pd.DataFrame(positions)
    df.index += start
    df['edit_fraction'] = df.apply(partial(get_ct_edit_fraction, strand=strand, etype=etype), axis=1)
    df['gene'] = region.name
    df['chrom'] = chrom
    """
    If we introduce a cleanup here, we sometimes hit FileNotFound exceptions, possibly due to race condition?
    If we don't, /tmp fills up and we start clogging our CWD. 
    Solution: Define a cache/temp directory that is safe from filling up, ideally in scratch or something?
    """
    # pybedtools.helpers.cleanup(verbose=False, remove_all=False)
    return df, edits
    
def main():
    

    parser = ArgumentParser(
        description="Creates windows surrounding edit sites and calculates the read coverage across edits and editable sites."
    )
    parser.add_argument(
        "--bam", 
        help="input file", 
        required=True, 
        default=None
    )
    parser.add_argument(
        "--regions", "--regions", 
        help="BED6 (stranded) file of regions, defaults to Gencode v40 Human GRCh38. Defaults to: /projects/ps-yeolab3/bay001/annotations/GRCh38/gencode_v40/GRCh38_v40_genes.bed", 
        required=False, 
        default='/projects/ps-yeolab3/bay001/annotations/GRCh38/gencode_v40/GRCh38_v40_genes.bed',
    )
    parser.add_argument(
        "--fasta", 
        help="fasta file of genome assembly, defaults to: /projects/ps-yeolab3/bay001/annotations/GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta", 
        required=False, 
        default='/projects/ps-yeolab3/bay001/annotations/GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta',
    )
    parser.add_argument(
        "--threads", 
        help="Number of threads to use. Default: 8", 
        required=False, 
        default=8,
    )
    parser.add_argument(
        "--outdir", 
        help="Output directory. Default: bam file parent directory", 
        required=False, 
        default=None,
    )
    parser.add_argument(
        "--all", 
        help="If specified, keep ALL positions. Default: False, keep only edited positions", 
        required=False, 
        default=False,
        action='store_true',
        dest='keepall',
    )
    parser.add_argument(
        "--etype", 
        help="Specifies edit type (can be either ct, ai). Default: ct", 
        required=False, 
        default='ct',
    )
    args = parser.parse_args()
    bam_file = args.bam
    regions_file = args.regions
    ref = args.fasta
    nthreads = args.threads
    keepall = args.keepall
    outdir = os.path.dirname(bam_file) if args.outdir is None else args.outdir
    etype = args.etype
    temp_root = os.path.join(os.path.dirname(bam_file), 'tmp')
    
    # setup logger
    logger = logging.getLogger(f'sailor2')
    ih = logging.FileHandler(f'sailor2.log')
    eh = logging.FileHandler(f'sailor2.err')
    sh = logging.StreamHandler(stream=sys.stdout)

    eh.setLevel(logging.ERROR)
    logger.addHandler(eh)
    formatter = logging.Formatter(
        '[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s')
    eh.setFormatter(formatter)

    logging.basicConfig(
        level=logging.INFO, 
        format='[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s',
        handlers=[ih, eh]
    )

    logger.info("Starting sailor2")
    logger.info(f"Genome fasta: {ref}")
    logger.info(f"Regions file: {regions_file}")
    logger.info(f"Edit type: {etype}")
    df = pd.read_csv(regions_file, sep='\t', names=['chrom','start','end','name','score','strand'])
    chromosomes = sorted(list(set(df['chrom'])))
    
    # One chromosome file open at a time to reduce memory footprint?
    for chromosome in chromosomes:
        temp_dir = os.path.join(temp_root, os.path.basename(bam_file), chromosome)
        Path(temp_dir).mkdir(parents=True, exist_ok=True)
        os.environ["TMPDIR"] = temp_dir
        
        regions = pybedtools.BedTool.from_dataframe(
            df[df['chrom']==chromosome]
        ).sort()  # there may be a better way of iterating over regions within each chromosome
        logger.info(f'Processing {bam_file}: chr: {chromosome} ({len(regions)} genes found)')
        logger.info(f'Temp folder: {temp_dir}')
        
        pool = multiprocessing.Pool(nthreads)
        header = True  # only write the header once, at the beginning.
        tasks = [(bam_file, interval, ref, etype) for interval in regions]
        jobs = [pool.apply_async(get_position_matrix, args = task) for task in tasks]
        progress = trange(len(tasks), position=0, leave=False, desc=chromosome)
        
        with open(os.path.join(outdir, os.path.splitext(os.path.basename(bam_file))[0] + f'_{chromosome}.pileup.tsv.gz'), 'w') as x:
            pass
        
        with open(os.path.join(outdir, os.path.splitext(os.path.basename(bam_file))[0] + f'_{chromosome}.reads.txt'), 'w') as o:
            for job, task in zip(jobs, tasks):
                edit_table, edit_reads = job.get()
                if not keepall:
                    edit_table = edit_table[edit_table['edit_fraction']!=-1]
                edit_table.to_csv(os.path.join(outdir, os.path.splitext(os.path.basename(bam_file))[0] + f'_{chromosome}.pileup.tsv.gz'), compression='gzip', mode='a', sep='\t', index=True, header=header)
                header = False
                for read, position in edit_reads.items():
                    for pos, count in position.items():
                        o.write(f"{read}\t{chromosome}:{pos}\t{count}\n")
                progress.update(1)
        pool.close()
    
    # For some reason, removing the tmp files 
    for chromosome in chromosomes:
        temp_dir = os.path.join(temp_root, os.path.basename(bam_file), chromosome)
        logger.info(f"Removing temp files for {bam_file}: chr: {chromosome}")
        shutil.rmtree(temp_dir, ignore_errors=True)
    logger.info("Finished sailor2")
if __name__ == '__main__':
    main()
