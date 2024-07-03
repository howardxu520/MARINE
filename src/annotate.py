import pybedtools
import pandas as pd

import pybedtools
import pandas as pd


def make_bedtool_from_final_sites(df):
    df_bed_cols = df[['contig', 'position', 'position', 'site_id', 'conversion', 'strand']]
    return pybedtools.BedTool.from_dataframe(df_bed_cols)


def get_strand_specific_conversion(r, reverse_stranded):
    ref_alt_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    ref = r.ref
    alt = r.alt
    mapped_strand = r.strand

    """
    if reverse_stranded:
        if mapped_strand == '+':
            mapped_strand = '-'
        elif mapped_strand == '-':
            mapped_strand = '+'
    """
    
    if mapped_strand == '+':
        return '{}>{}'.format(
            ref,
            alt
        )
    elif mapped_strand == '-':
        return '{}>{}'.format(
            ref_alt_dict.get(ref),
            ref_alt_dict.get(alt)
        )


def annotate_sites(sites_df, annotation_bedfile_path):    
    annotation_bedtool = pybedtools.BedTool(annotation_bedfile_path)
    sites_bedtool = make_bedtool_from_final_sites(sites_df)
    
    print("Annotating sites with GTF information from {}...".format(annotation_bedfile_path))
    annotation_intersect = sites_bedtool.intersect(annotation_bedtool, wb=True, loj=True, s=True).to_dataframe()
    #print(annotation_intersect)
    #print('num rows with annotation: {}'.format(len(annotation_intersect)))
    new_cols = ['contig', 'position', 'position', 'site_id', 'conversion', 'strand',
     'feature_chrom', 'feature_start', 'feature_end', 'feature_name', 'feature_type', 'feature_strand']
    annotation_intersect.columns = new_cols
    annotation_intersect = annotation_intersect[['site_id', 'feature_name', 'feature_strand', 'feature_type']]
    annotation_intersect = annotation_intersect.replace(-1, '.')
    
    annotation_intersect_agg = annotation_intersect.groupby('site_id').agg({'feature_name': ','.join, 'feature_type': ','.join, 'feature_strand': ','.join})
                                                             
    sites_with_annot_df = sites_df.merge(annotation_intersect_agg, on='site_id')

    return sites_with_annot_df
