import pybedtools
import pandas as pd

import pybedtools
import pandas as pd


def make_bedtool_from_final_sites(df):
    df_bed_cols = df[['contig', 'position', 'position', 'site_id', 'conversion', 'strand']]
    return pybedtools.BedTool.from_dataframe(df_bed_cols)


def get_feature_specific_conversion(r):
    ref_alt_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    conversion = r.conversion
    ref = r.ref
    alt = r.alt
    feature_strand = r.feature_strand

    if feature_strand == '+':
        return conversion
    elif feature_strand == '-':
        return '{}>{}'.format(
            ref_alt_dict.get(ref),
            ref_alt_dict.get(alt)
        )


def annotate_sites(sites_df, annotation_bedfile_path):
    annotation_bedtool = pybedtools.BedTool(annotation_bedfile_path)
    sites_bedtool = make_bedtool_from_final_sites(sites_df)

    print("Annotating sites with GTF information from {}...".format(annotation_bedfile_path))
    annotation_intersect = sites_bedtool.intersect(annotation_bedtool, wb=True, loj=True).to_dataframe()
    new_cols = ['contig', 'position', 'position', 'site_id', 'conversion', 'strand',
     'feature_chrom', 'feature_start', 'feature_end', 'feature_name', 'feature_strand', 'feature_type']
    annotation_intersect.columns = new_cols
    annotation_intersect = annotation_intersect[['site_id', 'feature_name', 'feature_strand', 'feature_type']]
    sites_with_annot_df = sites_df.merge(annotation_intersect, on='site_id')

    print("Adjusting conversion to feature strand...")
    sites_with_annot_df['feature_conversion'] = sites_with_annot_df.apply(get_feature_specific_conversion, axis=1)
    return sites_with_annot_df
