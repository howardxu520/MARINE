import pandas as pd
import sys
import os
from glob import glob


test_name_to_expectations = {
    "edge_case_test": {
        "folder": "singlecell_tests",
        "expectations": [{
            "contig": "10",
            "barcode": 	"AAACGAATCATTCATC-1",
            "position": 132330481,
            "num_rows": 1,
            "count": 1,
            "coverage": 1,
            "conversion": "T>G",
            "strand_conversion": "A>C",
            "strand": "-",
            "feature_name": "STK32C",
            "feature_strand": "-"
        },
        {
            "contig": "10",
            "barcode": 	"AAACGAATCATTCATC-1",
            "position": 132330438,
            "num_rows": 1,
            "count": 1,
            "coverage": 1,
            "conversion": "T>C",
            "strand_conversion": "A>G",
            "strand": "-",
            "feature_name": "STK32C",
            "feature_strand": "-"
        }
                        ]
    },
    "unstranded_pair_test": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "Citrine.dna",
            "position": 435,
            "count": 22,
            "coverage": 22,
            "num_rows": 1,
            "strand_conversion": "C>T",
            "strand": "+"
        }]
    },
    "pair_test": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "18",
            "position": 49491556,
            "count": 2,
            "coverage": 2,
            "num_rows": 1,
            "strand_conversion": "C>T",
            "strand": "-",
            "feature_name": "RPL17-C18orf32,RPL17",
            "feature_strand": "-,-"
        },
        {
            "contig": "18",
            "position": 49567494,
            "count": 2,
            "coverage": 2,
            "num_rows": 1,
            "strand_conversion": "C>T",
            "strand": "+",
            "feature_name": "LIPG",
            "feature_strand": "+"
        }]
    },
    
    "F1R2_pair_test": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "chr17",
            "position": 43044352,
            "count": 1,
            "coverage": 1,
            "conversion": "G>A",
            "num_rows": 1,
            "conversion": "G>A",
            "strand_conversion": "C>T",
            "strand": "-",
            "feature_name": "BRCA1",
            "feature_strand": "-"
        }]
    },
    
    "F1R2_pair_test-single_end_mode_sailor": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "chr17",
            "position": 43044352,
            "count": 1,
            "coverage": 2,
            "conversion": "G>A",
            "num_rows": 1,
            "conversion": "G>A",
            "strand_conversion": "C>T",
            "strand": "-",
            "feature_name": "BRCA1",
            "feature_strand": "-"
        }]
    },

    "F1R2_pair_test-single_end_mode": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "chr17",
            "position": 43044352,
            "count": 1,
            "coverage": 2,
            "conversion": "G>A",
            "num_rows": 1,
            "conversion": "G>A",
            "strand_conversion": "C>T",
            "strand": "-",
            "feature_name": "BRCA1",
            "feature_strand": "-"
        }]
    },

     "F2R1_end_second_in_pair_test": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "chr17",
            "position": 43001716,
            "count": 1,
            "coverage": 1,
            "conversion": "G>A",
            "strand_conversion": "G>A",
            "strand": "+",
            "feature_name": "RPL27"
        }]
    },
    "same_pos_dif_reads_test": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "chr17",
            "position": 83199872,
            "count": 9,
            "coverage": 9,
            "conversion": "C>G",
            "strand_conversion": "C>G",
            "strand": "+",
            "feature_name": "AC139099.2"
        }]
    },
    "tax1bp3_chr17_3665556_read_test": {
        "folder": "strandedness_tests",
        "expectations": [{
            "contig": "chr17",
            "position": 3665556,
            "num_rows": 1,
            "count": 1,
            "coverage": 1,
            "conversion": "G>A",
            "strand_conversion": "G>A",
            "strand": "+",
            #"feature_name": "AC139099.2"
        }]
    },

    "only_5_cells_test": {
        "folder": "singlecell_tests",
        "expectations": [{
            "contig": "9",
            "barcode": 	"GGGACCTTCGAGCCAC-1",
            "position": 3000524,
            "num_rows": 1,
            "count": 1,
            "coverage": 12,
            "conversion": "C>A",
            "strand_conversion": "G>T",
            "strand": "-"
        },
        {
            "contig": "9",
            "barcode": 	"GGGACCTTCGAGCCAC-1",
            "position": 3000525,
            "num_rows": 1,
            "count": 1,
            "coverage": 12,
            "conversion": "C>T",
            "strand_conversion": "G>A",
            "strand": "-"
        },
        {
            "contig": "9",
            "barcode": 	"GATCCCTCAGTAACGG-1",
            "position": 3000525,
            "num_rows": 1,
            "count": 1,
            "coverage": 4,
            "conversion": "C>G",
            "strand_conversion": "G>C",
            "strand": "-"
        }]
    },
    "long_read_sc_test": {
        "folder": "singlecell_tests",
        "expectations": [{
            "contig": "6",
            "barcode": 	"ENSMUST00000203816-AACGTGTTGGAGAGGG-16-G",
            "position": 115807969,
            "num_rows": 1,
            "count": 1,
            "coverage": 1,
            "conversion": "A>C",
            "strand_conversion": "T>G",
            "strand": "-",
            "feature_name": "Rpl32"
        },
        {
            "contig": "6",
            "barcode": 	"ENSMUST00000081840-AAGTCGTACCAGGCTC-40-C",
            "position": 115805653,
            "num_rows": 1,
            "count": 1,
            "coverage": 1,
            "conversion": "G>A",
            "strand_conversion": "C>T",
            "strand": "-",
            "feature_name": "Rpl32"
        },
        {
            "contig": "6",
            "barcode": 	"ENSMUST00000081840-AACGTGTTGGAGAGGG-40-G",
            "position": 115807015,
            "num_rows": 1,
            "count": 1,
            "coverage": 8,
            "conversion": "C>T",
            "strand_conversion": "G>A",
            "strand": "-",
            "feature_name": "Rpl32"
        }]
    }

}

print("Current directory: {}".format(os.getcwd()))

# Check results of each test
failures = 0
for test_name, info in test_name_to_expectations.items():        
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nChecking results for {}".format(test_name))
    
    expectations_list = info.get("expectations")
    for expectations in expectations_list:
        print("\tExpecting: {}".format(expectations))
    
              
        folder = info.get("folder")
        
        contig = expectations.get("contig")
        barcode = expectations.get("barcode", None)
        
        position = expectations.get("position")
        
        final_filtered_site_info_annotated = "{}/{}/final_filtered_site_info_annotated.tsv".format(folder, test_name)
        final_filtered_site_info_annotated_df = pd.read_csv(final_filtered_site_info_annotated, sep='\t', index_col=0)
    
        row_of_interest = final_filtered_site_info_annotated_df[
            (final_filtered_site_info_annotated_df['position'] == position) &\
            (final_filtered_site_info_annotated_df['contig'].astype(str) == contig)
        ]
    
    
        if barcode:
            row_of_interest = row_of_interest[row_of_interest.barcode == barcode]
    
        failure = False
        try:
            assert(len(row_of_interest) == expectations.get("num_rows", 1))
        except Exception as e:
            print("Num rows expected: {}, was {}".format(expectations.get("num_rows", 1), len(row_of_interest)))
            failure = True
            
        for attribute in list(expectations.keys()):
            if attribute in ['count', 'coverage', 'conversion', 'strand', 'feature_name']:
                attribute_expectation = expectations.get(attribute)
                try:
                    assert(row_of_interest[attribute].iloc[0] == attribute_expectation)
                except Exception as e:
                    print("Exception: {} was {}".format(attribute, row_of_interest[attribute].iloc[0]))
                    failure = True
        if not failure:
            print("\n\t >>> {} passed! <<<\n".format(test_name))
        else:
            print("\n\t ~~~ {} FAILED! ~~~\n".format(test_name))
            failures += 1



tests_dir = sys.argv[1]
print('tests dir is {}'.format(tests_dir))
# Single-cell vs bulk processing check for same single-cell dataset 
sc_5_cells_path = "singlecell_tests/only_5_cells_test/final_filtered_site_info.tsv"
bulk_5_cells_path = "singlecell_tests/only_5_cells_bulk_mode_test/final_filtered_site_info.tsv"

print("Current directory: {}".format(os.getcwd()))
print("\tMatching files: {}".format(glob('singlecell_tests/only_5*/*final_filtered_site_info.tsv')))

sc_5_cells = pd.read_csv(sc_5_cells_path, sep='\t').sort_values(['position', 'strand_conversion'])
bulk_5_cells = pd.read_csv(bulk_5_cells_path, sep='\t').sort_values(['position', 'strand_conversion'])

print("Checking that analyzing a single-cell dataset in 'bulk' mode (i.e. not specificying the 'CB' barcode) yields the exact same positions and base changes, but with counts and coverages aggregated rather than at a single-cell resolution")
grouped_sc = pd.DataFrame(sc_5_cells.groupby(['contig', 'position', 'strand_conversion']).agg({'count': sum, 'strand_conversion': 'unique'}))
grouped_sc.index.names = ['contig', 'position', 'c']
grouped_sc['strand_conversion'] = [i[0] for i in grouped_sc['strand_conversion']]
grouped_sc = grouped_sc.sort_values(['position', 'strand_conversion'])

grouped_sc_rows = []
for r in grouped_sc.iterrows():
    grouped_sc_rows.append(r[0])

bulk_rows = []
for r in bulk_5_cells.iterrows():
    r = r[1]
    bulk_rows.append((r['contig'], r['position'], r['strand_conversion']))

try:
    assert(grouped_sc_rows == bulk_rows)
    for bulk_item, grouped_sc_item in zip(bulk_rows, grouped_sc_rows):
        assert(bulk_item == grouped_sc_item)
    assert(len(grouped_sc_rows) == len(bulk_rows))
    print("\n\t >>> single-cell and bulk on same dataset comparison passed! <<<\n")
except Exception as e:
    print("Exception: {}".format(e))
    print("\n\t ~~~ single cell vs bulk modes on sc dataset equivalency test FAILED! ~~~\n")
    failures += 1


print("There were {} failures".format(failures))
if failures > 0:
    sys.exit(1)
else:
    sys.exit(0)
