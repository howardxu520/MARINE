import unittest
import os
import sys

directory_path = os.path.abspath(os.path.join('../src/'))
if directory_path not in sys.path:
    sys.path.append(directory_path)
    
from read_process import get_contig_lengths_dict,\
incorporate_replaced_pos_info,incorporate_insertions_and_deletions,\
get_positions_from_md_tag,reverse_complement


class TestReadProcessFunctions(unittest.TestCase): 
    def test_reverse_complement(self):
        # Test reverse complement function
        print("Testing reverse_complement")
        test_seq = 'ACTGAC'
        expected_seq = 'GTCAGT'
        self.assertEqual(reverse_complement(test_seq), expected_seq)
    
    def test_incorporate_replaced_pos_info_insertions(self):
        # Test incorporating insertions (1 is insertion)
        print("Testing incorporate_replaced_pos_info_insertions...")
        
        test_aligned_seq = \
        'GATGCTTATATAGGGAACAAAATGGTCCCTACACCATTTTTTTTTTCTGGAGTGCATAATGGATACATTTGATGACTTTTACCCTCTTATCTAAATCTAAA'
        fixed_test_aligned_seq = \
        'GATGCTTATATAGGGAACAAAATGGTCCCTACACCATTTTTTTTTCTGGAGTGCATAATGGATACATTTGATGACTTTTACCCTCTTATCTAAATCTAAA'
        test_cigar_tuples = [(0, 44), (1, 1), (0, 56)]
        self.assertEqual(incorporate_insertions_and_deletions(test_aligned_seq, test_cigar_tuples),
                         fixed_test_aligned_seq)
        
    def test_incorporate_replaced_pos_info_deletions(self):
        # Test incorporating deletions (2 is deletion, 4 is soft-clipping)
        print("Testing incorporate_replaced_pos_info_deletions...")
        
        test_aligned_seq = \
        'TCTTTGATAGAGCCACCAAGATGCTTATATAGGGAACAAATGGTCCCTACACCATTTTTTTTCCTGGAGTGCCCCATGTACTCTGCGTTGATACCACTGCT'
        fixed_test_aligned_seq =\
        'TCTTTGATAGAGCCACCAAGATGCTTATATAGGGAAC*AAATGGTCCCTACACCATTTTTTTTCCTGGAGTGC'
        test_cigar_tuples = [(0, 37), (2, 1), (0, 35), (4, 29)]
        self.assertEqual(incorporate_insertions_and_deletions(test_aligned_seq, test_cigar_tuples),
                         fixed_test_aligned_seq)

    def test_get_positions_from_md_tag(self):
        print("Testing get_positions_from_md_tag...")
        md_tags_and_expectations = {
            '11T63': [11, 75]
        }
        
        for md_tag, expected in md_tags_and_expectations.items():
            print('\tTesting {}...'.format(md_tag))
            positions = get_positions_from_md_tag(md_tag)
            self.assertEqual(expected, positions)
            
    def test_incorporate_replaced_pos_info(self):
        print("Testing incorporate_replaced_pos_info...")
        indicated_sequence, bases_at_pos = incorporate_replaced_pos_info('ACTAGACA', [0, 3, 6])
        self.assertEqual('ActAgaCa', indicated_sequence)
        self.assertEqual(['A', 'A', 'C'], bases_at_pos)
        
        
unittest.main()