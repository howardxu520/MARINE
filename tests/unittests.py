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
    def test_incorporate_replaced_pos_info(self):
        self.assertEqual(5, 4)
        
unittest.main()