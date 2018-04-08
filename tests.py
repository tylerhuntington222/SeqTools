"""
tests.py

PURPOSE:
Unit tests for Viterbi and Forward-Backward Algorithm classes and main programs.

Tyler Huntington, Nathan Holeman
CS68 lab08
April, 2018
"""

# imports
import unittest
from ViterbiClass import ViterbiClass
from viterbi import parse_cl_args, load_observed_data, load_params

#############################################################################
# --------------------------Viterbi Main Program Tests ---------------------#
#############################################################################
class ViterbiMainProgram(unittest.TestCase):

    def test_load_observed_data_function_works(self):
        fasta_file = "data/test_case.fasta"
        expected = "1001010011000101001000011010000001000011"

        res = load_observed_data(fasta_file)
        self.assertEqual(expected, res)




if (__name__ == '__main__'):
    unittest.main()
