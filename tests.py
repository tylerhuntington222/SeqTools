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

    def test_load_param_function_works(self):
        params_file = "data/initial_parameters_mu.txt"

        p_init, p_trans, p_emit = load_params(params_file)

        # check initial probabilities
        self.assertEqual(0.603154, p_init[0])
        self.assertEqual(0.357419, p_init[1])
        self.assertEqual(0.0388879, p_init[2])
        self.assertEqual(0.000539295, p_init[3])

        # check transition probabilities
        expected_transitions = [
            0.999916,       0.0000760902,    8.27877e-6,  1.14809e-7,
            0.000128404,    0.999786,       0.0000848889, 1.17723e-6,
            0.000128404,    0.000780214,    0.999068,     0.0000235507,
            0.000128404,    0.000780214,    0.00169821,   0.997393,
        ]

        counter = 0
        # iterates through transition matrix in row major order
        for k in range(3):
            for l in range(3):
                self.assertEqual(expected_trans[counter], p_trans[k][l])
                counter += 1


if (__name__ == '__main__'):
    unittest.main()
