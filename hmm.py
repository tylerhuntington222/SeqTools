"""
Top-level comment
Note: feel free to modify the starter code below
"""

import optparse
from math import log
import numpy as np
import sys

TIMES = np.array([0.32, 1.75, 4.54, 9.40])

def parse_args():
    """Parse and return command-line arguments"""

    parser = optparse.OptionParser(description='HMM for Tmrca')
    parser.add_option('-f', '--fasta_filename', type='string', help='path to input fasta file')
    parser.add_option('-p', '--param_filename', type='string', help='path to input parameter file')
    parser.add_option('-t', '--truth_filename', type='string', help='path to file of true values')
    parser.add_option('-i', '--num_iter', type='int', default=15, help='number of Baum-Welch iterations')
    parser.add_option('-o', '--out_folder', type='string', help='path to folder for output files')
    parser.add_option('-s', '--suffix', type='string', help='suffix string to include in output files')
    (opts, args) = parser.parse_args()

    mandatories = ['fasta_filename','param_filename','truth_filename','out_folder','suffix']
    for m in mandatories:
        if not opts.__dict__[m]:
            print('mandatory option ' + m + ' is missing\n')
            parser.print_help()
            sys.exit()

    return opts

def parse_params(param_filename):
    """ Parse initial parameter file to extract initial, transision, and
    emission probabilities.
    Authors: Andrew H. Chan, Sara Mathieson
    """
    param_file = open(param_filename,'r')
    param_lines = param_file.readlines()
    param_file.close()
    K = 4

    # parse init state
    init = np.array([float(l.split()[1]) for l in param_lines[2:2+K]])
    # parse transition matrix
    tran = np.array([[float(x) for x in l.split()] for l in param_lines[11:11+K]])
    # parse in emit distribution
    emit = np.array([[float(l.split()[1]), float(l.split()[2])] for l in param_lines[19:19+K]])

    # convert to log-space
    log_initial = np.array([log(x) for x in init])
    log_transition = np.array([[log(x) for x in row] for row in tran])
    log_emission = np.array([[log(x) for x in row] for row in emit])

    return log_initial, log_transition, log_emission

def display_params(params):
    """Create a parameter string to write to a file (for Part 2).
    Authors: Andrew H. Chan, Sara Mathieson
    """
    init, tran, emit = params

    param_str = '# Initial probabilities\n'
    for i in range(len(init)):
        param_str += str("%.6e" % init[i]) + '\n'

    param_str += '\n# Transition Probabilities\n'
    for i in range(len(tran)):
        row = ''
        for j in range(len(tran)):
            row += str("%.6e" % tran[i][j]) + ' '
        param_str += row + '\n'

    param_str += '\n# Emission Probabilities\n'
    for i in range(len(emit)):
        row = ''
        for j in range(2):
            row += str("%.6e" % emit[i][j]) + ' '
        param_str += row + '\n'

    return param_str

def main():

    # parse commandline arguments
    opts = parse_args()
    param_filename = opts.param_filename

    # read parameter file
    init_params = parse_params(opts.param_filename)
    print(init_params)

if __name__ == "__main__":
  main()
