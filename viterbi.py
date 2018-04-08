"""
Top level comment: be sure to include the purpose/contents of this file
as well as the author(s)
"""
import optparse
import sys
import os

def main():
    # parse command line args
    (fasta_file,
    param_file,
    true_tmrca,
    train_iters,
    out_dir,
    suffix) = parse_cl_args()

    # load observed data file

    # load paramaters file

    # pass observed data and params to new Viterbi object

    # calculate highest probability path

def parse_cl_args():

    # set up command line argument parser
    desc   = "A program to perform Viterbi's Algorithm."
    parser = optparse.OptionParser(description=desc)

    parser.add_option('-f', '--observed_data_fasta_file', type='string',
    help='Filepath of observed data in FASTA file format.')

    parser.add_option('-p', '--input_parameter_file',       type='string',
    help='Filepath of input parameter file.')

    parser.add_option('-t', '--true_tmrca_vals_file',       type='string',
    help='Filepath of true TMRCA values file.')

    parser.add_option('-i', '--num_training_iterations',       type='int',
    help='Number of training iterations for parameter estimation.')

    parser.add_option('-o', '--output_directory',       type='string',
    help='Target directory for output files.')

    parser.add_option('-s', '--output_file_suffix',       type='string',
    help='Suffix for output files.')

    (opts, args) = parser.parse_args()

    mandatories  = ['observed_data_fasta_file',
                    'input_parameter_file',
                    'true_tmrca_vals_file',
                    'num_training_iterations',
                    'output_directory',
                    'output_file_suffix',
                    ]

    for m in mandatories:
        if not opts.__dict__[m]:
            print('mandatory option ' + m + ' is missing\n')
            parser.print_help()
            sys.exit()

    f = opts.observed_data_fasta_file
    p = opts.input_paramater_file
    t = opts.true_tmrca_vals_file
    i = opts.n_training_iterations
    o = os.path.join(opts.output_directory)
    s = os.path.join(opts.output_file_suffix)

    return (f, p, t, n, o, s)

def load_observed_data(fasta_file):
    pass

def load_params(param_file):
    pass

if __name__ == '__main__':
    main()
