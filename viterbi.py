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
    params_file,
    true_tmrca,
    train_iters,
    out_dir,
    suffix) = parse_cl_args()

    # load observed data file
    #observed = load_observed_data(fasta_file)


    # load paramaters file
    p1, p2, p3 = load_params("data/initial_parameters_mu.txt")

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

    # mandatories  = ['observed_data_fasta_file',
    #                 'input_parameter_file',
    #                 'true_tmrca_vals_file',
    #                 'num_training_iterations',
    #                 'output_directory',
    #                 'output_file_suffix',
    #                 ]
    #TEMP: disable mandatory CL args
    mandatories = []

    for m in mandatories:
        if not opts.__dict__[m]:
            print('mandatory option ' + m + ' is missing\n')
            parser.print_help()
            sys.exit()

    f = opts.observed_data_fasta_file
    p = opts.input_parameter_file
    t = opts.true_tmrca_vals_file
    i = opts.num_training_iterations
    o = opts.output_directory
    s = opts.output_file_suffix

    return (f, p, t, i, o, s)

def load_observed_data(fasta_file):

    seqs = []
    cur_seq = ""

    # load the sequences into a list
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if cur_seq:
                    seqs.append(cur_seq)
                    cur_seq = ""
            else:
                cur_seq += line.strip()
    # append the last sequence
    seqs.append(cur_seq)

    # collapse the two sequences into string of 1's (matches) and 0's (mismatch)
    res = ""
    for i in range(len(seqs[0])):
        res += "1" if seqs[0][i] == seqs[1][i] else "0"

    return(res)



def load_params(params_file):
    p_init = {}
    p_trans = {}
    p_emit = {}

    with open(params_file, 'r') as f:
        for line in f:
            # print(line)
            if line.strip() == "# Initial probabilities":
                next(f)
                record = f.readline().strip().split()
                while record != []:
                    print(record)
                    print(record[0])
                    print(record[1])
                    p_init[record[0]] = record[1]
                    record = f.readline().strip().split()
            if line.strip() == "# Transition Probabilities":
                next(f)
                next(f)
                row = 0
                record = f.readline().strip().split()
                while record != []:
                    print(record)
                    p_trans[row] = {}
                    for i in range(len(record)):
                        p_trans[row][i] = record[i]
                    record = f.readline().strip().split()
                    row += 1
            if line.strip() == '# Emission Probabilities':
                next(f)
                record = f.readline().strip().split()
                while record != []:
                    print(record)
                    p_emit[row] = {}
                    for i in range(len(record)):
                        p_emit[row][i] = record[i]
                    record = f.readline().strip().split()
                    row += 1
            print(p_emit)


    return (1,1,1)


if __name__ == '__main__':
    main()
