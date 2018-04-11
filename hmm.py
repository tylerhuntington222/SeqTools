"""
Top-level comment
Note: feel free to modify the starter code below
"""

import optparse
from math import log
from viterbi import ViterbiClass
from baum_welch import ForwardBackwardClass
import sys
import matplotlib.pyplot as plt

# TIMES = np.array([0.32, 1.75, 4.54, 9.40]) # replaced with decodings dict
                                             # in ViterbiClass.py

def parse_args():

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

def parse_params(param_filename):
    """ Parse initial parameter file to extract initial, transision, and
    emission probabilities.
    Authors: Andrew H. Chan, Sara Mathieson
    """
    p_init = {}
    p_trans = {}
    p_emit = {}

    with open(param_filename, 'r') as f:
        for line in f:
            # print(line)
            if line.strip() == "# Initial probabilities":
                next(f)
                record = f.readline().strip().split()
                while record != []:
                    p_init[int(record[0])] = float(record[1]  )
                    record = f.readline().strip().split()
            if line.strip() == "# Transition Probabilities":
                next(f)
                next(f)
                row = 0
                record = f.readline().strip().split()
                while record != []:
                    p_trans[row] = {}
                    for i in range(len(record)):
                        p_trans[row][i] = float(record[i])
                    record = f.readline().strip().split()
                    row += 1
            if line.strip() == '# Emission Probabilities':
                next(f)
                record = f.readline().strip().split()
                while record != []:
                    p_emit[int(record[0])] = {}
                    for i in range(1, len(record)):
                        p_emit[int(record[0])][i-1] = float(record[i])
                    record = f.readline().strip().split()
    return (p_init, p_trans, p_emit)

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
        res += "0" if seqs[0][i] == seqs[1][i] else "1"

    return(res)

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

def write_decoding_to_file(step_data, outfile):
    with open(outfile, 'w') as f:
        f.write("# Viterbi_decoding posterior_decoding posterior_mean")
        f.write("\n")
        for vit_dec, post_dec, post_mean in step_data:
            line = " ".join([str(vit_dec), str(post_dec), \
            "%.6f" % (round(post_mean, 6))])
            f.write(line)
            f.write("\n")

def make_graph(truth,viterbi,decoding,mean,suffix):
    plt.plot([i for i in range(len(truth))], truth, '#000000', linewidth=2)
    plt.plot([i for i in range(len(viterbi))], viterbi, '#006600')
    plt.plot([i for i in range(len(decoding))], decoding, '#0000FF')
    plt.plot([i for i in range(len(mean))], mean, '#CC0000')
    plt.legend(["truth", "Viterbi", "decoding", "mean"], loc=9)
    plt.title("Initial Decodings, %s" %suffix)
    plt.xlabel("locus")
    plt.ylabel("TMRCA")
    plt.savefig("output/plot_initial_%s.png" %suffix)


def main():

    # parse command line args
    (fasta_file,
    params_file,
    true_tmrca,
    train_iters,
    out_dir,
    suffix) = parse_args()

    # parse truth
    with open(true_tmrca, 'r') as f:
        truth = f.readlines()
        truth = list(map(lambda x: float(x.strip()), truth))


    # fasta_file = "example/sequences_4mu.fasta"
    # load observed data file
    observed = load_observed_data(fasta_file)

    # load paramaters file
    # p_init, p_trans, p_emit = load_params("example/initial_parameters_4mu.txt")
    p_init, p_trans, p_emit = parse_params(params_file)

    # pass observed data and params to new Viterbi object
    viterbi = ViterbiClass(observed, p_init, p_trans, p_emit)

    # calculate highest probability path
    path = viterbi.compute_path()

    # pass observed data and params to new ForwardBackward object
    fb = ForwardBackwardClass(observed, p_init, p_trans, p_emit)
    fb.compute_fb()
    decodings,post_means = fb.get_posteriors()

    # print(len(path), len(decodings), len(post_means))
    step_triples = list(zip(path,decodings,post_means))

    outfile = out_dir+"decodings_initial_"+suffix+".txt"
    write_decoding_to_file(step_triples, outfile)
    make_graph(truth, path, decodings, post_means, suffix)

    # test our data with expected data
    # with open(outfile, 'r') as f:
    #     our_data = f.readlines()
    #
    # with open("example/decodings_initial_4mu.txt", 'r') as f:
    #     expected = f.readlines()
    #
    # for i in range(len(our_data)):
    #     assert our_data[i] == expected[i]




if __name__ == "__main__":
  main()
