"""
hmm.py

A program for estimating the hidden state sequence of DNA sequence data
using Hidden Markov Models. Both Viterbi's Algorithm and the Forward Backward
algorithm are implemented.

Authors:
Nathan Holeman
Tyler Huntington

April 16, 2018
"""

import optparse
from math import log, e
from viterbi import ViterbiClass
from baum_welch import ForwardBackwardClass
import sys
import matplotlib.pyplot as plt
from collections import defaultdict


def parse_args():
    """
    Parses the command line arguments for running the program

    Params:
        None
    Returns:
        Tuple containing the following items:

        observed_data_fasta_file
        input_parameter_file
        true_tmrca_vals_file
        num_training_iterations
        output_directory
        output_file_suffix
    """

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
    states = list(map(lambda x: int(x), p_init.keys()))
    return (p_init, p_trans, p_emit, states)

def load_observed_data(fasta_file):
    """
    Loads the observed data in .fasta file format for use in Hidden Markov
    Model.

    Params:
        fasta_file: filepath to fasta file containing two sequences of
            from different individuals over the same region of the genome.

    Returns: string of "0"s and "1"s in which each numeral corresponds to
        a match ("0") or mismatch ("1") between the two input sequences
        at a particular site
    """

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
    """
    Writes the Viterbi decoding, posterior decoding and posterior mean
    of the hidden sequence to an output file.

    Params:
        step_data: a tuple of Viterbi decodings, posterior decodings and
        posterior means structured as lists.

        outfile: file to which the decodings should be written

    Returns:
        None

    Side Effects:
        Writes the three types of input decodings to the output file
        passed by caller.
    """

    with open(outfile, 'w') as f:
        f.write("# Viterbi_decoding posterior_decoding posterior_mean")
        f.write("\n")
        for vit_dec, post_dec, post_mean in step_data:
            line = " ".join([str(vit_dec), str(post_dec), \
            "%.6f" % (round(post_mean, 6))])
            f.write(line)
            f.write("\n")

def make_graph(truth,viterbi,decoding,mean,suffix,tag):
    """
    Creates a plot of the true state sequence compared to the
    estimated state sequence determined by Viterbi decoding, posterior
    decoding and posterior mean calculations.

    Params:
        truth: a list of the true TMRCA vals (i.e. the hidden sequence)
        viterbi: list of the viterbi decodings
        decodings: listo of the posterior decoding vals
        mean: list of the posterior means
        suffix: a suffix to be attached to the output graph file
        tag: "initial" or "estimated" indicating the nature of the parameters
            used to generate the estimated sequences.
    Returns:
        None
    Side Effects:
        Creates plot and exports it to /ouput directory.
    """

    plt.plot([i for i in range(len(truth))], truth, '#000000', linewidth=2)
    plt.plot([i for i in range(len(viterbi))], viterbi, '#006600')
    plt.plot([i for i in range(len(decoding))], decoding, '#0000FF')
    plt.plot([i for i in range(len(mean))], mean, '#CC0000')
    plt.legend(["truth", "Viterbi", "decoding", "mean"], loc=9)
    plt.title("%s decodings, %s" %(tag, suffix))
    plt.xlabel("locus")
    plt.ylabel("TMRCA")
    plt.savefig("output/plot_%s_%s.png" %(tag,suffix))
    plt.clf()


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
    xs = list(map(int, observed))

    # load paramaters file
    # p_init, p_trans, p_emit = load_params("example/initial_parameters_4mu.txt")
    p_init, p_trans, p_emit, states = parse_params(params_file)

    # pass observed data and params to new Viterbi object
    viterbi = ViterbiClass(observed, p_init, p_trans, p_emit)

    # calculate highest probability path
    path = viterbi.compute_path()

    # pass observed data and params to new ForwardBackward object
    fb = ForwardBackwardClass(observed, p_init, p_trans, p_emit)
    fb.compute_fb()
    fw_bw_table = fb.get_fw_bw_table()
    decodings, post_means = fb.get_posteriors()

    step_triples = list(zip(path,decodings,post_means))

    outfile = out_dir+"decodings_initial_"+suffix+".txt"
    write_decoding_to_file(step_triples, outfile)
    make_graph(truth, path, decodings, post_means, suffix, "initial")


    # Part 2: Baum-Welch Training
    print("Baum-Welch")

    # initialize values for first training iteration
    curr_fb = fb
    curr_fw_bw_table = fw_bw_table
    est_p_init = p_init
    est_p_trans = p_trans
    est_p_emit = p_emit

    likelihoods = [curr_fb.p_xbar]

    for step in range(train_iters):
        print("Iteration: %d" %step)
        print(curr_fb.p_xbar)
        est_p_init = estimate_p_init(curr_fw_bw_table, curr_fb, states)
        est_p_trans = estimate_p_trans(curr_fw_bw_table,curr_fb,states, \
                                        xs,est_p_trans,est_p_emit)
        est_p_emit = estimate_p_emit(curr_fw_bw_table, curr_fb, states, xs)
        curr_fb = ForwardBackwardClass(observed, est_p_init, est_p_trans, \
                                        est_p_emit)
        curr_fb.compute_fb()
        curr_fw_bw_table = curr_fb.get_fw_bw_table()

    likelihoods.append(curr_fb.p_xbar)
    write_log_liklihoods("output/likelihoods_%s.txt" %suffix,likelihoods)

    params = (est_p_init, est_p_trans, est_p_emit)
    write_estimated_params_file("output/estimated_parameters_%s.txt" %suffix, params)


    # pass observed data and params to new Viterbi object
    viterbi = ViterbiClass(observed, est_p_init, est_p_trans, est_p_emit)

    # calculate highest probability path
    path = viterbi.compute_path()

    # pass observed data and params to new ForwardBackward object
    fb = ForwardBackwardClass(observed, est_p_init, est_p_trans, est_p_emit)
    fb.compute_fb()
    fw_bw_table = fb.get_fw_bw_table()
    decodings, post_means = fb.get_posteriors()

    step_triples = list(zip(path,decodings,post_means))

    outfile = "output/decodings_estimated_"+suffix+".txt"
    write_decoding_to_file(step_triples, outfile)
    make_graph(truth, path, decodings, post_means, suffix, "estimated")

def write_log_liklihoods(outfile, likelihoods):
    """
    Writes the log likelihoods to an output file.

    Params:
        outfile: output file name
        liklihoods: list of two elements where the first element is the
            initial log liklihood of the estimated sequence given the
            observed data and the second in the same probability
            calculated when using estimated paramters following 16 iterations
            of Baum Welch training
        returns:
            None
        Side-Effects:
            Writes liklihoods to output file
    """

    with open(outfile, 'w') as f:
        line_1 = "# Likelihood under {initial, estimated} parameters"
        f.write(line_1)
        f.write("\n"+str(likelihoods[0]))
        f.write("\n"+str(likelihoods[1]))

def write_estimated_params_file(outfile, params):
    """
    Writes the estimated params to an output file.

    Params:
        outfile: output file name
        params: tuple of:
            estimated initial probabilities
            estimated transition probabilities
            estimated emission probabilities

        Returns:
            None
        Side-Effects:
            Writes estimated params to output file
    """
    string = display_params(params)
    with open(outfile, 'w') as f:
        f.write(string)



def estimate_p_init(fw_bw_table, fb, states):
    """
    Estimates initial probabilities based on a run of the forward backward
    algorithm.

    Params:
        fw_bw_table: the DP table generated by running the forward backward
            algorith.
        fb: an object of the ForwardBackwardClass from which estimated params
            should be determined
        states: a list of possible states
    Returns:
        estimated initial probabilities for each state.

    """

    # calculate pi
    pi_vals = {}
    for k in states:
        f_val = fw_bw_table[0][k]["forward"]
        b_val = fw_bw_table[0][k]["backward"]
        p_val = fb.p_xbar
        pi_vals[k] = e**((f_val + b_val) - p_val)
    return pi_vals

def estimate_p_trans(fw_bw_table, fb, states, xs, p_trans, p_emit):
    """
    Estimates transition probabilities based on a run of the forward backward
    algorithm.

    Params:
        fw_bw_table: the DP table generated by running the forward backward
            algorithm.
        fb: an object of the ForwardBackwardClass from which estimated params
            should be determined
        states: a list of possible states
        xs: the observed state sequence as a list of ints
        p_trans: the initial transition probabilities
        p_emit: the inital emission probabilities

    Returns:
        estimated transition probabilities for each state.

    """
    A_kl_table = []
    # for i in range(train_iters):
    for k in states:
        A_kl_col = []
        for l in states:
            A_kl_cur = 0
            for step in range(0, len(xs)-1):
                # calculate P -- returned as non-log value
                P = calc_expected_transition(fb, k, l, step, xs, \
                                             fw_bw_table, p_trans, p_emit)
                # increment A_kl_cur accumulator
                A_kl_cur += P
            A_kl_col.append(A_kl_cur)
        # append this column of A_kl vals to table
        A_kl_table.append(A_kl_col)

    # estimated transition probs normalization step
    a_kl_table = defaultdict(dict)
    for k in states:
        for l in states:
            acc = 0
            for l_prime in states:
                acc += A_kl_table[k][l_prime]

            a_kl_table[k][l] = A_kl_table[k][l] / acc
    return a_kl_table

def estimate_p_emit(fw_bw_table,fb,states,xs):
    """
     Estimates emission probabilities based on a run of the forward backward
    algorithm.

    Params:
        fw_bw_table: the DP table generated by running the forward backward
            algorithm.
        fb: an object of the ForwardBackwardClass from which estimated params
            should be determined
        states: a list of possible states
        xs: the observed state sequence as a list of ints

    Returns:
        estimated emission probabilities for each state.

    """
    # calculate emissions
    E_kb_table = defaultdict(dict)
    for k in states:
        for b in range(2):
            acc = ""
            for i in range(len(xs)):
                if xs[i] == b:
                    f_val = fw_bw_table[i][k]["forward"]
                    b_val = fw_bw_table[i][k]["backward"]
                    p_val = fb.p_xbar
                    res = (f_val + b_val) - p_val
                    if acc == "":
                        acc = res
                    else:
                        acc = ForwardBackwardClass.log_summer(acc, res)
            E_kb_table[k][b] = acc

    # emissions normalization step
    e_kb_table = defaultdict(dict)
    for k in states:
        for b in range(2):
            acc = ""
            for b_prime in range(2):
                res = E_kb_table[k][b_prime]
                if acc == "":
                    acc = res
                else:
                    acc = ForwardBackwardClass.log_summer(acc, res)
            e_kb_table[k][b] = e**(E_kb_table[k][b] - acc)
    return(e_kb_table)

def calc_expected_transition(fb, k, l, i, xs, fw_bw_table, p_trans, p_emit):
    """
    Calculates the expected transition term for estimating transition
        probality of a given state and step in the observed sequence.

    Params:
        fb: an object of the ForwardBackwardClass from which estimated params
            should be determined
        k: the state from which we are transitioning out of
        l: the state we are transitioning to
        xs: the observed state sequence as a list of ints

        fw_bw_table: the DP table generated by running the forward backward
            algorithm.
        fb: an object of the ForwardBackwardClass from which estimated params
            should be determined
        p_trans: the initial transition probabilities
        p_emit: the inital emission probabilities

    Returns:
        Expected transition term for a given k, l, and i in our observed data.

    """
    f_k_i = fw_bw_table[i][k]["forward"]
    a_k_l = log(p_trans[k][l])
    e_l_term = log(p_emit[l][xs[i+1]])
    b_l_term = fw_bw_table[i+1][l]["backward"]
    p_xbar = fb.p_xbar

    non_log_res = e**(((f_k_i + a_k_l + e_l_term + b_l_term) - p_xbar))
    return(non_log_res)






























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
