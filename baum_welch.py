"""
baum_welch.py

Main file for implementing our ForwardBackwardClass. Given observations,
initial decodings, transition values, and emission values, this class can
create a Forward-Backward table, which is a table of dictionaries for all
forward and backward probabilities. This class can also compute posterior,
decoding, and mean values.

Author: Tyler Huntington
Author: Nathan Holeman

Date: 4/16/18
"""
from collections import defaultdict
from math import log, e

class ForwardBackwardClass:
    def __init__(self, observations, p_init, p_trans, p_emit):
        """
        @param observations: observed sequence
        @param p_init: initial decodings
        @param p_trans: transition values
        @param p_emit: emission values
        @return: None

        side-effects: initialize observations, inits, transitions, emissions,
                      state-mapping, a dynamic programming table for forward/
                      backward probabilities, log likelihood, posteriors,
                      decodings, and means.
        """

        self.x = [int(i) for i in observations]
        self.p_inits = p_init
        self.p_trans = p_trans
        self.p_emits = p_emit
        self.states = list(sorted(self.p_trans.keys()))
        self.state_dict = {0:0.32, 1:1.75, 2:4.54, 3:9.40}
        self.dp_table = defaultdict(dict)
        self.p_xbar = -1
        self.posteriors = []
        self.post_decodings = []
        self.post_means = []

    def get_posteriors(self):
        """
        @return: tuple of posterior decodings and means
        """
        return self.post_decodings, self.post_means

    def get_fw_bw_table(self):
        """
        @return: forward-backward table
        """
        return self.dp_table

    def compute_fb(self):
        """
        @return: None

        side-effects: method from which to run Baum-Welch algorithm: compute
                      forward probabilities, then backward probabiliites, then
                      posterior probabilities, then decodings, then means.
        """
        self.compute_forward()
        self.compute_backward()
        self.compute_posterior_probs()
        self.compute_posterior_decoding()
        self.compute_posterior_means()

    def compute_forward(self):
        """
        @return: None

        side-effects: Compute forward probabilities and update forward-backward
                      table.
        """

        # initialize first column of dp_table
        for i in range(len(self.p_inits)):
            emit_term = log(self.p_emits[i][self.x[0]])
            init_term = log(self.p_inits[i])
            f_k = emit_term + init_term
            self.dp_table[0][i] = {"forward": f_k, "backward": float("-inf")}


        # populate DP table
        for i in range(1, len(self.x)):
            for k in self.states:
                emit_term = log(self.p_emits[k][self.x[i]])
                log_sum = 0
                counter = 0
                for s in self.states:
                    f_l_prev = self.dp_table[i-1][s]["forward"]
                    prev = f_l_prev + log(self.p_trans[s][k])
                    if counter > 0:
                        log_sum = self.log_summer(log_sum, prev)
                    else:
                        log_sum = prev
                    counter += 1

                f_k_i = emit_term + log_sum

                # put this forward value into DP table
                self.dp_table[i][k] = \
                    {"forward": f_k_i, "backward": float("-inf")}

        last = list(self.dp_table[len(self.x)-1].values())
        vals = list(map(lambda x: x["forward"], last))

        self.p_xbar = self.array_log_summer(vals)

    @staticmethod
    def log_summer(p, q):
        """
        @param p: float of log value
        @param q: fload of log value
        @return: float of the log-sum of two log values
        """
        return(p + log(1 + e**(q-p)))


    @staticmethod
    def rec_array_log_summer(array, acc):\
        """
        @param array: list of float values (all log values)
        @param acc: float value to which all elements of the array are accumulated
        @return: accumulation of all values in array log-added together
        """
        if len(array) == 1:
            return ForwardBackwardClass.log_summer(array[0], acc)
        else:
            return ForwardBackwardClass.rec_array_log_summer(array[1:],\
             ForwardBackwardClass.log_summer(array[0], acc))

    @staticmethod
    def array_log_summer(array):
        """
        @param array: list of float values
        @return: return value of rec_array_log_summer -- accumulation of all values
                 in array log-added together.
        """
        if len(array) == 0:
            raise ValueError
        else:
            return(ForwardBackwardClass.rec_array_log_summer(array[1:], array[0]))

    def compute_backward(self):
        """
        @return: None

        side-effects: update forward-backward table with backward probabilities
        """
        for i in range(len(self.p_inits)):
            self.dp_table[len(self.dp_table)-1][i]["backward"] = log(1)

        for i in range(len(self.x)-2, -1, -1):
            for k in self.states:
                log_sum = 0
                counter = 0
                for s in self.states:
                    emit_term = log(self.p_emits[s][self.x[i+1]])
                    b_l_prev = self.dp_table[i+1][s]["backward"]
                    prev = b_l_prev + log(self.p_trans[k][s]) + emit_term
                    if counter > 0:
                        log_sum = self.log_summer(log_sum, prev)
                    else:
                        log_sum = prev
                    counter += 1

                b_k_i = log_sum
                self.dp_table[i][k]["backward"] = b_k_i

    def compute_posterior_probs(self):
        """
        @return: None

        side-effects: compute posterior probabilities and update posteriors
        """
        for i in range(len(self.x)):
            post_col = []
            for k in range(len(self.states)):
                f = self.dp_table[i][k]["forward"]
                b = self.dp_table[i][k]["backward"]
                p_i_k = ((f + b) - self.p_xbar, k)
                post_col.append(p_i_k)
            self.posteriors.append(post_col)

    def compute_posterior_decoding(self):
        """
        @return: None

        side-effects: compute posterior decodings and update post_decodings
        """
        dec =(list(map(lambda x:(max(x,key=lambda y:y[0]))[1],self.posteriors)))
        self.post_decodings =  list(map(lambda x: self.state_dict[x], dec))

    def compute_posterior_means(self):
        """
        @return: None

        side-effects: compute means and update post_means
        """
        for i in range(len(self.x)):
            mean = 0
            for k in range(len(self.states)):
                non_log_p = e**(self.posteriors[i][k][0])
                p_i_k = non_log_p * self.state_dict[k]
                mean += p_i_k
            self.post_means.append(mean)


if __name__ == "__main__":

    fb = ForwardBackward
