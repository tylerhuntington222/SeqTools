"""
ForwardBackwardClass.py

Definition of the Viterbi class for building a Hidden Markov
Model for determining TMRCA from DNA sequence data.
"""
from collections import defaultdict
from math import log, e

class ForwardBackwardClass:
    def __init__(self, observations, p_init, p_trans, p_emit):

        self.x = [int(i) for i in observations]
        self.p_inits = p_init
        self.p_trans = p_trans
        self.p_emits = p_emit
        self.states = list(sorted(self.p_trans.keys()))
        self.decodings = {0:0.32, 1:1.75, 2:4.54, 3:9.40}
        self.dp_table = defaultdict(dict)

    def compute_forward(self):


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
                q_hat = log_sum

                f_k_i = emit_term + log_sum

                # put this forward value into DP table
                self.dp_table[i][k] = \
                    {"forward": f_k_i, "backward": float("-inf")}

        print(self.dp_table[len(self.x)-1])
        last = list(self.dp_table[len(self.x)-1].values())
        vals = list(map(lambda x: x["forward"], last))
        print(vals)
        res = self.log_summer(self.log_summer(vals[0], vals[1]), \
        self.log_summer(vals[2], vals[3]))

        print(self.array_log_summer(vals))

        print(res)

    @staticmethod
    def log_summer(p, q):
        return(p + log(1 + e**(q-p)))


    @staticmethod
    def rec_array_log_summer(array, acc):
        if len(array) == 1:
            return ForwardBackwardClass.log_summer(array[0], acc)
        else:
            return ForwardBackwardClass.rec_array_log_summer(array[1:],\
             ForwardBackwardClass.log_summer(array[0], acc))

    @staticmethod
    def array_log_summer(array):
        return(ForwardBackwardClass.rec_array_log_summer(array[1:], array[0]))

    def compute_backward(self):
        pass

if __name__ == "__main__":

    fb = ForwardBackward
