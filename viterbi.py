"""
ViterbiClass.py

Definition of the Viterbi class for building a Hidden Markov
Model for determining TMRCA from DNA sequence data.
"""
from collections import defaultdict
from math import log

class ViterbiClass:
    def __init__(self, observations, p_init, p_trans, p_emit):

        self.x = [int(i) for i in observations]
        self.p_inits = p_init
        self.p_trans = p_trans
        self.p_emits = p_emit
        self.states = list(sorted(self.p_trans.keys()))
        self.decodings = {0:0.32, 1:1.75, 2:4.54, 3:9.40}
        self.best_path = []

    def compute_path(self):
        dp_table = defaultdict(dict)

        # initialize first column of dp_table
        for i in range(len(self.p_inits)):
            v_k = log(self.p_inits[i]) + log(self.p_emits[i][self.x[0]])
            dp_table[0][i] = {v_k: -1}

        # populate DP table
        for i in range(1, len(self.x)):
            for k in self.states:
                max_prev = float("-inf")
                for s in self.states:
                    v_l_prev = list(dp_table[i-1][s].keys())[0]
                    prev = v_l_prev + log(self.p_trans[s][k])
                    if prev > max_prev:
                        max_prev = prev
                        best_s = s
                v_k_i = log(self.p_emits[k][self.x[i]]) + max_prev
                dp_table[i][k] = {v_k_i: best_s}

        # initialize most probable path
        path = []

        last_col = dp_table[len(self.x)-1]
        # print(last_col)
        max_final_prob = float("-inf")
        for i in self.states:
            cur = list(last_col[i].keys())[0]
            # print(cur)
            if cur > max_final_prob:
                max_final_prob = cur
                best_final_state = i

        # print("Max final prob: ", max_final_prob)
        # print("Highes prob final state: ", best_final_state)

        path.append(best_final_state)


        for i in range(len(self.x)-2, -1, -1):
            backpointer = path[-1]

            prev_state = list(dp_table[i+1][backpointer].values())[0]
            path.append(prev_state)

        # reverse the path
        path = list(reversed(path))

        path = list(map(lambda x: self.decodings[x], path))
        self.best_path = path
        return(self.best_path)

    def get_best_path(self):
        return self.best_path









if __name__ == '__main__':
    assert 1==1
