"""
ForwardBackwardClass.py

Definition of the Viterbi class for building a Hidden Markov
Model for determining TMRCA from DNA sequence data.
"""
from collections import defaultdict
from math import log

class ForwardBackwardClass:
    def __init__(self, observations, p_init, p_trans, p_emit):

        self.x = [int(i) for i in observations]
        self.p_inits = p_init
        self.p_trans = p_trans
        self.p_emits = p_emit
        self.states = list(sorted(self.p_trans.keys()))
        self.decodings = {0:0.32, 1:1.75, 2:4.54, 3:9.40}
        self.best_path = []

    def compute_forward(self):
        pass


    def compute_backward(self):
        pass
