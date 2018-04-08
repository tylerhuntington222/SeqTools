"""
ViterbiClass.py

Definition of the Viterbi class for building a Hidden Markov
Model for determining TMRCA from DNA sequence data.
"""

class Viberbi:

    def __init__(p_init, p_trans, p_emit):

        self.p_inits = p_init
        self.p_trans = p_trans
        self.p_emits = p_emit


    def compute_path(self):
