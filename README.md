## CS 68 Lab 8 - Hidden Markov Models

Name 1: Nathan Holeman

Name 2: Tyler Huntington

userId1: nholema1

userId2: thuntin1

Number of Late Days Using for this lab: 0

---

### Analysis Questions

1. What general observations do you make about the differences/similarities between our 3 ways of estimating
the hidden state sequence (Viterbi, decoding, mean)? Which would you choose if you could only pick one?
One similarity we observed was that all three ways of estimating the hidden state
sequence perform better at higher mutation rates.
Both Viterbi and decoding exhibit a stepwise pattern, whereas the mean is more
curved and has many local variations.
Generally, decoding follows the truth closely. Viterbi seems to miss instances
of local variation in the truth; that is, Viterbi appears to settle on a state
and not transition out of that state, at the expense of deviating from the truth.
The mean roughly approximates the true sequence, but rarely finds the exact
values of the truth.

3. For the two true Tmrca sequences (`true_tmrca_test.txt` and `true_tmrca.txt`), how did the estimated
state sequences change between using the "initial" parameters and using the "estimated" parameters?
Why might we observe this difference? How did the log-likelihoods change for these two datasets?
The estimated state sequences became closer to the truth after multiple EM
iterations. We suspect that each round of Baum-Welch training tuned given parameters
to fit the model more closely to the training data. As a result, we observed an increase
in log likelihood of the estimated sequence from the initial parameter values to
the estimated values after 16 iterations of Baum-Welch.


4. For the `true_tmrca.txt` dataset, what was the difference in your output plots between mu, 2mu, and 5mu?
Why might we observe this difference as the mutation rate increases?
At low mutation rates, all three methods of estimating the state sequence where
less accurate. They appear to jointly increase in accuracy at higher mutation
rates. At higher mutation rates, there is more variation in the data, which can
make it easier for our algorithms to detect changes in the change sequence. For
a given locus, if there is a low mutation rate and a mutation is observed, the
algorithms are less likely to estimate a transition to another state.

---

### Lab Questionnaire

(None of your answers below will affect your grade; this is to help refine lab assignments in the future)

1. Approximately, how many hours did you take to complete this lab (provide your answer as a single integer on the line below).
16
2. How difficult did you find this lab?  (1-5, with 5 being very difficult and 1 being very easy)
4
3. Describe the biggest challenge you faced on this lab:
Working into and out of log space; debugging incorrect emission values.
