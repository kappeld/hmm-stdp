# Getting started

This folder contains the souce code to reproduce the experimental results reported in (Kappel et al. 2014).

The experiments were run in Matlab R2007a (7.4.0.287) on a Linux machine (Fedora, Core release 6
(Zod)).

To get started run the
[run_experiment.m](https://github.com/kappeld/hmm-stdp/blob/master/run_experiment.m) script (you need writing permission on the current working
directory to seccesfully run). This will perform a whole training session and plot the results.
This experimental setup is equal to experiment 1 of (Kappel et al. 2014), to create the reults in Figure 2.
The results are saved to the `results/` subfolder of your current working directory that is created while
running the script.

# Code snoopies

The network dynamics and learning rules are implemented in the file [snn_sample_ctrf.m
](https://github.com/kappeld/hmm-stdp/blob/master/hmmsem/snn_sample_ctrf.m).
All updates are updated using continous time updates using a variant of the Gillespie algorithm.
For further information see (Kappel et al. 2014).

### References

D. Kappel, B. Nessler, and W. Maass. *STDP Installs in Winner-Take-All Circuits an Online
   Approximation to Hidden Markov Model Learning*. PLOS Computational Biology, 2014.
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003511
