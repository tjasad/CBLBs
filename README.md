# CBLBs: Configurable (Bio)Logic Blocks

This repository is supplementing the paper Configurable (bio)logic blocks for distributed biological computing [1]. 

## Folders
The project is organized in the following folders
* [```clb```](/clb/): deterministic and stochastic models of a CLBL and its modules.
* [```robustness```](/robustness/): analysis of robustness.
* [```latch```](/latch/): model of a mutlicellular latch, i.e. toggle switch, as described in [2].
* [```mux```](/mux/): model of a mutlicellular multiplexer and derivation of its kinetic parameters as described in [3].


## Main files
The main files that can be used to reproduce the results reported in the paper [1] are as follows
* [```clb/models.py```](/clb/models.py): the implementation of deterministic and stochastic models.
* [```clb/run_ode_model_clb.py```](blob/master/clb/run_ode_model_clb.py): deterministic simulation of a CBLB.
* [```clb/run_ssa_model_clb.py```](blob/master/clb/run_ssa_model_clb.py): stochastic simulation of a CBLB.
* [```robustness/run_solver.py```](blob/master/robustness/run_solver.py): perform the exploration of viable parameter space.
* [```robustness/robustness_analysis.py```](blob/master/robustness/robustness_analysis.py): visualize and analyse the results of the obtained viable parameter space.

## References:

[1] Moškon M, et al. Configurable (bio)logic blocks for distributed biological computing. In review.

[2] Urrios A, Macia J, Manzoni R, Conde N, Bonforti A, de Nadal E, Posas F, Solé R. A synthetic multicellular memory device. ACS synthetic biology. 2016 Aug 19;5(8):862-73.

[3] Macia J, Manzoni R, Conde N, Urrios A, de Nadal E, Solé R, Posas F. Implementation of complex biological logic circuits using spatially distributed multicellular consortia. PLoS computational biology. 2016 Feb 1;12(2):e1004685.