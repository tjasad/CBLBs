# CBLBs: Configurable (Bio)Logic Blocks

This repository is supplementing the paper Field-programmable biological circuits based on the configurable (bio)logic blocks for distributed biological computin [1]. 

## Folders
The project is organized in the following folders
* [```clb```](/clb/): deterministic and stochastic models of a CLBL and its modules.
* [```robustness```](/robustness/): analysis of robustness.
* [```latch```](/latch/): model of a multicellular latch, i.e. toggle switch, as described in [2].
* [```mux```](/mux/): model of a multicellular multiplexer and derivation of its kinetic parameters using the transfer-function models and values from [3].


## Main files
The main files that can be used to reproduce the results reported in the paper [1] are as follows
* [```clb/models.py```](/clb/models.py): the implementation of deterministic and stochastic models.
* [```clb/run_ode_model_clb.py```](/clb/run_ode_model_clb.py): deterministic simulation of a CBLB.
* [```clb/run_ssa_model_clb.py```](/clb/run_ssa_model_clb.py): stochastic simulation of a CBLB.
* [```robustness/run_solver.py```](/robustness/run_solver.py): perform the exploration of viable parameter space as described in [4].
* [```robustness/robustness_analysis.py```](/robustness/robustness_analysis.py): visualize and analyse the results of the obtained viable parameter space.

## Requirements
Before running the code the following libraries should be installed
* `DEAP`
* `multiprocessing`
* `SciPy`
* `NumPy`
* `seaborn`
* `matplotlib`
* `pandas`
* `pickle`

## References:

[1] Moškon M, et al. Field-programmable biological circuits based on configurable (bio)logic blocks for distributed biological computin. In preparation.

[2] Urrios A, Macia J, Manzoni R, Conde N, Bonforti A, de Nadal E, Posas F, Solé R. A synthetic multicellular memory device. ACS synthetic biology. 2016 Aug 19;5(8):862-73.

[3] Macia J, Manzoni R, Conde N, Urrios A, de Nadal E, Solé R, Posas F. Implementation of complex biological logic circuits using spatially distributed multicellular consortia. PLoS computational biology. 2016 Feb 1;12(2):e1004685.

[4] Pušnik Ž, Mraz M, Zimic N, Moškon M. Computational analysis of viable parameter regions in models of synthetic biological systems. Journal of biological engineering. 2019 Dec 1;13(1):75.