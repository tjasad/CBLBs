# CBLBs: Configurable (Bio)Logic Blocks

This repository is supplementing the paper Configurable (bio)logic blocks for distributed biological computing [1]. 

## Folders
The project is organized in the following folders
* [```clb```](/clb/): deterministic and stochastic models of a CLBL and its modules.
* [```robustness```](/robustness/): analysis of robustness.
* [```latch```](/latch/): model of a mutlicellular latch, i.e. toggle switch, as described in [2].
* [```mux```](/muc/): model of a mutlicellular multiplexer and derivation of its kinetic parameters as described in [3].


## Main files
The main files that can be used to reproduce the results reported in the paper are as follows
* [```clb/models.py```](blob/master/clb/models.py): the implementation of deterministic and stochastic models

## References:

[1] Moškon M, et al. Configurable (bio)logic blocks for distributed biological computing. In review.

[2] Urrios A, Macia J, Manzoni R, Conde N, Bonforti A, de Nadal E, Posas F, Solé R. A synthetic multicellular memory device. ACS synthetic biology. 2016 Aug 19;5(8):862-73.

[3] Macia J, Manzoni R, Conde N, Urrios A, de Nadal E, Solé R, Posas F. Implementation of complex biological logic circuits using spatially distributed multicellular consortia. PLoS computational biology. 2016 Feb 1;12(2):e1004685.