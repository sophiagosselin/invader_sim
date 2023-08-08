# INVADER SIM

Simulates the evolution of intein invaded extein sequences based on simulated phylogenies and using a Monte-Carlo like process.
Intein Invasion Sequence Simulator v1.0.0\n

Dependencies:


    perl (5.30.1)
      -BioPerl (1.7.8)
      -Bio::TreeIO

    R (4.3.1)
      -ape (5.7)

    paml (4.1.6)



This program will simulate extein and intein sequences according to user parameters.
Then using a monte carlo based approach, invade the extein sequences with the intein sequences.
Returns simulated samples of extein nucleotide sequences invaded by intein sequences.

## Usage Example:

    perl invader_sim.pl

The program can be run in naiive or prepared mode.
If the user wishes to be guided through the simulation process with prompts, then simply launch the program as above.


Else, a parameter file (is.param) can be provided using the following parameters:

These parameters apply to phylogeny simulations

	-sp -> number of species.
	-b -> birth rate.
	-d -> death rate.
	-sf -> sample fraction.
	-m -> mutation rate.\n


These parameters apply to sequence simulations

	-nn -> number of nucleotides.
	-a -> alpha parameter for gamma distribution
	-cat -> number of rate categories for the gamma distribution
	-m -> model number: (0: Poisson, 1: Proportional)


These parameters are general to the code

	-sn -> number of simulations to create
	-ws -> window size for MC chain (how many tips away can an intein jump)
	-evo -> string to call paml-evolver (likely paml-evolver, or evolver)

Note that the is.param must be formated as:
extein phylogeny simulation parameters
intein phylogeny simulation parameters
global params
extein sequence simulation parameters
intein sequence simulation parameters

An example is.param file:

	-sp 50 -b 1 -d 1 -sf .1 -m .01
	-sp 10 -b 1 -d 1 -sf .1 -m .1
	-sn 100 -ws 10 -evo evolver
	-nn 100 -a .5 -cat 4 -m 0
	-nn 50 -a 1 -cat 4 -m 1
