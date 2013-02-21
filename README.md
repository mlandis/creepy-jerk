
g++ -O3 -lgsl *.cpp -o creepy-jerk
./creepy-jerk -printStdOut=True -printFreqStdOut=1000 -printFreqJumps=1000 -printFreqMH=1000 -numCycles=2000000 -useSteppingStone=False -tuningBM=0.9995 -outputFilePath=/User    s/mlandis/data/expr_phylo/output/ -inputFilePath=/Users/mlandis/data/expr_phylo/input/


Flags:

MCMC/Model settings
    -modelType       Integer indicating model type (1: alpha-stable, 3: compound Poisson with normally distributed jumps, 4: variance Gamma, 5: pure Brownian motion).
    -numCycles       Number of MCMC cycles
    -seed            Random number generator seed
    -tuningBM        (default: 0.9995)
    -useJumpKernel   (for testing: disables jumps)
    -sigmaJumpProposal (default: sigma=1.0)
    -printFreqMH     MCMC sample frequency of parameters
    -printFreqJump   MCMC sample frequency of jump sizes assigned to branches (leading to nodes)

Input/output
    -treeFileName     Name of file containing Newick formatted string with taxa labels and branch lengths
    -taxaFileName     Name of file containing list taxa present in tree
    -exprFileName    Name of file containing continuous triats, ordered according to taxaFileName
    -inputFilePath    Full filepath of input files (default local dir)
    -outputFilePath   Full filepath of output files (default local dir)
    -printStdOut      Enable printing of MCMC parameter state to console
    -printFreqStdOut  Frequency of printing to console
