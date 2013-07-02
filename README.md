# Hello
Interested in modeling continuous trait evolution on a phylogeny as a Lévy process? This git repository contains the most recent source code for the method we described in Systematic Biology, found <a href="http://sysbio.oxfordjournals.org/content/62/2/193.full">here</a>. We whipped together this guide to respond to questions received thus far. Please feel welcome to email us with questions about the method or the model: mlandis (at) berkeley (dot) edu.

Our best,

Michael, Josh, and Mason

# Update
* July 2, 2013
** Method now produces TreeFig compatible .jump_snr.txt output file to show which branches have strong signal for jumps (described in paper).
** SNR is now scaled by inverse square root of the branch length instead of just the branch length.
** Improved jump normal (modelType=3) performance.


# Installation
1. From the command line, navigate into the <code>src</code> folder.
2. Before compiling, verify you have the GNU Scientific Library (GSL) installed. For Mac OS X, you can install GSL using <a href="http://mxcl.github.com/homebrew/">Homebrew</a<> as follows: <code>brew install gsl</code>. For Ubuntu, you can install GSL as follows: <code>apt-get install libgsl0-dev</code>. Advice for how to accomplish on Windows is welcome.
3. To compile, issue the command: <code>g++ -O3 -lgsl *.cpp -o creepy-jerk</code>
4. You should now find an executable called <code>creepy-jerk</code> in your current directory.

# Examples

- Copy the <code>creepy-jerk</code> executable into th <code>examples</code> folder. Here, you will find the tree file, the taxa file, and the three data files used in Landis, Schraiber, and Liang (2013).
- To run an example, execute <code>./primates.mass.sh</code>, <code>./primates.ecv.sh</code>, or <code>./primates.mass_ecv_ratio.sh</code>.
- You may find it easier to simply modify these text files than to construct new command lines (next section). NOTE: all examples are configured to run with modelType=3, the compound Poisson process with normally distributed jumps.

# Going wild
- The command-line syntax is as follows: <code>./creepy-jerk -flag_name1=flag_value1 -flag_name2=flag_value2 ... </code>. If a flag is not called, the application runs under the default values given below.
- Example command-line string: ```./creepy-jerk -modelType=3 -simName=example -printStdOut=True -printFreqStdOut=1000 -printFreqJump=1000 -printFreqMH=1000 -numCycles=2000000 -outputFilePath=./ -inputFilePath=./ -dataFileName=data.txt -taxaFileName=taxa.txt -treeFileName=tree.txt</code>```

####MCMC/Model settings

*modelType*

Analysis model type: 1=alpha-stable, 3=compound Poisson with normally distributed jumps, 4=variance Gamma, 5=pure Brownian motion (default 3)

*numCycles*

Number of MCMC cycles for the analysis (default 10000000)

*printFreqMH*

MCMC parameter sample frequency (default: 10000)

*printFreqJump*

MCMC per-branch jump size sample frequency, indexed by branch's immediately descendant node (default: 10000)

*snrBurnIn*

Burn-in period for computing the branch-length-normalized signal-to-noise ratio tree figure (default: 0.25 * numCycles)

*seed*

Random number generator seed (default: system time, i.e. effectively random)

*sigmaJumpProposal*

Standard deviation of the normally distributed jump proposal, adjustable in case of poor mixing (default: 1.0)

*tuningBM*

Briefly, this is the fraction of the branch length converted into Brownian motion, and is used to prohibit singularities in our formulation of the Levy pdf. Default setting recommended (default: 0.9995).

####Input/output settings

*treeFileName*

Name of file containing Newick formatted string with taxa labels and branch lengths (default: "")

*taxaFileName*

Name of file containing list taxa present in tree (default: "")

*dataFileName*

Name of file containing continuous traits, ordered according to taxaFileName (default: "")

*inputFilePath*

Full filepath of input files (default: "./")

*outputFilePath*

Full filepath of output files (default: "./")

*printStdOut*

Enable printing of MCMC parameter state to console (default: true)

*printFreqStdOut*

Console print frequency (default: 1000)

# Output files

####*output_filename*.parameters.txt

To interpret your MCMC posterior, we recommend using <a href="http://tree.bio.ed.ac.uk/software/tracer/">Tracer</a>. The parameters file contains the MCMC chain states sampled at the interval defined by *printFreqMH*. All fields are tab-delimited, with each column corresponding to a state variable, with those states being (in order):

*Cycle*

The MCMC sample's discrete time index

*lnL*

The full model log likelihood (i.e. lnL = lnB + lnJ, see below)

*lnB*

The pure-diffusion model log likelihood

*lnJ*

The pure-jump model log likelihood

*var*

The variance per unit time, which is a function of the *modelType* variable.

*kurt*

The kurtosis per unit time, which is a function of the *modelType* variable

*parameterName-modelType*

The model parameters, where sigma-modelType is the Brownian motion component rate parameter, and the remaining parameters (if any) are the jump component parameters



####*output_filename*.jumps.txt and *output_filename*.jump_snr.txt

The jumps file contains the MCMC chain states sampled at the interval defined by *printFreqJump*. Each row contains a Newick string with node annotations reporting the sampled sum of trait change drawn from the jump measure (hereafter, jump values). Although the jump values are associated with nodes, they correspond to the jumps occurring along the branch from the node's ancestor to the node itself. Note, the topology of the tree remains constant in the current implementation of creepy-jerk.

The jump_snr file contains the final branch-length-normalized signal-to-noise ratio of the tree's posterior jump values. The color scale is defined such that a value of 0.0 is gray, and the intensity of colors indicates the deviation in this value from 0.0. Finally, these values exclude the burn-in period as defined by the *snrBurnIn* flag (default: 0.25 * numCycles).

In the convention of the New Hampshire eXtended (.nhx) format, jump values for branches are given as [&j=value] immediately following the taxon label or divergence event.

To quickly view the posterior jumps.txt and jumps_snr.txt output, we recommend using the tree visualization software, <a href="http://tree.bio.ed.ac.uk/software/figtree/">FigTree 1.4</a>. Load the jump file by selecting File -> Open.

