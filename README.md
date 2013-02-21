# Hello
Interested in modeling continuous trait evolution on a phylogeny as a Lévy process? This git repository contains the most recent source code for method we described in Systematic Biology, found <a href="http://sysbio.oxfordjournals.org/content/62/2/193.full">here</a>. We whipped together this guide to respond to questions received thus far and plan to release an official manual once we've collected feedback from users. Please feel welcome to email us with questions about the method or the model.

Our best,

Michael, Josh, and Mason

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

Number of MCMC cycles for the analysis (default 2000000)

*printFreqMH*

MCMC parameter sample frequency (default: 1000)

*printFreqJump*

MCMC per-branch jump size sample frequency, indexed by branch's immediately descendant node (default: 1000)

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
