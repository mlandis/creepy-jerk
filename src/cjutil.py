import os
import sys
import time
import numpy
from numpy import *

def ssExec(	execFp = '/Users/mlandis/Documents/code/creepy-jerk/',
		execCmd = './creepy-jerk',
		outFp = "/Users/mlandis/data/expr_phylo/output/",
		simName = "",
		timeName = True,
		numBeta = 20,
		argStr = ""
	):
	# create output directory
	if simName != "":
		simName += "."
	simName += time.strftime("%y%m%d%H%M%S")
	simName += ".numBeta" + str(numBeta)
	simFp = outFp + simName + "/"
	simDir = os.path.dirname(simFp)
	if not os.path.exists(simDir):
		os.makedirs(simDir)
	# generate and store beta values
	alpha = 0.3
	beta = []
	for k in range(0,numBeta+1):
		beta.append(pow(float(k)/numBeta, 1/alpha))
	f = open(simFp + "beta.txt", 'w')
	for b in beta:
		f.write(str(b) + "\t")
	f.write("\n")
	f.close()
	# construct argStr and save settings
	printStdOutStr = " -printStdOut=True"
	argStr +=  ' -printFreqMH=1000'
	argStr += ' -numCycles=1000000'
	argStr += ' -printFreqJump=1000001'
	argStr += ' -outputDirPath=' + simFp
	#argStr += ' -seed=6'
	argStr += ' -modelType=3'
	argStr += ' -useJumpKernel=True'
	argStr += ' -tuningBm=0.9995'
	args = [a for a in argStr.split(" ") if a != ""]
	f = open(simFp + "args.txt", 'w')
	for a in args:
		f.write(str(a) + "\n")
	f.close()
	# execute MCMC under beta values
	for k in range(0,numBeta):
		if k == numBeta:
			printStdOutStr = " -printStdOut=True"
		kStr = str(k).zfill(len(str(numBeta-1))) 
		os.system(execFp + execCmd + argStr + " -simName=" + simName + "." + kStr + " -useSteppingStone=True -betaSteppingStone=" + str(beta[k]) + printStdOutStr + " &")

def ssBatch(execFp = '/Users/mlandis/Documents/code/creepy-jerk/', execCmd = './creepy-jerk', outFp = "/Users/mlandis/data/expr_phylo/output/", simName = "", timeName = True, numBetaVec = [2,3,4,5,6,8,10,12,14,16,20]):
	for b in numBetaVec:
		ssExec(execFp=execFp, execCmd=execCmd, outFp=outFp, simName=simName, timeName=timeName, numBeta=b, argStr = " -seed=" + str(b))

def cjBatch(	execFp = "/Users/mlandis/Documents/code/creepy-jerk/",
		execCmd = "./creepy-jerk",
		argsFp = "/Users/mlandis/data/expr_phylo/input/",
		defaultFn = "default.settings",
		batchFn = "batch.settings",
		outFp = "/Users/mlandis/data/expr_phylo/output/",
		timeName = True,
		argStr = ""
	):
	# read in default settings
	defaultFile = open(argsFp + defaultFn, "r")
	for line in defaultFile:
		argStr += line.rstrip("\n\r")
	print("Default:\t" + argsFp + defaultFn + "\t" + argStr + "\n")
	defaultFile.close()	
	# read in batch file and execute
	batchFile = open(argsFp + batchFn, "r")
	for line in batchFile:
		tempArgStr = argStr
		tempArgStr += " " + line.rstrip("\n\r")
		print("Job:\t" + argsFp + batchFn + "\t" + tempArgStr + "\n")
		os.system(execFp + execCmd + " " + tempArgStr + " &")
