library(ape)
library(phylobase)

auto_cor =	function(	data,
				lag
			)
{
	indices = 1:(length(data)-lag)
	return(cor(data[indices],data[indices+lag]))
}

plotMcmc = 	function(	fn,
				b = 100,
				trueparm = c(),
				plot_auto = F,
				max_lag = 50,
				other_func = F,
				func_name = "",
				true_func,
				func,
				...
			)
{
	
	x = read.table(fn, header=T)
	xx = x[b:nrow(x),]
	dev.new()
	parm_names = colnames(xx)[5:ncol(xx)]
	num_plots = length(parm_names)
	
	if (other_func) 
	{
		num_plots = num_plots + 1
	}
	if (plot_auto)
	{
		par(mfrow=c(num_lots,4),oma=c(0,0,2,0))
	}
	else
	{
		par(mfrow=c(num_plots,3),oma=c(0,0,2,0))
	}
	
	for (i in 1:length(parm_names)) {
		plot(density(xx[,4+i], from=0),xlab=parm_names[i],main="Posterior")
		abline(v=trueparm[i])
		plot(xx[,4+i],xx$lnL,xlab=parm_names[i],ylab="lnL",main="Likelihood")
		plot(xx[,4+i],type="l",ylab=parm_names[i],main="Trace")
		if (plot_auto)
		{
			auto = sapply(1:max_lag,function(i) {auto_cor(xx[,4+i],i)})
			plot(auto,ylab="autocorrelation",main="Autocorrelation")	
		}
	}
	
	if (other_func)
	{
		func_vals = func(xx,...)
		plot(density(func_vals),xlab=func_name,main="Posterior")
		abline(v=true_func)
		plot(func_vals,xx$lnL,xlab=func_name,ylab="lnL",main="Likelihood")
		plot(func_vals,type="l",ylab=func_name,main="Trace")
		if (plot_auto)
		{
			auto = sapply(1:max_lag,function(i) {auto_cor(func_vals,i)})
			plot(auto,ylab="autocorrelation",main="Autocorrelation")
		}
	}
	
	
#	titleName = paste("SigBM=", trueparm[1], ", LamJN=", trueparm[2], ", SigJN=", trueparm[3], sep="")
	title(main=fn, outer=T)
	
#	dev.new()
#	plot(xx$lnL)
	
	return(xx)
}

plotJumps =	function(
				parmFp = "",
				jumpFp = "",
				treeFp = "~/data/expr_phylo/input/primates.eastman.isler_pruned.tree.txt",
				taxaFp = "~/data/expr_phylo/input/primates.eastman.isler_pruned.taxa.txt",
				b = 100,
				plot_auto = F,
				max_lag = 50
			)
{

	# read taxa
	taxa = scan(file=taxaFp,what="") 

	# read trees
	tRaw = reorder(phylo4(read.tree(treeFp)), "postorder")

	edgeLength(tRaw)[which(edgeLength(tRaw) == 0.0)] = 10^-4
	tJump = tRaw
	tRawJump = tRaw
	tVarJump = tRaw
	tAbsJump = tRaw
	tBoth = tRaw
	tBoth2 = tRaw

	# read in parameters
	dp = read.table(parmFp, header=T)
	dp = dp[b:nrow(dp),]

	# get jump statistics
	meanDrift = mean(dp$Sig.BM)
	varDrift = var(dp$Sig.BM)
	medDrift = median(dp$Sig.BM)

	# read in jumps
	dj = read.table(jumpFp, header=T)
	dj = dj[b:nrow(dj),]

	# get Newick string from dj.names()
	nodes = as.numeric(substr(names(dj)[-1],2,nchar(names(dj)[-1]))) + 1
	#newickStr = getNewickStr(nodes, nrow(tRaw$edge))
	
	# get jump statistics
	meanJumps = apply(dj[,2:ncol(dj)], 2, mean)
	meanJump2 <<- meanJumps
	varJumps = apply(dj[,2:ncol(dj)], 2, var)
	medJumps = apply(dj[,2:ncol(dj)], 2, median)
	meanTipJumps = apply(dj[,2:(length(tipLabels(tRaw))+1)], 2, mean)
	varTipJumps = apply(dj[,2:(length(tipLabels(tRaw))+1)], 2, var)
	medTipJumps = apply(dj[,2:(length(tipLabels(tRaw))+1)], 2, median)
	meanNodeJumps = apply(dj[,(length(tipLabels(tRaw))+1):ncol(dj)], 2, mean)
	varNodeJumps = apply(dj[,(length(tipLabels(tRaw))+1):ncol(dj)], 2, var)
	medNodeJumps = apply(dj[,(length(tipLabels(tRaw))+1):ncol(dj)], 2, median)
	
	# augment trees
	edgeLength(tAbsJump) = c(as.vector(abs(meanJumps)), NA) / edgeLength(tRaw)
	edgeLength(tRawJump) = c(as.vector(abs(meanJumps)), NA)
	edgeLength(tVarJump) = c(as.vector(abs(varJumps)), NA)
	edgeLength(tBoth) = edgeLength(tRaw) * meanDrift + c(abs(as.vector(meanJumps)), NA)
	edgeLength(tBoth2) = (edgeLength(tRaw) * meanDrift + c(abs(as.vector(meanJumps)), NA))/edgeLength(tRaw)

	# plot trees
	dev.new()
	par(mfrow=c(1,5), oma=c(0,0,2,0))
	dev.new();plot(tRaw, main="raw")
	#dev.new();plot(tAbsJump, main="abs(mean(jump))")
	#dev.new();plot(tRawJump, main="abs(mean(jump))")
	#plot(tVarJump, main="var(jump)")
	#dev.new();plot(tBoth, main="both")
	dev.new();plot(tBoth2, main="both2")

	return(dj)
	
}

plotJumpVsBranch =	function(
				jj,
				tt	
		)

{

}

getIsler =	function(

		)
{
	common = intersect(as.character(dd$Species),tt)
	which(dd$Species%in%common)
	ddd = dd[which(dd$Species%in%common),]
	rn = rownames(ddd) = sapply(rownames(ddd), as.numeric)
	cn = colnames(ddd) = sapply(colnames(ddd), as.name)
	write.table(ddd,"~/data/expr_phylo/input/isler.table.txt", row.names=rn, col.names=cn)
}

mergeDataSS =	function(	
				dirFp = "",
				outFile = "",
				burnIn = 10 
			)
{
	# read in files from directory into list, indexed by beta values
	fileList = list.files(path=dirFp, pattern="*.p")
	k = length(fileList)
	x = c()
	minLength = Inf 
	
	# correct for unequal # samples between files
	for (i in 1:length(fileList))
	{
		file = read.table(paste(dirFp,fileList[i],sep=""), header=T)
		if (nrow(file) < minLength)
		{
			minLength = nrow(file)
		}
	}

	# save data from burnIn:minLength
	for (i in 1:length(fileList))
	{
		file = read.table(paste(dirFp,fileList[i],sep=""), header=T)
	#	x = c(x, file$lnL[burnIn:minLength])
		x = c(x, file$kb[burnIn:minLength])
	}
	burnLength = minLength - burnIn 
	dim(x) = c(burnLength + 1, k)	
	
	# write to single file
	if (outFile == "")
	{
	#	outFile = strsplit(fileList[1], "\\.")[[1]][1]
		outFile = paste("mcmcmerge", "txt", sep=".")
		outFile = paste(dirFp, outFile, sep="")
	}
	write(t(x), file=outFile, ncolumns=ncol(x), sep="\t")
	return(x)
}

margLikeSS =	function(
					dirFp = "/Users/mlandis/data/expr_phylo/output/",
					mcmcFn = "mcmcmerge.txt",
					betaFn = "beta.txt"
				)
{
	# read in data
	beta = scan(paste(dirFp,betaFn,sep=""))
	x = read.table(paste(dirFp,mcmcFn,sep=""))
	xLen = nrow(x)
	k = length(beta)

	# get max marginal likelihood per beta
	maxLnL = rep(0, k-1)
	for (i in 1:(k-1))
	{
		maxLnL[i] = max(x[,i])
	}
	
	# compute marginal likelihood (non-log scale)
	r = c()	
	for (i in 1:(k-1))
	{
		betaDiff = beta[i+1] - beta[i]
		rK = betaDiff*maxLnL[i] + log(sum(exp(betaDiff*(x[,i] - maxLnL[i])))/xLen)
		r = c(r, rK)
	}

	return(r)
}

mergeMargLikeSS =	function(
				outFp = "/Users/mlandis/data/expr_phylo/output/",
				simName = "",
				mergeFn = "",
				burnIn = 100
			)
{
	dirStr = paste(simName, ".*", sep="")
	dirList = dir(path=outFp, pattern=dirStr)
	r = list() 
	for (i in 1:length(dirList))
	{
		dirFp = paste(outFp,dirList[i],"/",sep="")
		mergeDataSS(dirFp=dirFp,"",burnIn=burnIn)
		fileList = list.files(path=dirFp, pattern="*.p")
		numBeta = length(fileList)
		r[[numBeta]] = margLikeSS(dirFp)
	}
	return(r)
}

sumList = function(vect)
{
	ret = c()
	for (elem in vect) 
	{
		x = sum(elem)
		if (x != 0)
		{
			ret = c(ret,x)
		}
	}
	return(ret)
}

hm = function(L,b=100)
{
	n = length(L)
	M = L[b:n]
	return(length(M) / sum(1/M))
}

