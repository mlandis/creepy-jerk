auto_cor =	function(	data,
				lag
			)
{
	indices = 1:(length(data)-lag)
	return(cor(data[indices],data[indices+lag]))
}

plotMcmc = 	function(	fn,
				b = 100,
				trueparm = c(0.5, 1.2, 0.9),
				plot_auto = F,
				max_lag = 50
			)
{
	x = read.table(fn, header=T)
	xx = x[b:nrow(x),]
	dev.new()
	if (plot_auto)
	{
		par(mfrow=c(4,4),oma=c(0,0,2,0))
	}
	else
	{
		par(mfrow=c(4,3),oma=c(0,0,2,0))
	}
	
	plot(density(xx$Sig.BM,from=0.0))
	abline(v=trueparm[1])
	plot(xx$Sig.BM, xx$lnL)
	plot(xx$Sig.BM,type="l")
	if (plot_auto)
	{
		auto = sapply(1:max_lag,function(i) { auto_cor(xx$Sig.BM,i) })
		plot(auto)	
	}

	plot(density(xx$Lam.JN,from=0.0))
	abline(v=trueparm[2])
	plot(xx$Lam.JN, xx$lnL)
	plot(xx$Lam.JN,type="l")
	if (plot_auto)
	{
		auto = sapply(1:max_lag,function(i) { auto_cor(xx$Lam.JN,i) })
		plot(auto)
	}
	
	plot(density(xx$Sig.JN,from=0.0))
	abline(v=trueparm[3])
	plot(xx$Sig.JN, xx$lnL)
	plot(xx$Sig.JN,type="l")
	if (plot_auto)
	{
		auto = sapply(1:max_lag, function(i) { auto_cor(xx$Sig.JN,i) })
		plot(auto)
	}	

	plot(density(xx$Sig.BM^2 + xx$Lam.JN*xx$Sig.JN^2, from=0.0))
	abline(v=trueparm[1]^2 + trueparm[2]*trueparm[3]^2)
	plot(xx$Sig.BM^2 + xx$Lam.JN*xx$Sig.JN^2, xx$lnL)
	plot(xx$Sig.BM^2 + xx$Lam.JN*xx$Sig.JN^2,type="l")
	if (plot_auto)
	{
		auto = sapply(1:max_lag, function(i) { auto_cor(xx$Sig.BM^2 + xx$Lam.JN*xx$Sig.JN^2,i) })
		plot(auto)
	}
		
	titleName = paste("SigBM=", trueparm[1], ", LamJN=", trueparm[2], ", SigJN=", trueparm[3], sep="")
	title(main=titleName, outer=T)
	
#	dev.new()
#	plot(xx$lnL)
	
	return(xx)
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

