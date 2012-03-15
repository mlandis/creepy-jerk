require(psych)
require(DAAG)

jump_density = function(sigma_bm, sigma_jn, lambda_jn,jumps,tau,tuning_bm=.9995) {
	#WATCH OUT
	return(runif(1))
	#WARNING ABOVE
	w = tau*(1.0-tuning_bm)
	LL = 0
	for (jump in jumps) {
		if (lambda_jn == 0.0) {
			LL = LL + dnorm(jump,0,sigma_jn*sqrt(w),log=T)
			next
		} else {
			for (n in 0:100) {
				LL = LL + dpois(n,lambda_jn*tau)*dnorm(jump,0,sqrt(n*sigma_jn^2+w*sigma_bm^2)) 
			}
		}
	}
	return(log(LL))
}

one_jump = function(sigma_bm,sigma_jn,lambda_jn,jump,tau,tuning_bm=.9995) {
	#WATCH OUT
	return(runif(1))
	#WARNING ABOVE
	w = tau*(1.0-tuning_bm)
	LL = 0
	if (lambda_jn == 0.0) {
		return(dnorm(jump,0,sigma_jn*sqrt(w),log=T))
	} else {
		for (n in 0:100) {
			LL = LL + dpois(n,lambda_jn*tau)*dnorm(jump,0,sqrt(n*sigma_jn^2+w*sigma_bm^2))
		}
	}
	return(log(LL))
}

shape = 5
rate = 1/9
tau = 1.3
#sigma_bm, sigma_jn, lambda_jn, jump
samples = matrix(nrow=50000,ncol=3)
jumps = matrix(nrow=nrow(samples),ncol=10)
#initialize data
samples[1,] = rgamma(3,shape,rate)
jumps[1,] = rnorm(ncol(jumps),0,.5)
#compute initial ratios
k_j = jump_density(samples[1,1],samples[1,2],samples[1,3],jumps[1,],tau)
prior = sum(dgamma(samples[1,1:3],shape,rate,log=T))
for (i in 2:nrow(samples)) {
	if (i %% 500 == 0) {
		print(i)
		print(samples[i-1,])
	}
	#propose sigma_bm
	new_sig = rgamma(1,5,5/samples[i-1,1])
	k_j_new = jump_density(new_sig,samples[i-1,2],samples[i-1,3],jumps[i-1,],tau)
	prior_new = sum(dgamma(c(new_sig,samples[i-1,2],samples[i-1,3]),shape,rate,log=T))
	prop = dgamma(samples[i-1,1],5,5/new_sig,log=T)-dgamma(new_sig,5,5/samples[i-1,1],log=T)
#	ratio = (k_j_new - k_j)+(prior_new-prior)
	ratio = (k_j_new - k_j) + dgamma(new_sig,shape,rate,log=T)-dgamma(samples[i-1,1],shape,rate,log=T) + prop
	u = log(runif(1))
#	print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,1],new_sig))
	if (u < ratio) {
		#accept
#		print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,1],new_sig))
		samples[i,1] = new_sig
		prior = prior_new
		k_j = k_j_new
	} else {
		#reject
		samples[i,1] = samples[i-1,1]
	}
	#propose sigma_jn
	new_sig = rgamma(1,5,5/samples[i-1,2])
	k_j_new = jump_density(samples[i,1],new_sig,samples[i-1,3],jumps[i-1,],tau)
	prior_new = sum(dgamma(c(samples[i,1],new_sig,samples[i-1,3]),shape,rate,log=T))
	prop = dgamma(samples[i-1,2],5,5/new_sig,log=T)-dgamma(new_sig,5,5/samples[i-1,2],log=T)
#	ratio = (k_j_new - k_j)+(prior_new-prior)
	ratio = (k_j_new - k_j) + dgamma(new_sig,shape,rate,log=T)-dgamma(samples[i-1,2],shape,rate,log=T) + prop
	u = log(runif(1))
#	print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,2],new_sig))
	if (u < ratio) {
		#accept
#		print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,2],new_sig))
		samples[i,2] = new_sig
		prior = prior_new
		k_j = k_j_new
	} else {
		#reject
		samples[i,2] = samples[i-1,2]
	}
	#propose lambda_jn
	new_lam = rgamma(1,5,5/samples[i-1,3])
	k_j_new = jump_density(samples[i,1],samples[i,2],new_lam,jumps[i-1,],tau)
	prior_new = sum(dgamma(c(samples[i,1],samples[i,2],new_lam),shape,rate,log=T))
	prop = dgamma(samples[i-1,3],5,5/new_lam,log=T)-dgamma(new_lam,5,5/samples[i-1,3],log=T)
#	ratio = (k_j_new - k_j)+(prior_new-prior)
	ratio = (k_j_new - k_j)+ dgamma(new_lam,shape,rate,log=T)-dgamma(samples[i-1,3],shape,rate,log=T) + prop
	u = log(runif(1))
#	print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,3],new_lam))
	if (u < ratio) {
		#accept
#		print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,3],new_lam))
		samples[i,3] = new_lam
		prior = prior_new
		k_j = k_j_new
	} else {
		#reject
		samples[i,3] = samples[i-1,3]
	}
	#propose jump
	for (j in 1:ncol(jumps)) {
		new_jump = rnorm(1,jumps[i-1,j],.5)
		ratio = one_jump(samples[i,1],samples[i,2],samples[i,3],new_jump,tau)-one_jump(samples[i,1],samples[i,2],samples[i,3],jumps[i-1,j],tau)
		u = log(runif(1))
#		print(c(ratio, u,jumps[i-1,j],new_jump))
		if (u < ratio) {
			#accept
#			print(c(k_j, k_j_new, prior, prior_new, ratio, u,samples[i-1,4],new_jump))
			jumps[i,j] = new_jump
		} else {
			#reject
			jumps[i,j] = jumps[i-1,j]
		}
	}
	#update k_j
	k_j = jump_density(samples[i,1],samples[i,2],samples[i,3],jumps[i,],tau)
} 

