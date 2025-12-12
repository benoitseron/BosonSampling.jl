bosonSampler = function(A, sampleSize, perm = FALSE){
#
# A complex matrix 
#
# out:  samples = list of boson sample mode sequences 
#       perms   = list of associated permanents 
#       pmfs    = list of associated pmfs    
#
# transpose A -- benefit since both R and Armadillo store matrices in column order
  A = t(A)
#	
	m = ncol(A); n = nrow(A)
	modeSeqList = permList = pmfList = c()
	for(k in 1:sampleSize){
#
# permute rows of A, if more than one
#
    if(n>1) A = A[sample(n),]	
#
# start with row 1 and sample the first mode c_1
# i.e. the first mode in an individual sample sequence
#
		firstMode = sample(1:m,size = 1,prob = (Mod(A[1,]))^2)  
		if(n == 1){
			modeSeqList = c(modeSeqList, firstMode)
			if(perm){
				permList = c(permList,A[1,])
			  pmfList = c(pmfList,(Mod(A[1,]))^2)
			}
			next
		}
#
# now sample remaining modes c_2,c_3,...,c_n in an individual sequence
#
		modeSeq = firstMode
		for(modeLimit in 2:n){   
			subPerm = cxPermMinors(matrix(A[1:modeLimit, modeSeq],ncol = modeLimit-1))
#
# permanents of submatrices using Laplace-type expansion
#
			permVector = t(subPerm) %*% A[1:modeLimit,]
#
# next mode in an individual sequence 
#
			nextMode = sample(1:m, size = 1, prob = (Mod(permVector))^2)		
# increment mode sequence
			modeSeq = c(modeSeq, nextMode)		
		}
# build list of sample sequences and associated permanents/pmfs
		modeSeqList = cbind(modeSeqList,modeSeq)
		if(perm){
			finalPerm = permVector[nextMode]
			permList = c(permList,finalPerm)
		  pmfList = c(pmfList,Mod(finalPerm)^2/ prod(factorial(tabulate(modeSeq))))
		}
	}
	list(values = matrix(modeSeqList,nrow = n),perms = permList, pmfs = pmfList) 
}

