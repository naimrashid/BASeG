ase_joint <- function(Y, Y1, X, X1, Z, output.tag, p.cut, link="log", offset=NULL, 
      adjZ=TRUE, cis.only=TRUE, cis.distance=2e+05, eChr=NULL, ePos=NULL,
      mChr=NULL, mPos=NULL, converge=5e-5, convergeGLM=1e-8, scoreTestP=0.05, 
      trace=1, maxit=1000,
			dChr=NULL, dPos=NULL,cis.distance2=2e+05, nhrule 
= 20
)
{
  
	#here Y is n1, Y1 is nA1, X is n2, X1 is nA2

  if(link=="log"){
    linkR = 2  
  }else if(link=="identity"){
    linkR = 3
    if(adjZ){
      adjZ = FALSE
      warning("adjZ is set to be FALSE for identity link\n")
    }
  }else{
    stop("only 'log' or 'identity' link functions are allowed\n")
  }
  
  ## 
  # check the NAs in X, Y, and Z
  #
  if(any(is.na(Y))){
    stop("NA values in Y\n")
  }
  
  if(any(is.na(X))){
    stop("NA values in X\n")
  }
  
  if(any(is.na(Z))){
    stop("NA values in Z\n")
  }
  
  ## 
  # check the dimensions of X, Y, nA1, and Z
  #
  if(is.vector(Y)){
    nY = 1
    N1 = length(Y)
    if(var(Y) < converge){
      stop("Variance of Y is close to 0\n")
    }
  }else{
    nY = ncol(Y)
    N1 = nrow(Y)
    
    varY = apply(Y, 2, var)
    wVar = which(varY < converge)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }

  if(is.vector(Y1)){
    nY1 = 1
    N12 = length(Y1)
    if(var(Y) < converge){
      stop("Variance of Y is close to 0\n")
    }
  }else{
    nY1 = ncol(Y1)
    N12 = nrow(Y1)
    
    varY1 = apply(Y1, 2, var)
    wVar2 = which(varY1 < converge)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Y1 have tiny variances", length(wVar2))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar2))], collapse=", "))
      
      if(length(wVar2) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }

  
  if(is.vector(X)){
    nX = 1
    N2 = length(X)

    if(var(X) < converge){
      stop("Variance of X is 0\n")
    }
  }else{
    nX = ncol(X)
    N2 = nrow(X)
    
    varX = apply(X, 2, var)
    wVar = which(varX < converge)
    
    #if(length(wVar) > 0){
    if(length(wVar) > 2){
       stpStr = sprintf("%d columns in X have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }

  if(is.vector(X1)){
    nX1 = 1
    N22 = length(X1)

    if(var(X1) < converge){
      stop("Variance of X is 0\n")
    }
  }else{
    nX1 = ncol(X1)
    N22 = nrow(X1)
    
    varX1 = apply(X1, 2, var)
    wVar2 = which(varX1 < converge)
    
    #if(length(wVar) > 0){
    if(length(wVar) > 2){
      stpStr = sprintf("%d columns in X have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }

  
  if(is.vector(Z)){
    nZ = 1
    N3 = length(Z)
    
    if(var(Z) < converge){
      stop("Variance of X is 0\n")
    }
  }else{
    nZ = ncol(Z)
    N3 = nrow(Z)
    
    varZ = apply(Z, 2, var)
    wVar = which(varZ < converge)
    
    if(length(wVar) > 0){
      stpStr = sprintf("%d columns in Z have tiny variances", length(wVar))
      stpStr = sprintf("%s: %s", stpStr, paste(wVar[1:min(10, length(wVar))], collapse=", "))
      
      if(length(wVar) > 10){ stpStr = sprintf("%s ...", stpStr) }
      stop(stpStr, "\n")
    }
  }
    
  if(N1 != N2 || N1 != N3){
    stop("X, Y, and Z must have the same number of rows!")
  }
  
  ## 
  # check the length of offset
  #
  if(!is.null(offset)){
    if(length(offset) != N1){
      stop("length of offset does not match the sample size in Y\n")
    }
    useOffset = 1
  }else{
    useOffset = 0
    offset = rep(0, N1)
  }
  
  ## 
  # check the length of chromosome location information
  #
  if(cis.only){
    if(!is.numeric(eChr) || length(eChr) != nY){
      stop(sprintf("eChr must be a numeric vector of length %d\n", nY))
    }
    if(!is.numeric(ePos) || length(ePos) != nY){
      stop(sprintf("ePos must be a numeric vector of length %d\n", nY))
    }
    if(!is.numeric(dChr) || length(dChr) != nX){
      stop(sprintf("dChr must be a numeric vector of length %d\n", nY1))
    }
    if(!is.numeric(dPos) || length(dPos) != nX){
      stop(sprintf("dPos must be a numeric vector of length %d\n", nY1))
    }
    if(!is.numeric(mChr) || length(mChr) != nZ){
      stop(sprintf("mChr must be a numeric vector of length %d\n", nZ))
    }
    if(!is.numeric(mPos) || length(mPos) != nZ){
      stop(sprintf("mPos must be a numeric vector of length %d\n", nZ))
    }
  }

	library(gaussquad)
	hrule = hermite.h.quadrature.rules(nhrule)[[nhrule]]

  
  dims = numeric(6)
  dims[1] = nY
  dims[2] = nX
  dims[3] = nZ
  dims[4] = N1
  dims[5] = maxit
  dims[6] = useOffset
  dims[7] = nY1 
  dims[8] = nX1 
	dims[9] = nhrule

  succeed = 0
  
  output = character(2)
  output[1] = sprintf("%s_eqtl.txt", output.tag)
  output[2] = sprintf("%s_freq.txt", output.tag)

  # First fit a baseline model using only the confouding covariates 
  # yFailBaselineModel indicate whether the baseline model can be fitted
  
  yFailBaselineModel = numeric(nY)
  yFailBaselineModel2 = numeric(nY1) 

  
  # z is working version of y
  z1 = numeric(N1) 
  z2 = numeric(N2) 

  W = .C("glmEQTL_joint_ase2", as.integer(dims), as.double(Y), as.double(X), as.double(Z), 
         as.double(z1), as.integer(linkR), as.double(offset), 
         as.integer(adjZ), as.character(output), as.double(p.cut), 
         as.integer(cis.only), as.integer(cis.distance), 
         as.integer(eChr), as.integer(ePos), as.integer(mChr), 
         as.integer(mPos), as.double(converge), as.double(convergeGLM), 
         yFailBaselineModel = as.integer(yFailBaselineModel), 
         as.double(scoreTestP), as.integer(trace), 
         succeed=as.integer(succeed), 
				 as.double(Y1), as.double(X1), as.integer(dChr), as.integer(dPos), as.integer(cis.distance2), 
				 as.integer(yFailBaselineModel2), as.double(z2), as.double(hrule$x), as.double(hrule$w) ,
			PACKAGE="BASeG")
    
  list(succeed = W[["succeed"]], yFailBaselineModel=W[["yFailBaselineModel"]])
  
}

