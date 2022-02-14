# A few tests suggest that the Eigen implementation makes the code 10-20% faster

setRefClass(
  
  "GEDI",
  
  fields = list(
    params="list", # This is where the main parameters of the model are stored
    hyperparams="list", # This is where the hyper-parameters of the model are stored
    aux="list", # This is where the auxiliary variables that are useful for model fitting are stored
    target="list", # This is where the target matrices (Y and/or M) are stored
    tracking="list" # This is where the tracking information for model fitting progress is stored
  ),
  
  methods = list(
    
    #######################
    # Function to create the matrices and prepare other variables for a GEDI model
    #######################
    
    setup = function(
      Samples, # A factor variable indicating the sample of origin of each cell
      Y = NULL, # The log-transformed (and possibly normalized) gene expression matrix
      M = NULL, # The raw read count matrix; overrides Y if provided. It can also be a list of two matrices, in which case they are considered as paired observations whose log-ratio must be modelled
      C = NULL, # The gene-level biological prior. If NULL, it means that there is no prior for Z
      H = NULL, # The gene-level prior for unwanted sources of variation. If NULL, there will be no prior for Qi
      K = 10, # The number of latent variables
      mode = "Bl2", # Three values are allowed: "BL2" [L2 norm of the entire B matrix is fixed], or "Bsphere" [L2 norms of B columns are fixed]
      adjustD = T, # Whether D should be adjusted based on the L2 norm of B rows (TRUE), or just use the default (FALSE)
      orthoZ = T, # Whether the columns of Z should be orthogonal
      Z_shrinkage = 1, # The shrinkage multiplier for Z
      A_shrinkage = 1, # The shrinkage multiplier for A
      Qi_shrinkage = 1, # The shrinkage multiplier for Qi
      Ri_shrinkage = 1, # The shrinkage multiplier for Ri
      oi_shrinkage = 1, # The shrinkage multiplier for oi
      o_shrinkage = 1, # The shrinkage multiplier for o
      si_shrinkage = 1 # The shrinkage multiplier for si
    ) {
      
      message("Setting up the GEDI model...")
      
      ##### initialize the observation type

      # Note: M will be stored in the same format as the one provided (e.g. dgCMatrix)
      #  However, Y is always stored as `matrix`
      if( is.null(M) & is.null(Y) ) { stop("Error: Either Y or M should be provided") }
      if( !is.null(M) ) {
        if( is.list(M) ) { # M is a list
          if( length(M)!=2 ) {
            stop("If M is a list, it must be a list of exactly two sparse matrices")
          } else {
            if( !identical(rownames(M[[1]]),rownames(M[[2]])) |
                !identical(colnames(M[[1]]),colnames(M[[2]])) ) {
              stop("If M is a list, it must be a list of exactly two sparse matrices with matching rows and columns")
            }
            Y <- as.matrix( log( (1+M[[1]])/(1+M[[2]]) ) ) # if M is a pair of matrices, Y will be initialized (or over-written) as the log ratio of the two matrices
            aux$obs.type <<- "M_paired"
            aux$geneIDs <<- rownames(M[[1]])
            aux$cellIDs <<- colnames(M[[1]])
          }
        } else { # M is provided but is not a list
          Y <- as.matrix( log(1+M) ) # if M is provided as a single matrix, Y will be initialized (or over-written) as log of M
          aux$obs.type <<- "M"
          aux$geneIDs <<- rownames(M)
          aux$cellIDs <<- colnames(M)
        }
      } else { # M is not provided, so use Y
        Y <- as.matrix(Y) # To ensure `matrix` type
        aux$obs.type <<- "Y"
        aux$geneIDs <<- rownames(Y)
        aux$cellIDs <<- colnames(Y)
      }
      
      # Also store the mode that will be used for dimensionality reduction
      aux$mode <<- mode
      aux$adjustD <<- adjustD
      aux$orthoZ <<- orthoZ
      
      ##### Store the model dimensions and check the compatibility of different matrices
      
      aux$J <<- nrow(Y)
      aux$J_vec <<- rep(1,aux$J)
      aux$N <<- ncol(Y)
      aux$N_vec <<- rep(1,aux$N)
      if( aux$N != length(Samples) ) { stop("Incompatible matrices.") }
      
      if( !is.null(C) ) {
        aux$inputC <<- as.matrix(C) # To ensure `matrix` type
        minCL2 <- min(colSums(aux$inputC^2)) # calculate the L2 of the smallest column of C
        svdC <- svd(aux$inputC) # Perform SVD on C to obtain uncorrelated components
        aux$P <<- sum( cumsum(rev(svdC$d^2)) > minCL2*0.1 ) # only consider the SVD components that have at least 10% of the variance of the smallest column of C
        aux$C <<- svdC$u[,1:aux$P] # The C matrix used by GEDI is in fact the U matrix of SVD(C)
        aux$C.rotation <<- svdC$v[,1:aux$P] %*% diag(1/svdC$d[1:aux$P],nrow=aux$P) # Rotation matrix, so that aux$inputC %*% aux$C.rotation == aux$C
        if( aux$J != nrow(aux$C) ) { stop("Incompatible matrices.") } 
      } else {
        aux$C <<- matrix(nrow=0,ncol=0)
        aux$P <<- 0
      }
      
      if( !is.null(H) ) {
        aux$inputH <<- as.matrix(H) # To ensure `matrix` type
        minHL2 <- min(colSums(aux$inputH^2)) # calculate the L2 of the smallest column of H
        svdH <- svd(aux$inputH) # Perform SVD on H to obtain uncorrelated components
        aux$L <<- sum( cumsum(rev(svdH$d^2)) > minHL2*0.1 ) # only consider the SVD components that have at least 10% of the variance of the smallest column of H
        aux$H <<- svdH$u[,1:aux$L] # The H matrix used by GEDI is in fact the U matrix of SVD(H)
        aux$H.rotation <<- svdH$v[,1:aux$L] %*% diag(1/svdH$d[1:aux$L],nrow=aux$L) # Rotation matrix, so that aux$inputH %*% aux$H.rotation == aux$H
        if( aux$J != nrow(aux$H) ) { stop("Incompatible matrices.") } 
      } else {
        aux$H <<- matrix(nrow=0,ncol=0)
        aux$L <<- 0
      }
      
      aux$K <<- K
      aux$diag_K <<- diag(K)

      # store additional auxiliary information
      
      aux$BcolL2 <<- 1 # The expected L2 norm of each column of B in the Bsphere mode.
      aux$BrowL2 <<- 1 # The expected L2 norm of each row of B in the Bnorm mode

      aux$Samples <<- unique( Samples ) # The list of unique samples
      aux$numSamples <<- length(aux$Samples)
      aux$cells <<- vector("list",aux$numSamples) # The indexes of cells originating from each sample
      aux$Ni <<- rep(NA,aux$numSamples) # The number of cells originating from each sample
      aux$Ni_vec <<- vector("list",aux$numSamples) # A list containing a vector of 1's for each sample
      for( i in 1:aux$numSamples ) {
        aux$cells[[i]] <<- which( Samples==aux$Samples[i] )
        aux$Ni[i] <<- length(aux$cells[[i]])
        aux$Ni_vec[[i]] <<- rep(1,aux$Ni[i])
      }
      
      aux$ite <<- 0 # the number of optimization iterations completed so far
      
      ##### initialize hyperparameters
      
      # initialize o_0 and s_0, the priors for gene-level and cell-level scaling factors
      if( is.null(M) ) # M is not provided, so s and o should be calculated directly from Y
      {
        s_0 <- eigenVecMatProduct( aux$J_vec, Y ) / aux$J # s is initialized as the column-mean of Y
        Yp <- Y - eigenVecVecProduct( aux$J_vec, s_0 ) # Yp is the library size-normalized version of Y
        hyperparams$o_0 <<- eigenMatVecProduct(Yp,aux$N_vec) / aux$N # o is initialized as the row-mean of Y, after normalizing for s
        
      } else if( is.list(M) ) { # M is a pair of matrices
        s1_0 <- eigenVecMatProduct( aux$J_vec, as.matrix(M[[1]]) ) / aux$J # s1 is the column-mean of M1
        M1p <- as.matrix(M[[1]]) / eigenVecVecProduct( aux$J_vec, s1_0 ) # M1p is the library size-normalized version of M1
        o1_0 <- eigenMatVecProduct(M1p,aux$N_vec) / aux$N # o1 is the row-mean of M1p

        s2_0 <- eigenVecMatProduct( aux$J_vec, as.matrix(M[[2]]) ) / aux$J # s2 is the column-mean of M2
        M2p <- as.matrix(M[[2]]) / eigenVecVecProduct( aux$J_vec, s2_0 ) # M2p is the library size-normalized version of M2
        o2_0 <- eigenMatVecProduct(M2p,aux$N_vec) / aux$N # o2 is the row-mean of M2p
        
        # Now calculate s_0 and o_0. Note that s_0 is initially saved as a local object, but later it is split by sample and stored in `params`
        s_0 <- log(s1_0/s2_0) # s is initialized as the log ratio of the column means of the paired measurements
        hyperparams$o_0 <<- log(o1_0/o2_0) # o is initialized as the log ratio of the row means of the paired measurements, after normalizing for library size
        
      } else { # M is provided as a single matrix
        s_0 <- eigenVecMatProduct( aux$J_vec, as.matrix(M) ) / aux$J # s is initialized as the column-mean of M
        Mp <- as.matrix(M) / eigenVecVecProduct( aux$J_vec, s_0 ) # Mp is the library size-normalized version of M
        hyperparams$o_0 <<- eigenMatVecProduct(Mp,aux$N_vec) / aux$N # o is initialized as the row-mean of M, after normalizing for s
        
        # convert s_0 and o_0 to log-scale, which is compatible with Y
        # Note that s_0 is initially saved as a local object, but later it is split by sample and stored in `params`
        s_0 <- log(s_0)
        hyperparams$o_0 <<- log(hyperparams$o_0)
      }
      # Here is where s is split by sample and stored
      hyperparams$si_0 <<- vector("list",aux$numSamples)
      for( i in 1:aux$numSamples ) {
        hyperparams$si_0[[i]] <<- s_0[aux$cells[[i]]]
      }
      
      # Initialize sigma2_0
      Yp <- Y - eigenVecVecProduct( aux$J_vec, s_0 ) - hyperparams$o_0
      hyperparams$sigma2_0 <<- sum(c(Yp*Yp)) / (aux$N*aux$J+1)
      ### This is just a test to see what happens
      #hyperparams$sigma2_0 <<- 0.01

      # Initialize other hyperparameters
      hyperparams$S_o <<- 1/aux$N/o_shrinkage # The prior relative variance for o
      hyperparams$S_si <<- 1/aux$J/si_shrinkage # The prior relative variance for si
      hyperparams$S_Z <<- 1/Z_shrinkage # The prior relative variance for Z
      hyperparams$S_A <<- 1/A_shrinkage # The prior relative variance for A
      hyperparams$S_Ri <<- 1/Ri_shrinkage # The prior relative variance for Ri
      hyperparams$S_oi <<- 1/aux$Ni/oi_shrinkage # The prior relative variance for oi
      hyperparams$S_Qi <<- aux$N/aux$Ni/Qi_shrinkage # The prior relative variance for Qi
      
      
      ##### Initialize model parameters and matrices
      
      params$sigma2 <<- hyperparams$sigma2_0 # Initialize sigma2 as its prior mean
      params$o <<- hyperparams$o_0  # Initialize o as its prior mean
      params$si <<- hyperparams$si_0 # Initialize si as its prior mean 
      
      # Initialize Z and D
      params$Z <<- matrix(0,nrow=aux$J,ncol=K)
      if( mode=="Bl2" ){
        # since the L2 norm of each row is fixed at 1, D becomes equivalent to the expected row-wise L2 norm
        params$D <<- rep(aux$BrowL2,K)
      } else if( mode=="Bsphere" ) {
        # since L2 norm of each column is fixed at 1, on average each entry is 1/sqrt(K), and the expected row-wise L2  norm without adjustment is N/K --> each entry must be multiplied by sqrt(K/N)
        params$D <<- sqrt(aux$BrowL2) * rep(sqrt(K/aux$N),K)
      } else {
        stop("Unrecognized mode")
      }
      # Intialize U and S. U will be an orthonormal matrix, and S a diagonal matrix, so that Z=US
      params$U <<- matrix(0,nrow=aux$J,ncol=K)
      params$S <<- rep(1,K)
      
      # Initialize A
      if( is.null(C) ) {
        params$A <<- matrix(nrow=0,ncol=0)
      } else {
        params$A <<- matrix(0,nrow=aux$P,ncol=K)
      }
      
      # Initialize sample-specific matrices
      target$Mi <<- vector("list",aux$numSamples)
      target$M1i <<- vector("list",aux$numSamples)
      target$M2i <<- vector("list",aux$numSamples)
      target$Yi <<- vector("list",aux$numSamples)
      params$oi <<- vector("list",aux$numSamples)
      params$Bi <<- vector("list",aux$numSamples)
      params$Qi <<- vector("list",aux$numSamples)
      params$Ri <<- vector("list",aux$numSamples)
      
      for( i in 1:aux$numSamples ) {
        # Copy Mi
        if( is.null(M) ) { # No M is provided
          target$Mi[[i]] <<- matrix(nrow=0,ncol=0)
          target$M1i[[i]] <<- matrix(nrow=0,ncol=0)
          target$M2i[[i]] <<- matrix(nrow=0,ncol=0)
        } else if( is.list(M) ) { # M is a pair of matrices
          target$Mi[[i]] <<- matrix(nrow=0,ncol=0)
          target$M1i[[i]] <<-  M[[1]][,aux$cells[[i]],drop=F]
          target$M2i[[i]] <<-  M[[2]][,aux$cells[[i]],drop=F]
        } else { # M is a single matrix
          target$Mi[[i]] <<-  M[,aux$cells[[i]],drop=F]
          target$M1i[[i]] <<- matrix(nrow=0,ncol=0)
          target$M2i[[i]] <<- matrix(nrow=0,ncol=0)
        }
        target$Yi[[i]] <<-  Y[,aux$cells[[i]],drop=F] # Copy Yi
        params$oi[[i]] <<- rep(0,aux$J) # Initialize oi with zero
        params$Bi[[i]] <<- matrix(0,nrow=K,ncol=aux$Ni[i]) # Initialize Bi with zero
        params$Qi[[i]] <<- matrix(0,nrow=aux$J,ncol=K) # Initialize Qi with zero
        aux$ZDBi[[i]] <<- matrix(0,nrow=aux$J,ncol=aux$Ni[i]) # ZDBi stores Z x Bi product
        aux$QiDBi[[i]] <<- matrix(0,nrow=aux$J,ncol=aux$Ni[i]) # QiDBi stores Qi x Bi product
        if( is.null(H) ) {
          params$Ri[[i]] <<- matrix(nrow=0,ncol=0)
        } else {
          params$Ri[[i]] <<- matrix(0,nrow=aux$L,ncol=K)
        }
      }
      
      # Initialize tracking vectors
      tracking$dZ <<- NA # This vector tracks the convergence of Z
      tracking$dA <<- NA # This vector tracks the convergence of A
      tracking$do <<- NA # This vector tracks the convergence of o
      
      tracking$sigma2 <<- NA # This vector tracks the value of sigma2 across iterations
      
      tracking$dsi <<- vector("list",aux$numSamples) # This list tracks the convergence of si for each sample
      tracking$doi <<- vector("list",aux$numSamples) # This list tracks the convergence of oi for each sample
      tracking$dBi <<- vector("list",aux$numSamples) # This list tracks the convergence of Bi for each sample
      tracking$dQi <<- vector("list",aux$numSamples) # This list tracks the convergence of Qi for each sample
      tracking$dRi <<- vector("list",aux$numSamples) # This list tracks the convergence of Ri for each sample
      for( i in 1:aux$numSamples ) {
        tracking$dsi[[i]] <<- NA
        tracking$doi[[i]] <<- NA
        tracking$dBi[[i]] <<- NA
        tracking$dQi[[i]] <<- NA
        tracking$dRi[[i]] <<- NA
      }
    },

    
    
    #######################
    # Functions for fitting individual components of the model (by block coordinate optimization)
    #######################
    
    
    ######
    # Solves for Bi, the low-dimensional projection of cells from sample i: Yi ~ o + oi + (Z+Qi)xBi
    #   This function is used in the Bl2 and Bsphere modes, which do not restrict the sign of B
    # `i`: the index of the sample for which Bi is being optimized
    solve.Bi = function( i ) {
      # Call the Rcpp function that performs the actual computation
      params$Bi[[i]] <<- solveBi(
        target$Yi[[i]],
        params$D,
        params$Z, params$Qi[[i]],
        aux$diag_K,
        params$si[[i]], params$o, params$oi[[i]] )
      
      #return(1)
    },
    
    
    ######
    # Normalized Bi for all samples, depending on the mode of the GEDI model
    normalize.B = function() {
      if( aux$mode == "Bsphere" ) {
        # Re-normarlize each Bi to get the desired L2 norm for each column (can be done separately for each Bi)
        for( i in 1:aux$numSamples ) {
          BcolL2 <- eigenColL2( params$Bi[[i]] ) # calculate the L2 norm of each column
          scaling_factors <- sqrt(aux$BcolL2/(BcolL2+1e-10)) # calculate the scaling factor for achieving the expected L2 norm
          params$Bi[[i]] <<- eigenRowMult( params$Bi[[i]], scaling_factors ) # re-normalize each column
        }
      }

      # Now, calculate the current L2 norm of each row
      BrowL2 <- rep(0,aux$K)
      for( i in 1:aux$numSamples ) {
        BrowL2 <-
          BrowL2 + eigenRowL2( params$Bi[[i]] )
      }
      # calculate the scaling factor that would bring the L2 norms to 1
      scaling_factors <- sqrt( aux$BrowL2/BrowL2 )
      
      if( aux$mode == "Bl2" ) { # In the Bl2 mode, the L2 norm of each row of B is fixed.
        # Therefore, re-normarlize each Bi to get the desired L2 norm for each row of the entire B
        for( i in 1:aux$numSamples ){
          params$Bi[[i]] <<- eigenMatProduct(
            diag(scaling_factors,nrow=aux$K ), params$Bi[[i]] ) # re-normalize each row
        }
      } else if(aux$adjustD) { # In other modes, the D diagonal matrix ensures that D%*%B has row-wise the desired L2 norm
        params$D <<- scaling_factors
      }
        
      #return(1)
    },
    
    
    ######
    # Update ZDBi and QiDBi for all samples
    update.ZDBi = function() {
      for( i in 1:aux$numSamples ) {
        aux$ZDBi[[i]] <<- eigenMatProduct(
          params$Z%*%diag(params$D,nrow=aux$K),
          params$Bi[[i]] )
      }
      #return(1)
    },
    update.QiDBi = function() {
      for( i in 1:aux$numSamples ) {
        aux$QiDBi[[i]] <<- eigenMatProduct(
          params$Qi[[i]]%*%diag(params$D,nrow=aux$K),
          params$Bi[[i]] )
      }
      #return(1)
    },
    
    ######
    # Solves for Z, the shared metagene matrix: Yi ~ o + oi + (Z+Qi) x Bi
    # This function does not enforce orthogonality of Z
    # Note that this function requires QiDBi to be up-to-date
    solve.Z = function() {
      # For each sample, calculate the residual of Y after considering everything
      #  other than ZDBi; then concatenate the residuals; also concatenate Bi's to get B
      Y_res <- NULL
      B <- NULL
      for( i in 1:aux$numSamples )
      {
        Y_res <- cbind( Y_res, Yi_resZ(
          target$Yi[[i]],
          aux$QiDBi[[i]],
          params$si[[i]], params$o, params$oi[[i]] ) )
        B <- cbind( B, params$Bi[[i]] )
      }
      
      # Now solve Z by regression on the residual
      lambda <- 1/hyperparams$S_Z
      if( aux$P==0 ) # no C matrix is provided
      {
        params$Z <<- solveZ_noC(
          Y_res, params$D, B, aux$diag_K, lambda )
      } else { # matrix C is provided
        params$Z <<- solveZ_wC(
          Y_res, params$D, B, aux$C, params$A, aux$diag_K, lambda )
      }
      #return(1)
    },
    
    ######
    # Solves for Z, the shared metagene matrix: Yi ~ o + oi + (Z+Qi) x Bi
    # This function enforces orthogonality of Z, by finding an orthonormal U and diagonal S so that Z=US
    # Note that this function requires QiDBi to be up-to-date
    solve.Z.orthogonal = function() {
      # For each sample, calculate the residual of Y after considering everything
      #  other than ZDBi; then concatenate the residuals; also concatenate Bi's to get B
      Y_res <- NULL
      B <- NULL
      for( i in 1:aux$numSamples )
      {
        Y_res <- cbind( Y_res, Yi_resZ(
          target$Yi[[i]],
          aux$QiDBi[[i]],
          params$si[[i]], params$o, params$oi[[i]] ) )
        B <- cbind( B, params$Bi[[i]] )
      }
      B <- eigenMatProduct(diag(params$D,nrow=aux$K),B)
      
      lambda <- 1/hyperparams$S_Z

      # If C is provided, use C%*%A as prior
      if( aux$P > 0 ) {
        Y_res <- cbind( Y_res, sqrt(lambda)*eigenMatProduct(aux$C,params$A) )
        B <- cbind( B, sqrt(lambda)*aux$diag_K )
      }
      
      # Now solve U by SVD
      params$U <<- solveU(Y_res,params$S,B)
      
      # Solve S, one element at a time
      YBtU <- eigenMatTcrossprod(Y_res,B) * params$U
      params$S <<- eigenVecMatProduct(aux$J_vec,YBtU) / (1+lambda)
      
      # Update Z based on the new U and S
      params$Z <<- eigenMatProduct( params$U, diag(params$S,nrow=aux$K) )
      #return(1)
    },
    
    
    ######
    # Solves for A, the matrix that connects the prior, C, to the shared metagenes Z: Z ~ C x A
    solve.A = function() {
      lambda <- 1/hyperparams$S_A
      params$A <<- eigenMatProduct(
        solve( eigenMatCrossprod(aux$C,aux$C) + lambda*diag(aux$P) ), # The code can be optimized by pre-calculating this part
        eigenMatCrossprod(aux$C,params$Z) )
      #return(1)
    },
    
    
    ######
    # Solves for Qi, the sample-specific matrix of metagenes for sample i: Yi ~ o + oi + (Z+Qi)xBi
    # `i`: the index of the sample for which Qi is being updated
    # Note that this function requires ZDBi to be up-to-date
    solve.Qi = function( i ) {
      lambda <- 1/hyperparams$S_Qi[i]
      if( aux$L==0 ) { # no matrix H is provided
        params$Qi[[i]] <<- solveQi_noH(
          target$Yi[[i]],
          aux$ZDBi[[i]],
          params$si[[i]], params$o, params$oi[[i]],
          params$D, params$Bi[[i]],
          aux$diag_K, lambda )
      } else { # matrix H is provided
        params$Qi[[i]] <<- solveQi_wH(
          target$Yi[[i]],
          aux$ZDBi[[i]],
          params$si[[i]], params$o, params$oi[[i]],
          params$D, params$Bi[[i]],
          aux$H, params$Ri[[i]],
          aux$diag_K, lambda )
      }
      #return(1)
    },
    
    
    ######
    # Solves for Ri, the matrix that connects the prior H to the sample-specific metagenes Qi: Qi ~ H x Ri
    # `i`: the index of the sample for which Ri is being updated
    solve.Ri = function( i ) {
      lambda <- 1/hyperparams$S_Ri
      params$Ri[[i]] <<- eigenMatProduct(
        solve( eigenMatCrossprod(aux$H,aux$H) + lambda*diag(aux$L) ), # The code can be optimized by pre-calculating this part
        eigenMatCrossprod(aux$H,params$Qi[[i]]) )
      #return(1)
    },
    
    
    ######
    # Solve for oi, the sample-specific offset: Yi ~ o + oi + (Z+Qi) x Bi
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZDBi and QiDBi to be up-to-date
    solve.oi = function( i )
    {
      lambda <- 1/hyperparams$S_oi[i]
      params$oi[[i]] <<- solveOi(
        target$Yi[[i]],
        aux$ZDBi[[i]], aux$QiDBi[[i]],
        aux$Ni_vec[[i]],
        params$si[[i]], params$o,
        aux$Ni[[i]], lambda )
      #return(1)
    },
    
    
    ######
    # Solves for si, the cell-specific library size
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZDBi and QiDBi to be up-to-date
    solve.si = function( i )
    {
      lambda <- 1/hyperparams$S_si
      params$si[[i]] <<- solveSi(
        target$Yi[[i]],
        aux$ZDBi[[i]], aux$QiDBi[[i]],
        aux$J_vec,
        params$o, params$oi[[i]],
        hyperparams$si_0[[i]],
        aux$J, lambda )
      #return(1)
    },
    
    
    ######
    # Solves for o, the shared global offset: Yi ~ o + oi + (Z+Qi) x Bi
    # Note that this function requires ZDBi and QiDBi to be up-to-date for all samples
    solve.o = function()
    {
      # For each sample, calculate the residual of Yi after considering everything
      #  other than o, calculate row-wise sum of the residuals, aggregate across all samples,
      #  and then use to calculate o (while considering the prior)
      params$o <<- rep(0,aux$J)
      for( i in 1:aux$numSamples )
      {
        # Calculate residual
        params$o <<- params$o + Yi_resO_rowSum(
          target$Yi[[i]],
          aux$ZDBi[[i]], aux$QiDBi[[i]],
          aux$Ni_vec[[i]],
          params$si[[i]], params$oi[[i]] )
      }
      # Now solve o
      lambda <- 1/hyperparams$S_o
      params$o <<- ( params$o + hyperparams$o_0*lambda ) / ( aux$N + lambda )
      #return(1)
    },
    
    
    ######
    # Solve for sigma2, the mean squared error of the model, for the case where M is observed
    # Note that this function requires ZDBi and QiDBi to be up-to-date for all samples
    solve.sigma2 = function()
    {
      # The denominator of sigma2
      N <-
        aux$J*(aux$N+aux$numSamples) + # Y and oi
        aux$N + # si
        aux$J + #o
        aux$numSamples*aux$K*(aux$J+aux$L) + # Qi and Ri
        aux$K*(aux$J+aux$P) # Z and A
      # For the numerator, first calculate the weighted sum of all coefficients and residuals
      S <- vecL2_wPrior(params$o,hyperparams$o_0) / hyperparams$S_o
      if( aux$P==0 ) # matrix C is not provided
      {
        S <- S + matL2_noPrior(params$Z) / hyperparams$S_Z
      } else { # matrix C is provided
        S <- S +
          matL2_wPrior(params$Z,aux$C,params$A) / hyperparams$S_Z +
          matL2_noPrior(params$A) / hyperparams$S_Z / hyperparams$S_A
      }
      for( i in 1:aux$numSamples ) # Add the sample-specific coefficients
      {
        S <- S + vecL2_noPrior(params$oi[[i]]) / hyperparams$S_oi[i]
        if( aux$L==0 ) # matrix H is not provided
        {
          S <- S + matL2_noPrior(params$Qi[[i]]) / hyperparams$S_Qi[i]
        } else { # matrix H is provided
          S <- S +
            matL2_wPrior(params$Qi[[i]],aux$H,params$Ri[[i]]) / hyperparams$S_Qi[i] +
            matL2_noPrior(params$Ri[[i]]) / hyperparams$S_Qi[[i]] / hyperparams$S_Ri
        }
        S <- S + vecL2_wPrior(params$si[[i]],hyperparams$si_0[[i]]) / hyperparams$S_si
        # Calculate the residual of Yi, and add it to S
        S <- S + Yi_SSE(
          target$Yi[[i]],
          aux$ZDBi[[i]], aux$QiDBi[[i]],
          params$si[[i]], params$o, params$oi[[i]],
          params$sigma2 )
      }
      params$sigma2 <<- S/N
      tracking$sigma2[aux$ite] <<- params$sigma2
      cat(params$sigma2,"\n")
    },

    
    ######
    # Solve for sigma2, the mean squared error of the model, for the case where M1 and M2 are observed
    # Note that this function requires ZDBi and QiDBi to be up-to-date for all samples
    solve.sigma2.paired = function()
    {
      # The denominator of sigma2
      N <-
        aux$J*(aux$N+aux$numSamples) + # Y and oi
        aux$N + # si
        aux$J + #o
        aux$numSamples*aux$K*(aux$J+aux$L) + # Qi and Ri
        aux$K*(aux$J+aux$P) # Z and A
      # For the numerator, first calculate the weighted sum of all coefficients and residuals
      S <- vecL2_wPrior(params$o,hyperparams$o_0) / hyperparams$S_o
      if( aux$P==0 ) # matrix C is not provided
      {
        S <- S + matL2_noPrior(params$Z) / hyperparams$S_Z
      } else { # matrix C is provided
        S <- S +
          matL2_wPrior(params$Z,aux$C,params$A) / hyperparams$S_Z +
          matL2_noPrior(params$A) / hyperparams$S_Z / hyperparams$S_A
      }
      for( i in 1:aux$numSamples ) # Add the sample-specific coefficients
      {
        S <- S + vecL2_noPrior(params$oi[[i]]) / hyperparams$S_oi[i]
        if( aux$L==0 ) # matrix H is not provided
        {
          S <- S + matL2_noPrior(params$Qi[[i]]) / hyperparams$S_Qi[i]
        } else { # matrix H is provided
          S <- S +
            matL2_wPrior(params$Qi[[i]],aux$H,params$Ri[[i]]) / hyperparams$S_Qi[i] +
            matL2_noPrior(params$Ri[[i]]) / hyperparams$S_Qi[[i]] / hyperparams$S_Ri
        }
        S <- S + vecL2_wPrior(params$si[[i]],hyperparams$si_0[[i]]) / hyperparams$S_si
        # Calculate the residual of Yi, and add it to S
        S <- S + Yi_SSE_paired(
          target$Yi[[i]],
          target$M1i[[i]], target$M2i[[i]],
          aux$ZDBi[[i]], aux$QiDBi[[i]],
          params$si[[i]], params$o, params$oi[[i]],
          params$sigma2 )
      }
      params$sigma2 <<- S/N
      tracking$sigma2[aux$ite] <<- params$sigma2
      cat(params$sigma2,"\n")
    },

    
    ######
    # Solve Yi using raw counts Mi
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZDBi and QiDBi to be up-to-date
    solve.Yi = function( i )
    {
      # Calculate the estimate of Yi based on model parameters
      target$Yi[[i]] <<- solveYi(
        target$Mi[[i]],
        target$Yi[[i]],
        aux$ZDBi[[i]], aux$QiDBi[[i]],
        params$si[[i]], params$o, params$oi[[i]],
        params$sigma2 )
      
      #return(1)
    },
    

    ######
    # Solve Yi using a pair of raw count matrices M1i and M2i
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZDBi and QiDBi to be up-to-date
    solve.Yi.paired = function( i )
    {
      # Calculate the estimate of Yi based on model parameters
      target$Yi[[i]] <<- solveYi_paired(
        target$M1i[[i]], target$M2i[[i]],
        target$Yi[[i]],
        aux$ZDBi[[i]], aux$QiDBi[[i]],
        params$si[[i]], params$o, params$oi[[i]],
        params$sigma2 )
      
      #return(1)
    },
    
    
    #######################
    # Function to initialize the Z, B, oi, and Qi matrices
    #######################
    
    # `randomSeed`: if a random seed number is provided, it will be used in
    #   set.seed before rsvd or Kmeans
    initialize.LVs = function( randomSeed=NA ) {
      
      message("Initializing LVs...")
      
      # Initialize oi and calculate the residual after taking into
      #   account o, oi, and si
      message("  Initializing oi...")
      Yp <- NULL
      for( i in 1:aux$numSamples ) {
        solve.oi(i)
        Yp <- cbind( Yp,
          target$Yi[[i]] -
            eigenVecVecProduct(aux$J_vec,params$si[[i]]) -
            params$o - params$oi[[i]] )
      }
      
      # Now use the residual to initialize Z, Qi, and Bi
      message("  Performing initial decompsition of Y...")
      
      if( aux$mode=="Bl2" | aux$mode=="Bsphere" ) {
        # Use SVD decomposition to obtain an initial estimate of Z
        if( !is.na(randomSeed) ) { set.seed(randomSeed) } # To ensure reproducible results
        svdY <- rsvd( Yp, k = aux$K, nu = aux$K, nv = aux$K )
        params$Z <<- eigenMatProduct( svdY$u, diag(svdY$d[1:aux$K],nrow=aux$K) ) / 2 # Z is divided by 2, and the other half is later copied to each Qi
        params$U <<- svdY$u
        params$S <<- svdY$d[1:aux$K]/2
        # Solve Bi given the initial Z (and Qi) estimates
        for( i in 1:aux$numSamples ) {
          params$Qi[[i]] <<- params$Z
          solve.Bi(i)
        }
      } else {
        stop("Unsupported mode.")
      }
      # Normalize B, and update the ZDBi and QiDBi matrices
      normalize.B()
      # In Bl2 mode, normalization is not expected to have an effect
      # Also in the Bsphere mode, since B had L2 norm of 1 in the beginning,
      # normalization would only bring it back to the original scale.
      
      update.ZDBi()
      update.QiDBi()
    },

    
    #######################
    # Function for tracking changes in the model parameters
    #######################
    
    # `ref`: the reference relative to which the parameter changes are calculated
    track_changes = function( ref ) {
      tracking$do[aux$ite] <<- vecRMSD(params$o,ref$o)
      tracking$dZ[aux$ite] <<- matRMSD(params$Z,ref$Z)
      if(aux$P>0) {
        tracking$dA[aux$ite] <<- matRMSD(params$A,ref$A)
      }
      for( i in 1:aux$numSamples ) {
        tracking$dsi[[i]][aux$ite] <<- vecRMSD(params$si[[i]],ref$si[[i]])
        tracking$doi[[i]][aux$ite] <<- vecRMSD(params$oi[[i]],ref$oi[[i]])
        tracking$dBi[[i]][aux$ite] <<- matRMSD(params$Bi[[i]],ref$Bi[[i]])
        tracking$dQi[[i]][aux$ite] <<- matRMSD(params$Qi[[i]],ref$Qi[[i]])
        if(aux$L>0) {
          tracking$dRi[[i]][aux$ite] <<- matRMSD(params$Ri[[i]],ref$Ri[[i]])
        }
      }
    },
    
    
    #######################
    # Function for iterative optimization of the model
    #######################
    
    # `iterations`: the number of iterations for this round of optimization
    # `track_interval`: the interval for tracking statistics to be calculated and stored
    optimize = function( iterations=50, track_internval=5 ) {
      
      message("Performing block coordinate descent optimization...")
      
      # iteratively optimize the model through block coordinate descent
      for( ite in 1:iterations ) {
        aux$ite <<- aux$ite+1
        message("  Iteration ",ite,"/",iterations," (total ",aux$ite,")...")
        # if tracking info needs to be updated, store the current model parameters
        if( aux$ite %% track_internval == 1 ) {
          prev_params <- params
        }
        # Optimize Bi for all samples
        if( aux$mode=="Bl2" | aux$mode=="Bsphere" ) {
          for( i in 1:aux$numSamples ) { solve.Bi(i) }
        } else {
          stop("Unsupported mode.")
        }
        normalize.B() # normalize B after solving all Bi matrices
        # Solve Z. Note that in order to solve Z, QiDBi needs to be updated
        update.QiDBi()
        if( aux$orthoZ ) { # Z columns should be orthogonal
          solve.Z.orthogonal()
        } else {
          solve.Z()
        }
        # Solve Qi. Note that ZDBi needs to be updated first.
        update.ZDBi()
        for( i in 1:aux$numSamples ) { solve.Qi(i) }
        # Since Qi has changed, QiDBi needs to be updated again for downstream steps
        update.QiDBi()
        # Solve oi for all samples
        for( i in 1:aux$numSamples ) { solve.oi(i) }
        # Solve si for all samples
        for( i in 1:aux$numSamples ) { solve.si(i) }
        # Solve o
        solve.o()
        # If necessary, solve Yi for all samples, and then update sigma2
        if( aux$obs.type=="M" ) {
          for( i in 1:aux$numSamples ) { solve.Yi(i) }
          solve.sigma2() # since sigma2 is only needed for solving Yi, it won't be updated if Y is observed
        } else if( aux$obs.type=="M_paired" ) {
          for( i in 1:aux$numSamples ) { solve.Yi.paired(i) }
          solve.sigma2.paired() # since sigma2 is only needed for solving Yi, it won't be updated if Y is observed
        }
        # Solve A, if C is provided
        if( aux$P > 0 ) {
          solve.A()
        }
        # Solve Ri for all samples, if H is provided
        if( aux$L > 0 ) {
          for( i in 1:aux$numSamples ) { solve.Ri(i) }
        }
        # Calculate and store the tracking statistics
        if( aux$ite %% track_internval == 1 ) {
          track_changes( prev_params )
        }
      }
    },
    
    #######################
    # Function for extrapolation of sigma2
    #######################

    # Extrapolate the value of sigma2 for iter -> Inf.
    #  This function works well if (a) other model components have converged, and
    #  (b) sigma2 has been decreasing in the past iterations.
    # `use.ite`: the number of the most recent iterations that should be used for extrapolation
    extrapolate.sigma2 = function( use.ite=10 )
    {
      # extact the most recent sigma2 values
      x <- seq(aux$ite-use.ite+1,aux$ite)
      y <- tracking$sigma2[x]
      # calculate the inverse
      x_inv <- 1/x
      y_inv <- 1/y
      message(
        "Correlation between 1/t and 1/sigma2: ", cor(x_inv,y_inv) )
      # When x -> Inf, 1/x -> 0. Therefore, the extrapolated sigma2 at x -> Inf
      #  is the inverse of the intercept of the line 1/y ~ 1/x.
      fit <- lm( y_inv ~ x_inv )
      if( coef(summary(fit))[2,1] > 0 ) {
        message("Warning: sigma2 has been increasing. Extrapolation may not be suitable.")
      }
      if( coef(summary(fit))[1,1] < 0 ) {
        message("Warning: Extrapolation leads to a negative number. Sigma2 was left unchanged.")
      } else {
        params$sigma2 <<- 1/coef(summary(fit))[1,1]
        cat(params$sigma2,"\n")
        message("Recalculating Y...")
        if( aux$obs.type=="M" ) {
          for( i in 1:aux$numSamples ) { solve.Yi(i) }
        } else if( aux$obs.type=="M_paired" ) {
          for( i in 1:aux$numSamples ) { solve.Yi.paired(i) }
        }
      }
    },
    
    
    #######################
    # Function for imputing Y given raw counts, after model fitting is complete
    #
    # The effect of sample-specific distortions and cell-specific library sizes
    # will be removed before returning the imputed values.
    #######################
    
    # `logScale`: whether the imputed values should be on the logarithmic scale
    # `rowCentre`: whether the global gene-specific offsets should be removed from imputed values
    imputeY = function( logScale=T, rowCentre=T ) {
      if( aux$obs.type=="Y" ) {
        message("Warning: observation type is not raw counts.")
      }
      # Remove the effects of sample-specific factors from Yi, and concatenate
      Y_res <- NULL
      for( i in 1:aux$numSamples ) {
        if( rowCentre ) {
          Y_res <- cbind( Y_res, Yi_resZ(
            target$Yi[[i]],
            aux$QiDBi[[i]],
            params$si[[i]], params$o, params$oi[[i]] ) )
        } else {
          Y_res <- cbind( Y_res, Yi_resZ(
            target$Yi[[i]],
            aux$QiDBi[[i]],
            params$si[[i]], rep(0,aux$J), params$oi[[i]] ) )
        }
      }
      rownames(Y_res) <- aux$geneIDs
      colnames(Y_res) <- aux$cellIDs
      # Return the imputed values
      if( logScale ) {
        return( Y_res ) # Y_res is itself on the log-scale
      } else if( aux$obs.type=="M" ) { # If observation is simply the counts, return exp of Y
        return( exp(Y_res) )
      } else if( aux$obs.type=="M_paired") { # If observation is a pair of count matrices, return logistic of Y
        return( 1/(1+exp(-Y_res)) )
      } else {
        stop("Unrecognized mode.")
      }
    },
    
    #######################
    # Function to return the SVD projection of the integrated data
    #  After decomposing the GEDI ZDB projection into U %*% S %*% t(V), returns S %*% t(V)
    #######################
    
    # `K`: the target rank. If not-specified, it will be set to the number of
    #   latent variables of the model
    svd.projection = function( K = NA, randomSeed = NA ) {
      # Determine the target rank of the SVD decomposition
      if( is.na(K) ) { K = aux$K }
      if( K > aux$K ) { K = aux$K }
      # Concatenate the ZxDxBi matrices
      ZDB <- NULL
      for( i in 1:aux$numSamples ) {
        ZDB <- cbind( ZDB, aux$ZDBi[[i]] )
      }
      # Calculate the rSVD
      if( !is.na(randomSeed) ) { set.seed(randomSeed) } # To ensure reproducible results
      svd_ZDB <- rsvd( ZDB, k=K, nu=K, nv=K ) # use rsvd for fast SVD decomposition
      SVt <- eigenMatTcrossprod( diag(svd_ZDB$d,nrow=K), svd_ZDB$v )
      # set the row and column names
      rownames(SVt) <- paste0("LV",1:K)
      colnames(SVt) <- aux$cellIDs
      # return the matrix X
      return(SVt)
    },
    
    
    #######################
    # Function to return the GEDI ZDB projection (i.e., ZxDxB matrix)
    #######################
    
    getZDB = function() {
      # Concatenate the ZxDxBi matrices
      ZDB <- NULL
      for( i in 1:aux$numSamples ) {
        ZDB <- cbind( ZDB, aux$ZDBi[[i]] )
      }
      # set the row and column names, and return the result
      rownames(ZDB) <- aux$geneIDs
      colnames(ZDB) <- aux$cellIDs
      return( ZDB )
    },
    
    
    #######################
    # Function to return the GEDI A and ADB projection (i.e., AxDxB matrix)
    #######################
    
    getADB = function() {
      if( aux$P > 0 ) {
        # Concatenate the AxDxBi matrices
        ADB <- NULL
        for( i in 1:aux$numSamples ) {
          ADB <- cbind( ADB, params$Bi[[i]] )
        }
        ADB <- eigenMatProduct(
          aux$C.rotation %*% params$A %*% diag(params$D,nrow=aux$K),
          ADB )
        # set the row and column names, and return the result
        rownames(ADB) <- colnames(aux$inputC)
        colnames(ADB) <- aux$cellIDs
        return( ADB )
      } else {
        return( aux$C )
      }
    },
    
    getA = function() {
      if( aux$P > 0 ) {
        A <- aux$C.rotation %*% params$A
        # set the row and column names, and return the result
        rownames(A) <- colnames(aux$inputC)
        colnames(A) <- paste0("LV:",1:aux$K)
        return( A )
      } else {
        return( aux$C )
      }
    },
    
    getRiDB = function( i ) {
      # Concatenate the AxDxBi matrices
      if( aux$L > 0 ) {
        RiDB <- NULL
        for( i in 1:aux$numSamples ) {
          RiDB <- cbind( RiDB, params$Bi[[i]] )
        }
        RiDB <- eigenMatProduct(
          aux$H.rotation %*% params$Ri[[i]] %*% diag(params$D,nrow=aux$K),
          RiDB )
        # set the row and column names, and return the result
        rownames(RiDB) <- colnames(aux$inputH)
        colnames(RiDB) <- aux$cellIDs
        return( RiDB )
      } else {
        return( aux$H )
      }
    },
    
    getRi = function( i ) {
      if( aux$L > 0 ) {
        Ri <- aux$H.rotation %*% params$Ri[[i]]
        # set the row and column names, and return the result
        rownames(Ri) <- colnames(aux$inputH)
        colnames(Ri) <- paste0("LV:",1:aux$K)
        return( Ri )
      } else {
        return( aux$H )
      }
    },

    #######################
    # Function to plot convergence metrics
    #######################
    
    plotTracking = function() {
      l <- length(tracking$dZ)
      # Create the plot for convergence of Z, A, and o
      ZAo <- data.frame(
        Iteration = rep(1:l), RMSD=tracking$dZ, Parameter="Z" )
      ZAo <- rbind( ZAo, data.frame(
        Iteration = rep(1:l), RMSD=tracking$do, Parameter="o" ) ) 
      if(aux$P>0) {
        ZAo <- rbind( ZAo, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dA, Parameter="A" ) ) 
      }
      ZAo <- ZAo[!is.na(ZAo$RMSD),] # Remove undocumented iterations
      # Create the plot for convergence of sample-specific parameters
      Bi <- Qi <- oi <- si <- Ri <- NULL
      for( i in 1:aux$numSamples ) {
        Bi <- rbind( Bi, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dBi[[i]], Bi=aux$Samples[[i]] ) ) 
        Qi <- rbind( Qi, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dQi[[i]], Qi=aux$Samples[[i]] ) ) 
        oi <- rbind( oi, data.frame(
          Iteration = rep(1:l), RMSD=tracking$doi[[i]], oi=aux$Samples[[i]] ) ) 
        si <- rbind( si, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dsi[[i]], si=aux$Samples[[i]] ) ) 
        if(aux$L>0) {
          Ri <- rbind( Ri, data.frame(
            Iteration = rep(1:l), RMSD=tracking$dRi[[i]], Ri=aux$Samples[[i]] ) ) 
        }
      }
      # Remove undocumented iterations
      Bi <- Bi[!is.na(Bi$RMSD),]
      Qi <- Qi[!is.na(Qi$RMSD),]
      oi <- oi[!is.na(oi$RMSD),]
      si <- si[!is.na(si$RMSD),]
      if(aux$L>0) { Ri <- Ri[!is.na(Ri$RMSD),] }
      # Create the ggplot objects
      ZAo <- ggplot(ZAo,aes(x=Iteration,y=RMSD,color=Parameter))+
        geom_line()+geom_point()+theme_bw()+scale_y_log10()
      Bi <- ggplot(Bi,aes(x=Iteration,y=RMSD,color=Bi))+
        geom_line()+geom_point()+theme_bw()+scale_y_log10()
      Qi <- ggplot(Qi,aes(x=Iteration,y=RMSD,color=Qi))+
        geom_line()+geom_point()+theme_bw()+scale_y_log10()
      oi <- ggplot(oi,aes(x=Iteration,y=RMSD,color=oi))+
        geom_line()+geom_point()+theme_bw()+scale_y_log10()
      si <- ggplot(si,aes(x=Iteration,y=RMSD,color=si))+
        geom_line()+geom_point()+theme_bw()+scale_y_log10()
      if(aux$L>0){
        Ri <- ggplot(Ri,aes(x=Iteration,y=RMSD,color=Ri))+
          geom_line()+geom_point()+theme_bw()+scale_y_log10()
      }
      # Now create the tracking plot for sigma2
      sigma2 <- NULL
      if( !is.na(tracking$sigma2[1]) ) {
        sigma2 <-
          ggplot( data.frame(
            Iteration=1:length(tracking$sigma2),
            MSE=tracking$sigma2 ),
            aes(x=Iteration,y=MSE) ) +
          geom_line()+geom_point()+theme_bw()+scale_y_log10()
      }
      return(list(ZAo=ZAo,Bi=Bi,Qi=Qi,oi=oi,si=si,Ri=Ri,sigma2=sigma2))
    }
  )
)

# Loading libraries and scripts that are needed for core functionality
message("Loading required libraries and scripts...")
library(Rcpp)
library(quadprog) # required for solve.QP
library(ClusterR) # required for KMeans_arma
library(Matrix) # required for nearPD, and also handling sparse matrices
library(rsvd) # required for rsvd
library(corpcor) # required for fast.svd
sourceCpp("./gedi/Rcpp/Eigen.v6.cpp")


################ Auxiliary functions

################################################################
## Functions to read files
################################################################

#' Reading data from Enrichr
readEnrichr <- function(filename)
{
    tmp <- fread(filename,data.table=F,header = F)
    genes <- unique( unlist( strsplit( unlist(tmp[,-1]), "," ) ) )
    genes <- genes[!is.na(genes)]
    
    mx <- matrix( 0, nrow=length(genes),ncol=nrow(tmp) )
    colnames(mx) <- unlist(tmp[,1])
    rownames(mx) <- genes
    for( i in 1:ncol(mx) )
    {
        mx[,i] <- ( genes %in% unlist(strsplit(tmp[i,-1],",")) )*1
    }
    return(mx)
}

#' Reading Human Cell Atlas file
readHCL <- function(filename)
{
    tmp <- fread(filename,data.table=F)
    tmp <- tmp[ , which( tmp[2,]=="n" ) ]
    
    genes <- unique(unlist(tmp[-(1:3),]))
    genes <- genes[!is.na(genes)]
    
    mx <- matrix( 0, nrow=length(genes),ncol=ncol(tmp) )
    colnames(mx) <- unlist(tmp[1,])
    rownames(mx) <- genes
    for( i in 1:ncol(mx) )
    {
        mx[,i] <- ( genes %in% tmp[-(1:3),i] )*1
    }
    return(mx)
}

# Loading libraries and scripts that are needed for auxiliary functions
library(ggplot2)
library(scales)

# Function for plotting
plot_umap <- function(umap_mat,colour)
{
  umap_obj <- data.frame(
    UMAP_1=umap_mat[,1],
    UMAP_2=umap_mat[,2],
    var_plot=colour )
  
  umap_obj <- umap_obj[ sample.int(nrow(umap_obj)), ]
  if(is.numeric(colour))
  {
    umap_obj$var_plot <- umap_obj$var_plot - mean(umap_obj$var_plot)
    lim <- quantile(abs(umap_obj$var_plot),0.99)
    ggplot(
      umap_obj, aes_string( "UMAP_1","UMAP_2", colour="var_plot"))+
      geom_point(size=0.05)+
      theme_minimal()+
      scale_color_gradientn( limits=c(-lim,lim), colours=c("blue","light grey","red"), oob=squish )
  } else {
    ggplot(
      umap_obj, aes_string( "UMAP_1","UMAP_2", colour="var_plot"))+
      geom_point(size=0.05)+
      theme_minimal()+guides(colour=guide_legend(override.aes=list(size=3)))
    
  }
}


#' Function to calculate alignment score based on Butler et al 2018
#' This score quantifies how well any group of data sets is aligned
#' First, the datasets are randomly sub-sampled to have the same number of cells
#' as the smallest dataset. Then, a nearest-neighbor graph is constructed based on
#' the cell's embedding low-dimensional space. For every cell, it calculates
#' how many of its k nearest neighbors belong to the same data set and averages this over all cells to obtain x_hat.
#' N is the number of datasets
#' The formula is: 1  - ((x_hat -k/N) / (k - k/N ) )
#' @param X Cell's embedding
#' @param meta Metadata that contains a column that has the batches groups
#' @param per Percentage of cells to use for estimating K
#' @param var_use Variable name in meta that has the batches groups. Default to Sample
#' @param seed Seed for the sub-sampling per batch
#' @return Alignment score. 
alignment_score<- function(X, meta, per=0.01, var_use="Sample",seed=43){
    require(BiocNeighbors)
    ## Subsample cells to the smallest data set
    Nsmall<- sort(unique(table(meta[,var_use])) )[1]
    l<- lapply( unique(meta[,var_use]), function(sample){
        sample_meta<- meta[meta[,var_use] == sample,]
        set.seed(seed)
        return( sample_meta[sample(rownames(sample_meta), Nsmall),])
    })
    meta<- do.call(rbind, l)
    N<- length(unique(meta[,var_use])) ## Number of datasets
    X<- X[,rownames(meta)]
    ## Choosing k based on a percentage of cells, with minimum 10
    k<- max( round(per*N), 10)
    ## Get Nearest Neighbors
    NNs<- findKNN( t(X), k=k)
    ## For every cell, we then calculate how many of its k nearest-neighbors belong to the same data set and average this over all cells 
    x_hat<- mean(sapply( 1:nrow(meta), function(i) sum( meta[,var_use][NNs$index[i,]] ==  meta[i,var_use]) ))
    as<- 1- ( (x_hat - k/N) / (k - k/N ) )    
    return(as)    
}

pheatmap.colorsymmetric <- function(x,...)
{
    require(pheatmap)
    lim <- max(abs(x), na.rm=TRUE)
    if( min(x, na.rm=TRUE) < 0 ){
        lim_down<- -lim
        col_palette<- colorRampPalette(c("blue","white","red"))(256)
    }else{        
        lim_down<- 0
        col_palette<- colorRampPalette(c("white","red"))(256)
    }
    pheatmap(
        x, color = col_palette,
        breaks=seq(lim_down,lim,length.out=255), ... )
}


# res: result list from a previous gedi call
# lv: the LV that should be examined
# C: the C matrix that was provided to gedi
# top_k: the number of top genes to be shown
# bottom_k: the number of bottom genes to be shown
# C_max: the maximum number of C entries to be shown
plot_top_genes <- function( res, lv="LV1",C,  top_k=20, bottom_k=20, C_max=10 )
{
    ## Select the top and bottom genes based on the Z column for the specified LV
    G <- nrow(res$Z)
    genes <- order(res$Z[,lv],decreasing = T)[c(1:top_k,(G-bottom_k+1):G)]
    Y_normalized<- do.call( cbind, lapply(res$datasets, function(slide_data) slide_data$Yi ) )
    Y_select <- Y_normalized[genes,]
    
    C_select <- t(t(C)*res$A[,lv]) # First, weight the columns of C by their A coefficient
    C_range <- max(apply(C_select,2,function(x)max(x)-min(x))) # store the range of C after weighting. This is later useful for visualization
    C_select <- C_select[genes,res$A[,lv]>0] # Then, select the entries that correspond to non-zero A coefficients as well as the top/bottom genes
    ## Sort the columns of C by their contribution to the top/bottom genes. Also keep at most C_max columns
    C_select <- C_select[ , order(
        colSums(C_select[1:top_k,]),decreasing = T)[1:min(ncol(C_select),C_max)] ]
    
    ## Order the columns of the heatmap by how the cells/samples are projected on the LV
    col_order <- order( res$Vt[lv,] )
    Y_select <- Y_select[ , col_order ]
    
    ## Set the colors of the C matrix so that they reflect the contribution of each column
    C_select <- as.data.frame(C_select)
    annot_colors<- list()
    for( i in colnames(C_select) )
    {
        maxC <- max(C[,i]) # The maximum value in the corresponding column of the original C matrix
        minC <- min(C[,i]) # The minimum value in the corresponding column of the original C matrix
        index <- (maxC-minC)*res$A[i,lv]/C_range*255
        annot_colors[[i]]<- c(rgb(0.9,0.9,0.9),colorRampPalette(c(rgb(0.9,0.9,0.9),rgb(0,0,0)))(255)[index])
    }
    names(annot_colors)<- colnames(C_select)


    ## determine the range of expression values in the Y matrix, and return the final heatmap
    Y_range <- quantile(abs(c(as.matrix(Y_select))),0.99)
    return( pheatmap(
        as.matrix(Y_select),        
        cluster_rows = F,
        cluster_cols = F,
        breaks = seq(-Y_range,Y_range,length.out = 255),
        color = colorRampPalette(c("blue","white","red"))(256),
        annotation_row=C_select,
        annotation_colors=annot_colors,
        annotation_legend=FALSE,
        show_colnames= FALSE ) )
}


#' Function to return object for ggplot2 plotting
#' @param size_text Size of the text to use for plotting
#' @param leg_pos Position of legend
#' @return ggplot theme object
fun_theme_plot <- function(size_text=15,
                           leg_pos="right",
                           size_text_titles=NULL,
                           classic=FALSE,
                           line_size=0.85,
                           line_size_color="black",
                           axis_ticks_length=0.25,
                           axis_ticks_color="black"
                           ){
    require(ggplot2)
    if( is.null(size_text_titles) ){
        size_text_titles<- size_text + 2
    }
    if( classic ){
        gg.plots<- theme_classic()
    }else{
        gg.plots<- theme_bw()
    }
    gg.plots <- gg.plots +
        theme(axis.text.x=element_text(size=size_text, angle=45, hjust=1, color=line_size_color),
              axis.line=element_line(size=line_size, colour=line_size_color),
              axis.ticks=element_line(colour=axis_ticks_color, size=line_size),
              axis.ticks.length=unit(axis_ticks_length, "cm"),
              axis.title.x=element_text(size=size_text_titles),
              legend.position=leg_pos,
              plot.title=element_text(hjust=0.5, size=size_text),
              axis.text.y=element_text(size=size_text, color="black"),
              axis.title.y=element_text(size=size_text_titles) )
}

#' Function to plot embedding
#'@param obj Could be a dataframe that contains coordinates in the 2 first columns or a SCE object
#'@param var_plot Variable to plot
#'@param meta Additional metadata with variables to plot
#' @param center Whether to center continuous variables (remove mean)
#' @return plot 
plot_dimred <- function(obj,                        
                        var_plot=NULL,
                        meta=NULL,
                        dimred_use="umap",
                        center=FALSE,
                        size_text=20,
                        size_dot=0.5,
                        leg_pos="right",
                        cols_use=NULL,
                        title_use="",                         
                        size_guide_legend=7,                        
                        use_ggrastr=FALSE
                        ){
    require(ggplot2)
    require(RColorBrewer)
    require(ggrastr)
    require(viridis)

    ## Checking if obj is a SCE object
    if( class(sce) == "SingleCellExperiment" ){
        dim_obj<- cbind(  reducedDim(sce, dimred_use), data.frame(colData(sce) , check.names=FALSE))
    }    
    col1<- colnames(dim_obj)[1]
    col2<- colnames(dim_obj)[2]
    if( is.null(meta) ){
        meta<- dim_obj[,-c(1:2), drop=FALSE]
    }
    if( !is.null(var_plot) ){
        colour<- meta[,var_plot]
        if(is.numeric(colour)){
            if( grepl("\\-", var_plot) ){ ## Changing - into . for the variable name
                var_plot<- gsub("\\-", "neg", var_plot)                
            }
            if( grepl("\\+", var_plot) ){ ## Changing - into . for the variable name
                var_plot<- gsub("\\+", "pos", var_plot)                
            }
            if( center ){
                dim_obj[,var_plot] <- colour - mean(colour)
                lim <- quantile(abs(dim_obj[,var_plot]),0.9999)
            }else{
                dim_obj[,var_plot] <- colour            
            }
            ggp<- ggplot(dim_obj, aes_string( col1, col2, colour=var_plot))
            if( use_ggrastr ){
                ggp<- ggp + rasterise(geom_point(size=size_dot), dpi = 72)
            }else{
                ggp<- ggp + geom_point(size=size_dot)
            }
            if( center) {
                ggp<- ggp +
                    scale_color_gradientn( limits=c(-lim,lim), colours=c("blue","#d8d8d8","red") )
            }else{
                ggp<- ggp +
                    scale_color_viridis()           
            }
            ggp<- ggp +
                fun_theme_plot(size_text=size_text, leg_pos=leg_pos)+            
                labs(color=var_plot)
            
        } else {
            dim_obj[,var_plot] <- colour
            ggp<- ggplot(dim_obj, aes_string(col1,col2, color=var_plot))
            if( use_ggrastr ){
                ggp<- ggp + rasterise(geom_point(size=size_dot), dpi = 72)
            }else{
                ggp<- ggp + geom_point(size=size_dot)
            }
            ggp<- ggp +
                fun_theme_plot(size_text=size_text, leg_pos=leg_pos)+        
                guides(colour = guide_legend(override.aes = list(size=size_guide_legend)))+
                labs(color=var_plot)
            if( !is.null(cols_use) ){
                ggp<- ggp +
                    scale_color_manual(values=c(cols_use))
            }
        }
    }else{
        ggp<- ggplot(dim_obj, aes_string( col1, col2)) 
        if( use_ggrastr ){
            ggp<- ggp + rasterise(geom_point(size=size_dot), dpi = 72)
        }else{
            ggp<- ggp + geom_point(size=size_dot)
        }
        ggp<- ggp +
            fun_theme_plot(size_text=size_text, leg_pos=leg_pos)                        
    }
    return(ggp)    
}

summary_markers<- function(markers_sce){
    lis_markers<- lapply(names(markers_sce), function(cluster_look) {
        temp_df<- markers_sce[[cluster_look]]
        temp_df$Gene<- rownames(temp_df)
        temp_df$cluster<- cluster_look
        temp_df$min_lfc<- apply(temp_df[,grep("logFC.", colnames(temp_df),value=TRUE)], 1, min)
        temp_df$neglog_FDR<- -log10(temp_df$FDR)
        temp_df<- temp_df[sort(rownames(temp_df)),]
        data.frame(temp_df[,c("Gene", "cluster", "p.value", "FDR", "neglog_FDR", "min_lfc")])
    })
    names(lis_markers)<- names(markers_sce)    
    return(lis_markers)    
}
