## usethis namespace: start
#' @import RcppEigen
#' @importFrom Rcpp sourceCpp
#' @useDynLib GEDI
## usethis namespace: end
NULL

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
      colData = NULL, # An optional data frame for storing the column metadata (cell metadata). It will not be used in the modeling
      C = NULL, # The gene-level biological prior. If NULL, it means that there is no prior for Z
      H = NULL, # The sample-level prior for sources of variation. If NULL, there will be no prior for Qi. Must be a matrix with column names corresponding to sample names
      K = 10, # The number of latent variables
      mode = "Bl2", # Three values are allowed: "BL2" [L2 norm of the entire B matrix is fixed], or "Bsphere" [L2 norms of B columns are fixed]
      adjustD = T, # Whether D should be adjusted based on the L2 norm of B rows (TRUE), or just use the default (FALSE)
      orthoZ = T, # Whether the columns of Z should be orthogonal
      Z_shrinkage = 1, # The shrinkage multiplier for Z
      A_shrinkage = 1, # The shrinkage multiplier for A
      Qi_shrinkage = 1, # The shrinkage multiplier for Qi
      Rk_shrinkage = 1, # The shrinkage multiplier for Rk
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

      aux$K <<- K
      aux$diag_K <<- diag(K)

      # store additional auxiliary information

      aux$BcolL2 <<- 1 # The expected L2 norm of each column of B in the Bsphere mode.
      aux$BrowL2 <<- 1 # The expected L2 norm of each row of B in the Bnorm mode
      aux$update_B <<- T

      aux$Samples <<- unique( Samples ) # The list of unique samples; all lists will have the same order of samples
      aux$numSamples <<- length(aux$Samples)
      aux$cells <<- vector("list",aux$numSamples) # The indexes of cells originating from each sample
      aux$Ni <<- rep(NA,aux$numSamples) # The number of cells originating from each sample
      aux$Ni_vec <<- vector("list",aux$numSamples) # A list containing a vector of 1's for each sample
      for( i in 1:aux$numSamples ) {
        aux$cells[[i]] <<- which( Samples==aux$Samples[i] )
        aux$Ni[i] <<- length(aux$cells[[i]])
        aux$Ni_vec[[i]] <<- rep(1,aux$Ni[i])
      }
      aux$cells.all <<- do.call(c,aux$cells) # store the order of the cells for all samples relative to the input data
      aux$cellIDs <<- aux$cellIDs[aux$cells.all] # reorder the cell IDs to match the order in which the matrices are stored
      aux$colData <<- colData
      if( !is.null(aux$colData) ) {
        aux$colData <<- aux$colData[aux$cells.all,,drop=F]
      }

      # If provided, reorder the columns of the H matrix to follow the sample order,
      #  and prepare the H matrix for fitting GEDI
      if( !is.null(H) ) {
        aux$inputH <<- as.matrix(H) # To ensure `matrix` type
        aux$inputH <<- aux$inputH[ ,match(aux$Samples,colnames(aux$inputH)), drop=F ] # reorder columns
        if( ncol(aux$inputH) != aux$numSamples ) { stop("Incompatible matrices.") } # check if all samples are present in H
        minHL2 <- min(rowSums(aux$inputH^2)) # calculate the L2 of the smallest row of H
        svdH <- svd(t(aux$inputH)) # Perform SVD on t(H) to obtain uncorrelated components
        aux$L <<- sum( cumsum(rev(svdH$d^2)) > minHL2*0.1 ) # only consider the SVD components that have at least 10% of the variance of the smallest column of H
        aux$H <<- t(svdH$u[,1:aux$L]) # The H matrix used by GEDI is in fact the U matrix of SVD(H)
        aux$H.rotation <<- t( svdH$v[,1:aux$L] %*% diag(1/svdH$d[1:aux$L],nrow=aux$L) ) # Rotation matrix, so that aux$H.rotation %*% aux$inputH == aux$H
      } else {
        aux$H <<- matrix(nrow=0,ncol=0)
        aux$L <<- 0
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

      aux$update_priors <<- F
      hyperparams$scaling <<- 1
      hyperparams$S_o <<- 1/aux$N/o_shrinkage # The prior relative variance for o
      hyperparams$S_si <<- 1/aux$J/si_shrinkage # The prior relative variance for si
      hyperparams$S_Z <<- 1/Z_shrinkage # The prior relative variance for Z
      hyperparams$S_A <<- 1/A_shrinkage # The prior relative variance for A
      hyperparams$S_R <<- 1/Rk_shrinkage # The prior relative variance for Rk
      hyperparams$S_oi <<- 1/aux$Ni/oi_shrinkage # The prior relative variance for oi
      hyperparams$S_oi_mean <<- aux$numSamples/aux$N/oi_shrinkage # The prior relative variance for oi
      hyperparams$S_Qi <<- aux$N/aux$Ni/Qi_shrinkage # The prior relative variance for Qi
      hyperparams$S_Qi_mean <<- aux$numSamples/Qi_shrinkage # The mean prior relative variance for Qi; this is in fact used for fitting Rk


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
      aux$oi_hat <<- vector("list",aux$numSamples)
      params$Bi <<- vector("list",aux$numSamples)
      params$Qi <<- vector("list",aux$numSamples)
      aux$Qi_hat <<- vector("list",aux$numSamples)

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
        aux$oi_hat[[i]] <<- rep(0,aux$J) # Initialize oi with zero
        params$Bi[[i]] <<- matrix(0,nrow=K,ncol=aux$Ni[i]) # Initialize Bi with zero
        params$Qi[[i]] <<- matrix(0,nrow=aux$J,ncol=K) # Initialize Qi with zero
        aux$Qi_hat[[i]] <<- matrix(0,nrow=aux$J,ncol=K) # Initialize Qi_hat with zero
        aux$ZDBi[[i]] <<- matrix(0,nrow=aux$J,ncol=aux$Ni[i]) # ZDBi stores Z x Bi product
        aux$QiDBi[[i]] <<- matrix(0,nrow=aux$J,ncol=aux$Ni[i]) # QiDBi stores Qi x Bi product
      }

      # Initialize factor-specific matrices
      params$Rk <<- vector("list",aux$K)
      for( i in 1:aux$K ) {
        if( is.null(H) ) {
          params$Rk[[i]] <<- matrix(nrow=0,ncol=0)
        } else {
          params$Rk[[i]] <<- matrix(0,nrow=aux$J,ncol=aux$L)
        }
      }
      # Ro predicts the sample-specific oi from sample variables
      if( is.null(H) ) {
        params$Ro <<- matrix(nrow=0,ncol=0)
      } else {
        params$Ro <<- matrix(0,nrow=aux$J,ncol=aux$L)
      }

      # Initialize tracking vectors
      tracking$dZ <<- NA # This vector tracks the convergence of Z
      tracking$dA <<- NA # This vector tracks the convergence of A
      tracking$do <<- NA # This vector tracks the convergence of o
      tracking$dRo <<- NA # This vector tracks the convergence of Ro

      tracking$sigma2 <<- NA # This vector tracks the value of sigma2 across iterations

      tracking$dsi <<- vector("list",aux$numSamples) # This list tracks the convergence of si for each sample
      tracking$doi <<- vector("list",aux$numSamples) # This list tracks the convergence of oi for each sample
      tracking$dBi <<- vector("list",aux$numSamples) # This list tracks the convergence of Bi for each sample
      tracking$dQi <<- vector("list",aux$numSamples) # This list tracks the convergence of Qi for each sample
      tracking$dRk <<- vector("list",aux$numSamples) # This list tracks the convergence of Rk for each sample
      for( i in 1:aux$numSamples ) {
        tracking$dsi[[i]] <<- NA
        tracking$doi[[i]] <<- NA
        tracking$dBi[[i]] <<- NA
        tracking$dQi[[i]] <<- NA
      }
      for( k in 1:aux$K ) {
        tracking$dRk[[k]] <<- NA
      }
    },

    # This function locks the prior scaling factor to the current value
    #  If another value is provided to this function, the scaling factor will be
    #  changed to that value
    lock.priors = function( scaling_factor = NA ) {
      aux$update_priors <<- F
      if( !is.na(scaling_factor) ) {
        hyperparams$scaling <<- scaling_factor
      }
    },

    # This function releases the prior scaling factor so that it's iteratively
    #  updated during optimization, starting from the current value. If another
    #  value is provided to this function, the scaling factor will be initialized
    #  to that value
    release.priors = function( scaling_factor = NA ) {
      aux$update_priors <<- T
      if( !is.na(scaling_factor) ) {
        hyperparams$scaling <<- scaling_factor
      }
    },

    # This function locks the Bi matrices to their current values so that they
    #  are not updated anymore during optimization. If a new B is provided, it
    #  is used to set the Bi matrices. Note that the columns in the provided B
    #  matrix must be in the same order as the original Y/M matrices provided
    #  to the setup function
    lock.B = function( B = NULL ) {
      aux$update_B <<- F
      if( !is.null(B) ) {
        for( i in 1:aux$numSamples ) {
          params$Bi[[i]] <<-  B[,aux$cells[[i]],drop=F]
        }
      }
    },

    # This function releases the Bi matrices so that they can be updated during
    #  optimization. If a new B is provided, it is used to reinitialize the Bi
    #  matrices. Note that the columns in the provided B matrix must be in the
    #  same order as the original Y/M matrices provided to the setup function
    release.B = function( B = NULL ) {
      aux$update_B <<- T
      if( !is.null(B) ) {
        for( i in 1:aux$numSamples ) {
          params$Bi[[i]] <<-  B[,aux$cells[[i]],drop=F]
        }
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
      if( aux$update_B ) {
        params$Bi[[i]] <<- solveBi(
          target$Yi[[i]],
          params$D,
          params$Z, params$Qi[[i]],
          aux$diag_K,
          params$si[[i]], params$o, params$oi[[i]] )
      }

      #return(1)
    },


    ######
    # Normalized Bi for all samples, depending on the mode of the GEDI model
    normalize.B = function() {
      if( aux$mode == "Bsphere" & aux$update_B ) {
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

      if( aux$mode == "Bl2" & aux$update_B ) { # In the Bl2 mode, the L2 norm of each row of B is fixed.
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
      lambda <- 1/(hyperparams$S_Z*hyperparams$scaling)
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

      lambda <- 1/(hyperparams$S_Z*hyperparams$scaling)

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
      lambda <- 1/hyperparams$S_A # No need to multiply by hyperparams$scaling
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
      lambda <- 1/(hyperparams$S_Qi[i]*hyperparams$scaling)
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
          aux$Qi_hat[[i]],
          aux$diag_K, lambda )
      }
      #return(1)
    },


    ######
    # Solves for Rk, the matrix that connects the prior H to the sample-specific metagenes Qi
    # `k`: the index of the factor for which Rk is being updated
    solve.Rk = function( k ) {
      # Build the Qk and Hp matrices required for solving Rk
      Qk <- matrix(0,nrow=aux$J,ncol=aux$numSamples)
      Hp <- aux$H
      for( i in 1:aux$numSamples ) {
        Qk[,i] <- params$Qi[[i]][,k]/sqrt(hyperparams$S_Qi[i])
        Hp[,i] <- Hp[,i]/sqrt(hyperparams$S_Qi[i])
      }
      # Solve Rk
      lambda <- 1/hyperparams$S_Qi_mean/hyperparams$S_R # No need to multiply by hyperparams$scaling
      params$Rk[[k]] <<- eigenMatProduct(
        eigenMatTcrossprod(Qk,Hp),
        solve( eigenMatTcrossprod(Hp,Hp) + lambda*diag(aux$L) ) ) # The code can be optimized by pre-calculating the second part (solve function), which remains constant as long as S_Qi[i] values are constant
      # Update Qi_hat for each sample, which will be used in the next round of optimizing Qi
      Qk_hat <- eigenMatProduct(params$Rk[[k]],aux$H)
      for( i in 1:aux$numSamples ) {
        aux$Qi_hat[[i]][,k] <<- Qk_hat[,i]
      }
      #return(1)
    },


    ######
    # Solve for oi, the sample-specific offset: Yi ~ o + oi + (Z+Qi) x Bi
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZDBi and QiDBi to be up-to-date
    solve.oi = function( i )
    {
      lambda <- 1/(hyperparams$S_oi[i]*hyperparams$scaling)
      if( aux$L==0 ) {
        params$oi[[i]] <<- solveOi_noH(
          target$Yi[[i]],
          aux$ZDBi[[i]], aux$QiDBi[[i]],
          aux$Ni_vec[[i]],
          params$si[[i]], params$o,
          aux$Ni[[i]], lambda )
      } else {
        params$oi[[i]] <<- solveOi_wH(
          target$Yi[[i]],
          aux$ZDBi[[i]], aux$QiDBi[[i]],
          aux$Ni_vec[[i]],
          params$si[[i]], params$o,
          aux$oi_hat[[i]],
          aux$Ni[[i]], lambda )

      }
      #return(1)
    },

    ######
    # Solves for Ro, the matrix that connects the prior H to the sample-specific metagenes oi
    solve.Ro = function() {
      # Build the O and Hp matrices required for solving Ro
      O <- matrix(0,nrow=aux$J,ncol=aux$numSamples)
      Hp <- aux$H
      for( i in 1:aux$numSamples ) {
        O[,i] <- params$oi[[i]]/sqrt(hyperparams$S_oi[i])
        Hp[,i] <- Hp[,i]/sqrt(hyperparams$S_oi[i])
      }
      # Solve Ro
      lambda <- 1/hyperparams$S_oi_mean/hyperparams$S_R # No need to multiply by hyperparams$scaling
      params$Ro <<- eigenMatProduct(
        eigenMatTcrossprod(O,Hp),
        solve( eigenMatTcrossprod(Hp,Hp) + lambda*diag(aux$L) ) ) # The code can be optimized by pre-calculating the second part (solve function), which remains constant as long as S_Qi[i] values are constant
      # Update oi_hat for each sample, which will be used in the next round of optimizing oi
      O_hat <- eigenMatProduct(params$Ro,aux$H)
      for( i in 1:aux$numSamples ) {
        aux$oi_hat[[i]] <<- O_hat[,i]
      }
      #return(1)
    },


    ######
    # Solves for si, the cell-specific library size
    # `i`: the index of the sample for which oi is being updated
    # Note that this function requires ZDBi and QiDBi to be up-to-date
    solve.si = function( i )
    {
      # lambda <- 1/(hyperparams$S_si*hyperparams$scaling)
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
      # lambda <- 1/(hyperparams$S_o*hyperparams$scaling)
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
      N1 <- aux$J*aux$N # Y
      # The denominator of scaling factor
      N2 <-
        aux$J*aux$numSamples + # oi
        aux$numSamples*aux$J*aux$K + # Qi
        (aux$K+1)*aux$J*aux$L + # Rk and Ro
        aux$K*(aux$J+aux$P) # Z and A
      # For the numerator, first calculate the weighted sum of all coefficients and residuals
      S1 <- 0
      S2 <- 0
      if( aux$P==0 ) # matrix C is not provided
      {
        S2 <- S2 + matL2_noPrior(params$Z) / hyperparams$S_Z
      } else { # matrix C is provided
        S2 <- S2 +
          matL2_wPrior(params$Z,aux$C,params$A) / hyperparams$S_Z +
          matL2_noPrior(params$A) / hyperparams$S_Z / hyperparams$S_A
      }
      for( i in 1:aux$numSamples ) # Add the sample-specific coefficients
      {
        if( aux$L==0 ) # matrix H is not provided
        {
          S2 <- S2 +
            vecL2_noPrior(params$oi[[i]]) / hyperparams$S_oi[i] +
            matL2_noPrior(params$Qi[[i]]) / hyperparams$S_Qi[i]
        } else { # matrix H is provided
          S2 <- S2 +
            vecL2_noPrior(params$oi[[i]]-aux$oi_hat[[i]]) / hyperparams$S_oi[i] +
            matL2_noPrior(params$Qi[[i]]-aux$Qi_hat[[i]]) / hyperparams$S_Qi[i]
        }
        # Calculate the residual of Yi, and add it to S
        if( aux$obs.type=="Y" ) {
          S1 <- S1 + Yi_SSE_fixed(
            target$Yi[[i]],
            aux$ZDBi[[i]], aux$QiDBi[[i]],
            params$si[[i]], params$o, params$oi[[i]] )
        } else if( aux$obs.type=="M" ) {
          S1 <- S1 + Yi_SSE_M(
            target$Yi[[i]],
            aux$ZDBi[[i]], aux$QiDBi[[i]],
            params$si[[i]], params$o, params$oi[[i]],
            params$sigma2 )
        } else if( aux$obs.type=="M_paired" ) {
          S1 <- S1 + Yi_SSE_M_paired(
            target$Yi[[i]],
            target$M1i[[i]], target$M2i[[i]],
            aux$ZDBi[[i]], aux$QiDBi[[i]],
            params$si[[i]], params$o, params$oi[[i]],
            params$sigma2 )
        } else { stop("Unrecognized mode.") }
      }
      # add the L2 norm of Rk's and Ro
      if( aux$L > 0 ) {
        for( k in 1:aux$K ) {
          S2 <- S2 + matL2_noPrior(params$Rk[[k]]) / hyperparams$S_Qi_mean / hyperparams$S_R
        }
        S2 <- S2 + matL2_noPrior(params$Ro) / hyperparams$S_oi_mean / hyperparams$S_R
      }
      # Normalize sigma2
      if( aux$update_priors ) {
        params$sigma2 <<- S1/N1
        hyperparams$scaling <<- S2/N2/params$sigma2+1 # The law of total variance
      } else {
        params$sigma2 <<- (S1+S2)/(N1+N2)
      }

      tracking$sigma2[aux$ite] <<- params$sigma2
      cat(params$sigma2,"\n")
      cat(hyperparams$scaling,"\n")

      # Update the hyperpriors for o and si
      hyperparams$o_0 <<- rep(mean(params$o),aux$J)
      hyperparams$S_o <<- vecL2_wPrior(params$o,hyperparams$o_0) / aux$J / params$sigma2
      tmp <- 0
      for( i in 1:aux$numSamples ) {
        tmp <- tmp + vecL2_wPrior(params$si[[i]],hyperparams$si_0[[i]])
      }
      hyperparams$S_si <<- tmp / aux$N / params$sigma2

      cat("Mean(o):",hyperparams$o_0[1],"; Var(o):",hyperparams$S_o,"; Var(si):",hyperparams$S_si,"\n")
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
        svdY <- rsvd::rsvd( Yp, k = aux$K, nu = aux$K, nv = aux$K )
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
      }
      if(aux$L>0) {
        for( k in 1:aux$K ) {
          tracking$dRk[[k]][aux$ite] <<- matRMSD(params$Rk[[k]],ref$Rk[[k]])
        }
        tracking$dRo[aux$ite] <<- matRMSD(params$Ro,ref$Ro)
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
        if( aux$ite %% track_internval == 1 | track_internval==1 ) {
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
          solve.sigma2() # since sigma2 is only needed for solving Yi, it won't be updated if Y is observed
        }
        # Solve A, if C is provided
        if( aux$P > 0 ) {
          solve.A()
        }
        # Solve Rk for all factors, if H is provided
        if( aux$L > 0 ) {
          for( k in 1:aux$K ) { solve.Rk(k) }
          solve.Ro()
        }
        # Calculate and store the tracking statistics
        if( aux$ite %% track_internval == 1 | track_internval==1 ) {
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
      Bi <- Qi <- oi <- si <- Rk <- NULL
      for( i in 1:aux$numSamples ) {
        Bi <- rbind( Bi, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dBi[[i]], Bi=aux$Samples[[i]] ) )
        Qi <- rbind( Qi, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dQi[[i]], Qi=aux$Samples[[i]] ) )
        oi <- rbind( oi, data.frame(
          Iteration = rep(1:l), RMSD=tracking$doi[[i]], oi=aux$Samples[[i]] ) )
        si <- rbind( si, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dsi[[i]], si=aux$Samples[[i]] ) )
      }
      if(aux$L>0) {
        for( k in 1:aux$K ) {
          Rk <- rbind( Rk, data.frame(
            Iteration = rep(1:l), RMSD=tracking$dRk[[k]], Rk=paste0("LV:",k) ) )
        }
        Rk <- rbind( Rk, data.frame(
          Iteration = rep(1:l), RMSD=tracking$dRo, Rk="o" ) )
      }

      # Remove undocumented iterations
      Bi <- Bi[!is.na(Bi$RMSD),]
      Qi <- Qi[!is.na(Qi$RMSD),]
      oi <- oi[!is.na(oi$RMSD),]
      si <- si[!is.na(si$RMSD),]
      if(aux$L>0) { Rk <- Rk[!is.na(Rk$RMSD),] }
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
        Rk <- ggplot(Rk,aes(x=Iteration,y=RMSD,color=Rk))+
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
      return(list(ZAo=ZAo,Bi=Bi,Qi=Qi,oi=oi,si=si,Rk=Rk,sigma2=sigma2))
    }
  )
)

# Loading libraries and scripts that are needed for core functionality
#message("Loading required libraries and scripts...")
#library(Rcpp)
#library(quadprog) # required for quantile
#library(ClusterR) # required for KMeans_arma
#library(Matrix) # required for nearPD, and also handling sparse matrices
#library(rsvd) # required for rsvd
#library(corpcor) # required for fast.svd


################ Auxiliary functions

# Loading libraries and scripts that are needed for auxiliary functions
library(ggplot2)
library(scales)


#' Retrieve stored colData
#'
#' Retrieve stored colData from GEDI object.
#'
#' @param object GEDI object
#'
#' @return colData
#' @export
#'
colData.gedi <- function( object ) {
  return(object$aux$colData)
}

#' Impute Y
#'
#' Impute Y given raw counts, after model fitting is complete
#' The effect of sample-specific distortions and cell-specific library sizes
#' will be removed before returning the imputed values.
#'
#' @param object GEDI object
#' @param logScale Whether the imputed values should be on the logarithmic scale
#' @param rowCentre Whether the global gene-specific offsets should be removed from imputed values
#'
#' @return Y imputed gene matrix
#' @export
#'
getY.gedi <- function( object, logScale=T, rowCentre=T ) {
  if( object$aux$obs.type=="Y" ) {
    message("Warning: observation type is not raw counts.")
  }
  # Remove the effects of sample-specific factors from Yi, and concatenate
  Y_res <- NULL
  for( i in 1:object$aux$numSamples ) {
    if( rowCentre ) {
      Y_res <- cbind( Y_res, Yi_resZ(
        object$target$Yi[[i]],
        object$aux$QiDBi[[i]],
        object$params$si[[i]], object$params$o, object$params$oi[[i]] ) )
    } else {
      Y_res <- cbind( Y_res, Yi_resZ(
        object$target$Yi[[i]],
        object$aux$QiDBi[[i]],
        object$params$si[[i]], rep(0,object$aux$J), object$params$oi[[i]] ) )
    }
  }
  rownames(Y_res) <- object$aux$geneIDs
  colnames(Y_res) <- object$aux$cellIDs
  # Return the imputed values
  if( logScale ) {
    return( Y_res ) # Y_res is itself on the log-scale
  } else if( object$aux$obs.type=="M" ) { # If observation is simply the counts, return exp of Y
    return( exp(Y_res) )
  } else if( object$aux$obs.type=="M_paired") { # If observation is a pair of count matrices, return logistic of Y
    return( 1/(1+exp(-Y_res)) )
  } else {
    stop("Unrecognized GEDI object mode.")
  }
}


#' Variance of imputed Y
#'
#' Return the variance of the imputed Y (log-scale)
#'
#' @param object GEDI object
#'
#' @return Matrix with variance of imputed Y
#' @export
#'
getY.var.gedi <- function( object ) {
  if( object$aux$obs.type=="Y" ) {
    message("Warning: observation type is not raw counts. Returning NULL.")
    return( NULL )
  }

  Y_var <- NULL
  for( i in 1:object$aux$numSamples ) {
    if( object$aux$obs.type=="M" ) { # If observation is simply the counts, return exp of Y
      Y_var <- cbind(
        Y_var,
        Yi_var( object$target$Yi[[i]], object$params$sigma2 ) )
    } else if( object$aux$obs.type=="M_paired") { # If observation is a pair of count matrices, return logistic of Y
      Y_var <- cbind(
        Y_var,
        Yi_var_paired(
          object$target$Yi[[i]],
          object$target$M1i[[i]], object$target$M2i[[i]],
          object$params$sigma2 ) )
    } else {
      stop("Unrecognized GEDI object mode.")
    }
  }
  rownames(Y_var) <- object$aux$geneIDs
  colnames(Y_var) <- object$aux$cellIDs
  return( Y_var )
}

#' Return GEDI DB
#'
#' Return the GEDI DB projection
#'
#' @param object GEDI object
#'
#' @return DxB matrix
#' @export
#'
getDB.gedi <- function( object ) {
  DB <- NULL
  for( i in 1:object$aux$numSamples ) {
    DB <- cbind( DB, object$params$Bi[[i]] )
  }
  DB <- eigenMatProduct( diag(object$params$D,nrow=object$aux$K), DB )
  rownames(DB) <- paste0("LV",1:object$aux$K)
  colnames(DB) <- object$aux$cellIDs
  return( DB )
}

#' Return GEDI ZDB
#'
#' Return the GEDI ZDB projection
#'
#' @param object GEDI object
#'
#' @return ZxDxB matrix
#' @export
#'
getZDB.gedi <- function( object ) {
  # Concatenate the ZxDxBi matrices
  ZDB <- NULL
  for( i in 1:object$aux$numSamples ) {
    ZDB <- cbind( ZDB, object$aux$ZDBi[[i]] )
  }
  # set the row and column names, and return the result
  rownames(ZDB) <- object$aux$geneIDs
  colnames(ZDB) <- object$aux$cellIDs
  return( ZDB )
}


#######################
# Functions to return the GEDI A and ADB projection (i.e., AxDxB matrix)
#######################

#' Return GEDI ADB
#'
#' Return the GEDI ADB projection
#'
#' @param object GEDI object
#'
#' @return AxDxB matrix
#' @export
#'
getADB.gedi <- function( object ) {
  if( object$aux$P > 0 ) {
    # Concatenate the AxDxBi matrices
    ADB <- NULL
    for( i in 1:object$aux$numSamples ) {
      ADB <- cbind( ADB, object$params$Bi[[i]] )
    }
    ADB <- eigenMatProduct(
      object$aux$C.rotation %*% object$params$A %*%
        diag(object$params$D,nrow=object$aux$K),
      ADB )
    # set the row and column names, and return the result
    rownames(ADB) <- colnames(object$aux$inputC)
    colnames(ADB) <- object$aux$cellIDs
    return( ADB )
  } else {
    return( object$aux$C )
  }
}

#' Return GEDI A
#'
#' Return the GEDI A projection
#'
#' @param object GEDI object
#'
#' @return A matrix
#' @export
#'
getA.gedi <- function( object ) {
  if( object$aux$P > 0 ) {
    A <- object$aux$C.rotation %*% object$params$A
    # set the row and column names, and return the result
    rownames(A) <- colnames(object$aux$inputC)
    colnames(A) <- paste0("LV:",1:object$aux$K)
    return( A )
  } else {
    return( object$aux$C )
  }
}


#' SVD projection of integrated data
#'
#' After decomposing the GEDI ZDB projection into U %*% S %*% t(V), returns S %*% t(V)
#'
#' @param object GEDI object
#'
#' @return list with SVD results ( d, u, v)
#' @export
#'
svd.gedi <- function( object ) {
  # Perform SVD of Z
  svd.Z <- svd(
    object$params$Z,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of DB
  svd.DB <- svd(
    getDB.gedi(object),
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of the middle SVD matrices
  svd.middle <- svd(
    diag(svd.Z$d,nrow=object$aux$K) %*% t(svd.Z$v) %*%
      svd.DB$u %*% diag(svd.DB$d,nrow=object$aux$K),
    nu=object$aux$K, nv=object$aux$K )
  # Calculate the SVD of ZDB
  u <- svd.Z$u %*% svd.middle$u
  v <- svd.DB$v %*% svd.middle$v
  d <- svd.middle$d
  # Set the dimension names
  names(d) <- colnames(u) <- colnames(v) <- paste0("LV",1:object$aux$K)
  rownames(u) <- object$aux$geneIDs
  rownames(v) <- object$aux$cellIDs

  return( list(d=d,u=u,v=v) )
}


#######################
# Functions to return the sample-variable effects on Qi and oi
#######################

#' Sample-variable effects on Qi
#'
#' Return the sample-variable effects on Qi
#'
#' @param object GEDI object
#' @param contrast a vector representing the sample-variable contrast to use for extracting effect on Qi; must be the same size as aux$L
#'
#' @return Matrix R
#' @export
#'
getDiffQ.gedi <- function( object, contrast ) {
  R <- NULL
  for( k in 1:object$aux$K ) {
    R <- cbind( R, object$params$Rk[[k]] %*% object$aux$H.rotation %*% contrast )
  }
  rownames(R) <- object$aux$geneIDs
  colnames(R) <- paste0("LV:",1:object$aux$K)
  return(R)
}

#' Sample-variable effects on oi
#'
#' Return the sample-variable effects on oi
#'
#' @param object GEDI object
#' @param contrast a vector representing the sample-variable contrast to use for extracting effect on oi; must be the same size as aux$L
#'
#' @return Matrix R
#' @export
getDiffO.gedi <- function( object, contrast ) {
  R <- object$params$Ro %*% object$aux$H.rotation %*% contrast
  rownames(R) <- object$aux$geneIDs
  colnames(R) <- "o"
  return( R )
}

#' Differential expression
#'
#'
#' @param object GEDI object
#' @param contrast a vector representing the sample-variable contrast to use for extracting effect on Qi; must be the same size as aux$L
#' @param include_O Logical, indicating whether the effect of contrast on the global offset O
#'
#' @return RDB matrix
#' @export
#'
getDiffExp.gedi <- function( object, contrast, include_O=F ) {
  RDB <- eigenMatProduct(
    getDiffQ.gedi( object, contrast ),
    getDB.gedi( object ) )
  if( include_O ) {
    RDB <- RDB + c( getDiffO.gedi( object, contrast ) )
  }
  rownames(RDB) <- object$aux$geneIDs
  colnames(RDB) <- object$aux$cellIDs
  return(RDB)
}

#' SVD for vector fields
#'
#' Calculate svd for cells and their vector fields.
#' @param object GEDI object
#' @param start.cond contrast vector for start condition
#' @param end.cond contrast vector for end condition
#' @param scale_cond_vector  a scaling factor that will be applied to the length of the vectors in the condition-associated vector field
#' 
#' @return list with SVD results ( d, u, v, embedding_indices, vectorField_indices )
#' @export
#'
svd.vectorField.gedi <- function(
    object,
    start.cond,
    end.cond,
    scale_cond_vector=1 ) {
  # Calculate the coordinate vectors for the start and end conditions
  startQ <- getDiffQ.gedi( object, start.cond ) + object$params$Z
  endQ <- getDiffQ.gedi(
    object, start.cond+(end.cond-start.cond)*scale_cond_vector ) +
    object$params$Z
  # Now let's map these vectors onto the Z coordinate system.
  #  In other words, we want to express these vector sets as Z %*% rot, where rot is a
  #  rotation/transformation matrix that needs to be calculated
  rotStartQ <-
    solve(eigenMatCrossprod(object$params$Z,object$params$Z)) %*%
    eigenMatCrossprod(object$params$Z,startQ)
  rotEndQ <-
    solve(eigenMatCrossprod(object$params$Z,object$params$Z)) %*%
    eigenMatCrossprod(object$params$Z,endQ)
  
  # Now, the start positions of the cells are Z %*% rotStartQ %*% DB, and end positions are Z %*% rotEndQ %*% DB
  #  We will now perform SVD on the concatenation of start and end coordinates
  projDB <- getDB.gedi( object )
  projDB <- cbind( rotStartQ %*% projDB, rotEndQ %*% projDB  )
  # Perform SVD of Z
  svd.Z <- svd(
    object$params$Z,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of projDB
  svd.projDB <- svd(
    projDB,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of the middle SVD matrices
  svd.middle <- svd(
    diag(svd.Z$d,nrow=object$aux$K) %*% t(svd.Z$v) %*%
      svd.projDB$u %*% diag(svd.projDB$d,nrow=object$aux$K),
    nu=object$aux$K, nv=object$aux$K )
  # Calculate the SVD of ZDB
  u <- svd.Z$u %*% svd.middle$u
  v <- svd.projDB$v %*% svd.middle$v
  d <- svd.middle$d
  # Set the dimension names
  names(d) <- colnames(u) <- colnames(v) <- paste0("LV",1:object$aux$K)
  rownames(u) <- object$aux$geneIDs
  rownames(v) <- c(
    paste0(object$aux$cellIDs,".start"),
    paste0(object$aux$cellIDs,".end") )
  
  return( list(
    d = d, u = u, v = v,
    embedding_indices = 1:object$aux$N,
    vectorField_indices = 1:(2*object$aux$N) ) )
}

#' Get Activity gradients
#'
#' Return the gradient of all pathways 
#' @param object GEDI object
#' 
#' @return gradient
#' @export
#' 
getActivityGradients.gedi <- function( object ) {
  A <- object$aux$C.rotation %*% object$params$A
  gradient <- t(A) # gradient is the same as the corresponding row in A
  # convert the gradient from the Z reference frame to the reference frame of genes
  gradient <- object$params$Z %*% gradient
  rownames(gradient) <- object$aux$geneIDs
  colnames(gradient) <- colnames(object$aux$inputC)
  return(gradient)
}

#' SVD activity Gradient
#'
#' Calculate the activity gradient
#' @param object GEDI object
#' @param C_index the index of the entry in the C matrix whose gradient should be calculated
#' @param scale_gradient a scaling factor that determines the magnitude of the gradient vectors.  If NA, the scaling factor will be automatically calculated so that the length of the gradient vector is 20% of the mean length of cell embedding vectors in the DB projection.
#'
#' @return list with SVD results ( d, u, v, embedding_indices, vectorField_indices )
#' @export
svd.activityGradient.gedi <- function(
    object,
    C_index,
    scale_gradient=NA ) {
  A <- object$aux$C.rotation %*% object$params$A
  gradient <- A[C_index,] # gradient is the same as the corresponding row in A
  
  # Now, the start positions of the cells are Z %*% DB, and end positions are Z %*% (DB+gradient)
  #  We will now perform SVD on the concatenation of start and end coordinates
  projDB <- getDB.gedi( object )
  if( is.na(scale_gradient) ) {
    scale_gradient <- sqrt( stats::var(c(projDB))/mean(gradient^2) ) * 0.2
    message(paste0("Gradient vectors will be scaled by a factor of ",scale_gradient,"."))
  }
  projDB <- cbind( projDB, projDB+gradient*scale_gradient  )
  # Perform SVD of Z
  svd.Z <- svd(
    object$params$Z,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of projDB
  svd.projDB <- svd(
    projDB,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of the middle SVD matrices
  svd.middle <- svd(
    diag(svd.Z$d,nrow=object$aux$K) %*% t(svd.Z$v) %*%
      svd.projDB$u %*% diag(svd.projDB$d,nrow=object$aux$K),
    nu=object$aux$K, nv=object$aux$K )
  # Calculate the SVD of ZDB
  u <- svd.Z$u %*% svd.middle$u
  v <- svd.projDB$v %*% svd.middle$v
  d <- svd.middle$d
  # Set the dimension names
  names(d) <- colnames(u) <- colnames(v) <- paste0("LV",1:object$aux$K)
  rownames(u) <- object$aux$geneIDs
  rownames(v) <- c(
    paste0(object$aux$cellIDs,".start"),
    paste0(object$aux$cellIDs,".end") )
  
  return( list(
    d = d, u = u, v = v,
    embedding_indices = 1:object$aux$N,
    gradient_indices = 1:(2*object$aux$N) ) )
}

#' SVD joint vector Field gradient
#'
#' create a umap embedding for cells and their joint vector fields both for sample condition and for pathway activity gradient
#' @param object GEDI object
#' @param start.cond contrast vector for start condition
#' @param end.cond contrast vector for end condition
#' @param C_index the index of the entry in the C matrix whose gradient should be calculated.
#' @param scale_cond_vector a scaling factor that will be applied to the length of the vectors in the condition-associated vector field
#' @param scale_gradient a scaling factor that determines the magnitude of the gradient vectors. If NA, the scaling factor will be automatically calculated so that the length of the gradient vector is 20% of the mean length of cell embedding vectors in the DB projection.
#'
#' @return list with SVD results ( d, u, v, embedding_indices, vectorField_indices, gradient_indices )
#'
#' @export
#' 
svd.joint_vectorField_gradient.gedi <- function(
    object,
    start.cond,
    end.cond,
    C_index,
    scale_cond_vector=1,
    scale_gradient=NA ) {
  A <- object$aux$C.rotation %*% object$params$A
  gradient <- A[C_index,] # gradient is the same as the corresponding row in A
  # Calculate the coordinate vectors for the start and end conditions
  startQ <- getDiffQ.gedi( object, start.cond ) + object$params$Z
  endQ <- getDiffQ.gedi(
    object, start.cond+(end.cond-start.cond)*scale_cond_vector ) +
    object$params$Z
  # Now let's map these vectors onto the Z coordinate system.
  #  In other words, we want to express these vector sets as Z %*% rot, where rot is a
  #  rotation/transformation matrix that needs to be calculated
  rotStartQ <-
    solve(eigenMatCrossprod(object$params$Z,object$params$Z)) %*%
    eigenMatCrossprod(object$params$Z,startQ)
  rotEndQ <-
    solve(eigenMatCrossprod(object$params$Z,object$params$Z)) %*%
    eigenMatCrossprod(object$params$Z,endQ)
  
  # Now, the start positions of the cells are Z %*% rotStartQ %*% DB, and end positions are Z %*% rotEndQ %*% DB
  #  We will now perform SVD on the concatenation of start and end coordinates
  projDB <- getDB.gedi( object )
  if( is.na(scale_gradient) ) {
    scale_gradient <- sqrt( stats::var(c(projDB))/mean(gradient^2) ) * 0.2
    message(paste0("Gradient vectors will be scaled by a factor of ",scale_gradient,"."))
  }
  
  projDB <- cbind(
    rotStartQ %*% projDB,
    rotEndQ %*% projDB,
    rotStartQ %*% (projDB+gradient*scale_gradient)  )
  # Perform SVD of Z
  svd.Z <- svd(
    object$params$Z,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of projDB
  svd.projDB <- svd(
    projDB,
    nu=object$aux$K, nv=object$aux$K )
  # Perform SVD of the middle SVD matrices
  svd.middle <- svd(
    diag(svd.Z$d,nrow=object$aux$K) %*% t(svd.Z$v) %*%
      svd.projDB$u %*% diag(svd.projDB$d,nrow=object$aux$K),
    nu=object$aux$K, nv=object$aux$K )
  # Calculate the SVD of ZDB
  u <- svd.Z$u %*% svd.middle$u
  v <- svd.projDB$v %*% svd.middle$v
  d <- svd.middle$d
  # Set the dimension names
  names(d) <- colnames(u) <- colnames(v) <- paste0("LV",1:object$aux$K)
  rownames(u) <- object$aux$geneIDs
  rownames(v) <- c(
    paste0(object$aux$cellIDs,".start"),
    paste0(object$aux$cellIDs,".condEnd"),
    paste0(object$aux$cellIDs,".gradEnd") )
  
  return( list(
    d = d, u = u, v = v,
    embedding_indices = 1:object$aux$N,
    vectorField_indices = 1:(2*object$aux$N),
    gradient_indices = c(1:object$aux$N,(2*object$aux$N+1):(3*object$aux$N)) ) )
  
}


#' Visualization of the vector field embedding
#'
#' @param embedding_mat Embedding from the vector field
#' @param colour vector of variable to plot
#' @param alpha alpha
#' @param randomize Logical. Whether to randomize data before plotting.
#' @param nbin number of bins
#' @param minNum minNum
#'
#' @return ggplot2 object
#' @export
#'
plot_vectorField <- function(embedding_mat,colour=1,alpha=1,randomize=T,nbin=50,minNum=10) {
  # The first half of the rows in the provided matrix are the start points,
  #  and the second half are the end points of the vector field arrows
    n <- nrow(embedding_mat)/2
  if(length(colour)==1) {
    colour=rep(colour,n)
  }
  if(length(alpha)==1) {
    alpha=rep(alpha,n)
  }
  Dim1 <- embedding_mat[1:n,1]
  Dim2 <- embedding_mat[1:n,2]
  To1 <- embedding_mat[(n+1):(n*2),1]
  To2 <- embedding_mat[(n+1):(n*2),2]
  
  # Assign the data points into elements of a grid, and calculate the average
  #  vector  per grid element
  xbinlims <- seq(
    min(Dim1), max(Dim1),
    length.out=nbin+1 )
  ybinlims <- seq(
    min(Dim2), max(Dim2),
    length.out=nbin+1 )
  xbins <- cut( Dim1, xbinlims, include.lowest = T)
  ybins <- cut( Dim2, ybinlims, include.lowest = T)
  bins <- list(paste0(xbins,ybins))
  
  embedding_obj <- data.frame(
    Dim1 = stats::aggregate(Dim1,bins,mean)[,2],
    Dim2 = stats::aggregate(Dim2,bins,mean)[,2],
    To1 = stats::aggregate(To1,bins,mean)[,2],
    To2 = stats::aggregate(To2,bins,mean)[,2],
    Alpha = stats::aggregate(alpha,bins,mean)[,2],
    n = stats::aggregate(Dim1,bins,length)[,2] )
  if( is.numeric(colour) ) {
    embedding_obj$Color <- stats::aggregate(colour,bins,mean)[,2]
  } else {
    embedding_obj$Color <- stats::aggregate(colour,bins,function(x)names(which.max(table(x))))[,2]
  }
  embedding_obj$deltaDim1 <- embedding_obj$To1 - embedding_obj$Dim1
  embedding_obj$deltaDim2 <- embedding_obj$To2 - embedding_obj$Dim2
  # Filter the grid based on the minimum required observations per grid element
  embedding_obj <- embedding_obj[ embedding_obj$n >= minNum, ]
  
  #embedding_obj <- embedding_obj[ sample(1:n,10000), ]

  # randomize the order of the objects
  if( randomize ) {
    embedding_obj <- embedding_obj[ sample.int(nrow(embedding_obj)), ]
  }
  # create the plots
  if(length(unique(colour))==1) {
    ggplot2::ggplot(
      embedding_obj, ggplot2::aes_string( x="Dim1", y="Dim2", dx="deltaDim1", dy="deltaDim2", alpha="Alpha"))+
      metR::geom_arrow(pivot=0)+
      ggplot2::theme_minimal()
  } else if(is.numeric(colour)) # the color variable is numeric
  {
    #embedding_obj$Color <- embedding_obj$Color - mean(embedding_obj$Color)
    lim <- stats::quantile(abs(embedding_obj$Color),0.99)
    ggplot2::ggplot(
      embedding_obj, ggplot2::aes_string( x="Dim1", y="Dim2", dx="deltaDim1", dy="deltaDim2", colour="Color", alpha="Alpha"))+
      metR::geom_arrow(pivot=0)+
      ggplot2::scale_color_gradientn( limits=c(-lim,lim), colours=c("blue","light grey","red"), oob=scales::squish )+
      ggplot2::theme_minimal()
  } else {
    ggplot2::ggplot(
      embedding_obj, ggplot2::aes_string( x="Dim1", y="Dim2", dx="deltaDim1", dy="deltaDim2", colour="Color", alpha="Alpha"))+
        metR::geom_arrow(pivot=0)+
      ggplot2::theme_minimal() #+guides(colour=guide_legend(override.aes=list(size=3)))
  }
}

#' Plot 2D embedding
#'
#' Plot a 2D representation (embedding) of cells
#' @param embedding_mat Embedding
#' @param colour vector of variable to plot
#' @param randomize Logical. Whether to randomize data before plotting.
#'
#' @return ggplot2 object
#' @export
#'
plot_embedding <- function(embedding_mat,colour,randomize=T) {
  # create a data frame that will have the embedding as well as the colors
  embedding_obj <- data.frame(
    Dim1=embedding_mat[,1],
    Dim2=embedding_mat[,2],
    Var=colour )

  # randomize the order of the objects
  if( randomize ) {
    embedding_obj <- embedding_obj[ sample.int(nrow(embedding_obj)), ]
  }
  # create the plots
  if(is.numeric(colour)) # the color variable is numeric
  {
    #embedding_obj$Var <- embedding_obj$Var - mean(embedding_obj$Var)
    lim <- stats::quantile(abs(embedding_obj$Var),0.99)
    ggplot2::ggplot(
      embedding_obj, ggplot2::aes_string( x="Dim1", y="Dim2", colour="Var"))+
      ggplot2::geom_point(size=0.05)+
      ggplot2::theme_minimal()+
      ggplot2::scale_color_gradientn( limits=c(-lim,lim), colours=c("blue","light grey","red"), oob=scales::squish )
  } else {
    ggplot2::ggplot(
      embedding_obj, ggplot2::aes_string( x="Dim1", y="Dim2", colour="Var"))+
      ggplot2::geom_point(size=0.05)+
      ggplot2::theme_minimal()+
      ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(size=3)))

  }
}


#######################
# Functions to calculate/plot the dispersion (variance of observation given the
#  expected value from the model)
#######################

#' Calculate dispersion
#'
#' Calculate the dispersion (variance of observation given the expected value from the model)
#'
#' @param object GEDI object
#' @param subsample subsample
#'
#' @return res dataframe
#' @export
#'
dispersion.gedi <- function( object, subsample=1e6 ) {
  res <- NULL
  for( i in 1:object$aux$numSamples ) {
    # If observation type is M:
    if( object$aux$obs.type=="M" ) {
      # Choose a subsample of the data to work with
      subsample <- min( subsample, length(object$target$Mi[[i]]) )
      samples <- sample.int(length(object$target$Mi[[i]]),subsample)
      # Get the predicted lambda based on the model
      predicted <- exp( predict_Yhat(
        object$aux$ZDBi[[i]], object$aux$QiDBi[[i]],
        object$params$si[[i]], object$params$o, object$params$oi[[i]] )[samples] )
      # Get the square of the difference between the predicted and observed values
      res.sq <- ( object$target$Mi[[i]][samples] - predicted )^2
    }
    # The rest is common across all observation types
    # Find the bins
    bins <- list( cut(
      predicted,
      stats::quantile(predicted,seq(0,1,length=sqrt(length(predicted)))),
      include.lowest=T ) )
    # For each bin, calculate the mean(predicted), which is the expected variance of
    #   Poisson, and mean((M-predicted)^2), which is the observed variance
    xy <- stats::aggregate(
      cbind( predicted, res.sq ), by=bins,
      FUN=mean )
    colnames(xy) <- c("Sample","Expected_Var","Observed_Var")
    xy$Sample <- object$aux$Samples[i]
    # Also, just calculate how many observations are in each bin
    xy$n <- stats::aggregate(
      predicted, by=bins,
      FUN=length )$x
    # Add to the results
    res <- rbind(res,xy)
  }
  return(res)
}


#' Plot dispersion
#'
#' @param dispersion Dispersion dataframe
#'
#' @return ggplot2 object
#' @export
#'
plot_dispersion <- function(dispersion) {
  ggplot2::ggplot( dispersion ) +
    ggplot2::geom_point( ggplot2::aes_string(x="Expected_Var",y="Observed_Var", color="Sample"), size=0.1 ) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype="dashed" ) +
    ggplot2::theme_minimal()
}
