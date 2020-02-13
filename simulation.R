library(robCompositions)
library(Rfast)
library(corpcor)
library(tidyverse)
library(coseq)
library(clv)

#####SFPCA algorithms
LQA <- function(X, XTX, Y, W, lambda1, lambda2, coefhot)
{
  ## Settings
  eps1 <- 1e-10
  eps2 <- 1e-3
  eps3 <- 1e-3
  num_max_iter <- 100
  tol_value <- 1e-6
  p <- ncol(W)
  
  ## Initialize coefficients
  XTY <- t(X) %*% Y
  #coef_hat <- as.vector(XTY / apply(X^2, 2, sum))
  coef_hat <- coefhot
  coef_hat[abs(coef_hat)<eps1] <- eps1
  diff_value <- 1e+10
  count <- 0
  
  while((count < num_max_iter) & (diff_value > tol_value))
  {
    print(count)
    ## Construct the Laplacian matrix for LQA
    tmpmat <- matrix(rep(coef_hat, p), p, p, byrow=F)
    tmpmat <- abs(tmpmat - W * tmpmat)
    tmpmat[tmpmat < eps1] <- eps1
    W1 <- W / tmpmat        
    D <- diag(apply(abs(W1), 1, sum))
    G <- D - W1
    V <- diag(1/abs(as.vector(coef_hat)))
    tmp1 <- XTX + lambda1 * V + lambda2 * G
    rcondnum <- 1/kappa(tmp1, exact=T)     
    if (rcondnum > 1e-15)
    {
      coef_new <- solve(tmp1, XTY)
      #coef_new <- spdinv(tmp1)%*%XTY
    }
    else
    {
      coef_new <- matrix(0, ncol(X), 1)
    }
    diff_value <- mean(abs(coef_hat - coef_new)) / mean(abs(coef_hat))
    coef_hat <- coef_new
    coef_hat[which(abs(coef_hat)<eps1),] <- eps1
    count <- count + 1
  }
  
  ## Adjust coefficients, shrink the similar loadings to the same
  idx_rest <- seq(1, p)
  size_idx <- length(idx_rest)
  while (size_idx > 1)
  {
    idx_query <- idx_rest[1]
    idx_hold <- idx_query
    for (u in idx_rest[-1])
    {
      if (abs(coef_hat[idx_query] - coef_hat[u]) < eps2)
      {
        idx_hold <- c(idx_hold, u)
      }
    }
    mm <- mean(coef_hat[idx_hold])
    if (abs(mm) < eps3)
    {
      coef_hat[idx_hold] <- 0
    }
    else
    {
      coef_hat[idx_hold] <- mm
    }
    idx_rest <- setdiff(idx_rest, idx_hold)
    size_idx <- length(idx_rest)
  }   
  
  output <- list()
  output$coef <- coef_hat
  return(output)
}
## Initialize A by ordinary PCA
initializeA <- function(X, K)
{
  covmat <- cov(X)
  A <- fast.svd(covmat)$v[, seq(1, K)]
  A_init <- as.matrix(A)
  return(A_init)
}    
## This function implement PCA with generalized elasticnet algorithm with estiamted covariance matrix of X
SFPCA <- function(X, lambda1, lambda2, K, centering=TRUE)
{
  ## Initialize the parameters controlling the loop
  tol1 <- 1e-6
  max_iter1 <- 100
  count <- 0
  diff_value <- 1e+10
  degree_of_freedom <- rep(0, K)
  
  
  ## If centering==TRUE, center the data
  X <- scale(X, center=TRUE, scale=FALSE)
  p <- ncol(X)  ## Number of variables
  n <- nrow(X)  ## Number of observations
  XTX <- t(X) %*% X
  
  ## Preprocess lambda1, lambda2. If lambda1 and lambda2 are scalars, extend them to vectors.
  if (length(lambda1)==1)
  {
    lambda1 <- rep(lambda1, K)
  }
  if (length(lambda2)==1)
  {
    lambda2 <- rep(lambda2, K)
  }
  
  ## Creat similarity matrix W
  covmat <- XTX / n
  ww <- diag(1/sqrt(diag(covmat)))
  similarity_mat <- ww %*% covmat %*% ww
  W <- similarity_mat - diag(diag(similarity_mat))
  
  ## Initialize A given B
  A <- initializeA(X, K)
  
  ## Initialize B given A
  B <- A
  B_new <- B
  for (i in seq(1, K))
  {
    alpha_value <- A[, i]
    trY <- X %*% alpha_value
    B[, i] <- LQA(X, XTX, trY, W, lambda1[i], lambda2[i], coefhot=B[, i])$coef
  }
  
  ## Loop, update A and B iteratively
  while ((diff_value >tol1) & (count < max_iter1))
  {
    ## Given B, optimize A by Procrustes rotation. Complexity O(p*K^2)
    temp1 <- XTX %*% B
    z <- fast.svd(temp1)
    A <- (z$u) %*% t(z$v)
    
    ## Given A, optimize B, solve a number of generalize elasticnet problems
    for (i in seq(1, K))
    {
      message('lambda1: ', lambda1[i],'; lambda2: ',lambda2[i])
      
      alpha_value <- A[, i]
      trY <- X %*% alpha_value
      obj_qa <- LQA(X, XTX, trY, W, lambda1[i], lambda2[i], coefhot=B[, i])
      B_new[, i] <- obj_qa$coef
    }
    
    ## Normalize beta_i
    normbeta <- sqrt(apply(B_new^2, 2, sum))
    normbeta[normbeta == 0] <- 1
    B_new <- t(t(B_new)/normbeta)
    
    ## Control the loop
    #diff_value <- sum(abs(B_new - B)) / sum(abs(B))
    diff_value <- max(abs(B_new - B))
    B <- B_new
    count <- count + 1
  }
  
  output <- list()
  output$A <- A
  output$B <- B
  Z <- X %*% B
  tmp1 <- qr(Z)
  R <- qr.R(tmp1)
  tmp2 <- diag(R)^2 / nrow(X)
  output$adjvar <- tmp2[seq(1, K)]
  return(output)
}


#GMPR
GMPR <- function (comm, intersect.no = 10, ct.min = 0, trace = TRUE) {
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {		
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))		
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n', 
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}


N<-6

#Generate Initial center composition 
snapper <- c(0.4,0.4,0.04,0.04,0.04,0.04,0.04)
chatter <- c(0.04,0.04,0.4,0.4,0.04,0.04,0.04)
snap_as_camera <- c(0.01,0.01,0.01,0.01,0.01,0.01,0.94)
snap_chat_viewer <- c(0.04,0.4,0.04,0.4,0.04,0.04,0.04)
story_viewer <- c(0.01,0.01,0.01,0.01,0.01,0.9,0.05)
story_poster <- c(0.01,0.01,0.01,0.01,0.6,0.01,0.35)

#Persona Composition
p <- c(0.25,0.25,0.05,0.1,0.25,0.1)
total<- 500000

total*p

generate_comp <- function(num,comp){
  tmp <- mvrnorm(num,comp,Sigma=diag(comp)/50) %>% as_tibble
  tmp[tmp<0]<-0
  tmp %<>% mutate(row_sum=rowSums(.)) %>% mutate_all(~./row_sum) %>% dplyr::select(-row_sum)
  #tmp %<>% normalize(.,byrow=TRUE) 
  return(tmp)
}

snapper_data <- generate_comp(total*p[1],snapper) %>% mutate(label='snapper')
chatter_data <- generate_comp(total*p[2],chatter) %>% mutate(label='chatter')
snap_as_camera_data <- generate_comp(total*p[3],snap_as_camera) %>% mutate(label='snap_as_camera')
snap_chat_viewer_data <- generate_comp(total*p[4],snap_chat_viewer) %>% mutate(label='snap_chat_viewer')
story_viewer_data <- generate_comp(total*p[5],story_viewer) %>% mutate(label='story_viewer')
story_poster_data <- generate_comp(total*p[6],story_poster) %>% mutate(label='story_poster')

total_comp <- rbind(snapper_data,chatter_data,snap_as_camera_data,snap_chat_viewer_data,story_viewer_data,story_poster_data)
total_eng <- rnbinom(total,mu=1000,size=0.1)+6

total_comp %<>% mutate(total_eng=total_eng)

total_count<- total_comp %>% mutate_at(c('V1','V2','V3','V4','V5','V6','V7'),~.*total_eng)
total_count_rounded<-total_count %>% mutate_at(c('V1','V2','V3','V4','V5','V6','V7'),round)
colnames(total_count_rounded)[1:7]<-c('Snap_Send','Snap_View','Chat_Send','Chat_View','Story_Post','Story_View','Filter_Swipe')

norm_factor<-GMPR(total_count_rounded[,1:7] %>% as.matrix)

#Transform the counts by GMPR scaling factor
total_count_transformed<-t(t(total_count_rounded[,1:7])/c(norm_factor))

#Avoid 
total_count_transformed<- total_count_transformed+0.1

#Transform back into compositions
comp_transformed<- total_count_transformed %>% as_tibble %>% mutate(row_sum=rowSums(.)) %>%  mutate_at(c(1,2,3,4,5,6,7),~./row_sum)

#Compositional Correlation Analysis
#comp_correlation <- corCoDa(comp_transformed[,1:7] %>% as.matrix)

#K-means clustering using slope heuristic

res<-coseq(comp_transformed[,1:7] %>% as.matrix,model='kmeans',2:12,transformation='clr',normFactors='none',parallel=TRUE)
conds<-colnames(comp_transformed[,1:7])
plot(res, graphs="boxplots",conds=conds, collapse_reps = "average")

#Clustering assignment
cluster_res<-clusters(res)

cluster_res[which(cluster_res==4)]<-'snap_as_camera'
cluster_res[which(cluster_res==2)]<-'story_viewer'
cluster_res[which(cluster_res==1)]<-'snap_chat_viewer'
cluster_res[which(cluster_res==5)]<-'snapper'
cluster_res[which(cluster_res==3)]<-'story_poster'
cluster_res[which(cluster_res==6)]<-'chatter'

total_comp %<>% mutate(res=cluster_res)

