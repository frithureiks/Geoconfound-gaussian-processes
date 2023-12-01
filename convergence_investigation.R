library(rstan)
library(ape)
library(MASS)
library(boot)
library(stringr)
library(phangorn)

set.seed(42)


stan_program <- "
functions{
  
  
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-(square(x[i,j])/sq_rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
  
  matrix cov_OU(matrix x, real sq_alpha, real rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    
    
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-(abs(x[i,j])/rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
  
  matrix matrx_conv(vector Dmat2_vec, int Nlangs) {
  matrix[Nlangs,Nlangs] Dmat2;
  int counter;
  counter = 1;
  for (i in 1:(Nlangs-1)) {
    Dmat2[i, i] = 0;
    for (j in (i + 1):Nlangs) {
      Dmat2[i, j] = Dmat2_vec[counter];
      Dmat2[j, i] = Dmat2[i, j];
      counter = counter+1;
    }
  }
  Dmat2[Nlangs, Nlangs] = 0;
  return Dmat2;
  }
  
  
}

data{
  int N;
  int Nlangs;
  int Nsites;
  int Ncomb;
  int C[N];
  int L[N];
  matrix[Nlangs,Nlangs] Dmat;
}
parameters{
  
  real<lower=0> etasq;
  real<lower=0> rhosq;
  
  real<lower=0> etasq2;

  vector<lower=0>[Ncomb] Dmat2_vec;
  
  vector[N] eta_k;
  
  vector[N] eta_r;
}
transformed parameters{
  

  vector[N] r;
  vector[N] k;
  vector[N] p; 
  
{
  

  matrix[Nlangs,Nlangs] SIGMA;
  matrix[Nlangs,Nlangs] SIGMA2;
  matrix[Nlangs, Nlangs] L_K;
  matrix[Nlangs, Nlangs] L_K2;
 

  
  SIGMA = cov_GPL2(Dmat, etasq, rhosq, 0.01);
  SIGMA2 = cov_OU(matrx_conv(Dmat2_vec, Nlangs), etasq2, 3, 0.01);

  
  L_K = cholesky_decompose(SIGMA);
  L_K2 = cholesky_decompose(SIGMA2);
  
  for (e in 1:Nsites){

    k[((e*Nlangs)-Nlangs+1):(e*Nlangs)] = L_K * eta_k[((e*Nlangs)-Nlangs+1):(e*Nlangs)];
    r[((e*Nlangs)-Nlangs+1):(e*Nlangs)] = L_K2 * eta_r[((e*Nlangs)-Nlangs+1):(e*Nlangs)];
  
  }
  
}
  p = inv_logit(r + k);
  
}

model{

  Dmat2_vec ~ exponential( 1 );
  etasq2 ~ exponential( 2 );
  rhosq ~ exponential( 2 );
  etasq ~ exponential( 2 );
  eta_r ~ normal(0,1);
  eta_k ~ normal(0,1);  
  
  C ~ binomial( 1 , p);
  
}

"

charGen <- function(phytree, nsites){
  C = c()
  counter = 1
  while (counter <= nsites) {
    innovedge = sample(1:nrow(phytree$edge), size=1, prob=phytree$edge.length/sum(phytree$edge.length))
    
    ones = rep(0, max(phytree$edge[,2]))
    if (phytree$edge[innovedge,2] > length(phytree$tip.label)){
      
      descs = allDescendants(phytree)[phytree$edge[innovedge,2]][[1]]
      ones[descs] =  1
      descedges = which(phytree$edge[,2] %in% descs)
      loseprobs = phytree$edge.length[descedges]*0.1
      for (i in 1:length(descs)){
        if( runif(1) < loseprobs[i]){
          losedescs = allDescendants(phytree)[descs[i]][[1]]
          ones[losedescs] =  0
          
        }
      }
    }else{
      ones[phytree$edge[innovedge,2]] =  1
    }
    
    ones = ones[1:length(phytree$tip.label)]
    if(length(unique(ones)) > 1){
      C = c(C, ones)
      counter = counter + 1
    }
  }
  return(C)
}


nlangs = c(5,10,20)
nsites = c(100, 300, 500)
covariant = c('yes', 'no')

for (run in 1:4){
  for (nlang in nlangs){
    for (nsite in nsites){
      for (covariance in covariant){
        
          phytree = rcoal(nlang)#random tree
          phydist = cophenetic(phytree)
          
          if (covariance == 'yes'){
            geodist = phydist*1.2
          } else {
            geoloc = cbind(runif(nlang, 0, 3),runif(nlang, 0, 3))
            geodist = as.matrix(dist(geoloc, diag = T, upper = T))
          }
        
          
          C = charGen(phytree, nsite)
          L = rep(c(1:nlang), nsite)

          
          dat_phylo = list(
            C = C,
            L = L,
            Nsites = nsite,
            Nlangs = nlang,
            N = length(C),
            Ncomb = nlang*nlang,
            Dmat = geodist
          )
          mod <- stan(model_code = stan_program, data = dat_phylo, iter=4000, chains = 4, cores = 4, init = 0)
          extr = extract(mod, permuted = F)
          
          rhats = c()
          for (i in 1:dim(extr)[3]){
            rhats = c(rhats, Rhat(extr[,,i]))
          }
          
          extr = extract(mod, permuted = F, pars = c('Dmat2_vec'))
          rhats_PhyMat = c()
          for (i in 1:dim(extr)[3]){
            rhats_PhyMat = c(rhats_PhyMat, Rhat(extr[,,i]))
          }
          
          smpls = extract(mod, pars = c('Dmat2_vec'))
          
          posterior_draws = as.data.frame(smpls)
          posterior_trees = c()
    
          for(r in 1:8000){
            mymatrix = matrix(0, nrow = nlang, ncol = nlang)
            counter = 1
            for(i in 1:nlang){
              for(e in 1:nlang){
                if(e>i){
                  mymatrix[i,e] = posterior_draws[r,counter]
                  mymatrix[e,i] = mymatrix[i,e]
                  counter = counter+1
                }
              }
            }
            cluster = as.phylo(hclust(as.dist(mymatrix), method = 'average'))
            if(r == 1){
              prevclust = cluster
            }else if (r == 2){
              posterior_trees = c(prevclust,cluster)
            }else{
              
              posterior_trees = c(posterior_trees, cluster)
            }
            
            
          }
          posterior_trees2 = c()
          
          taxa = posterior_trees[[1]]$tip.label
          
          for (i in 1:length(posterior_trees)){
            cluster = rotateConstr(posterior_trees[[i]], taxa)
            if(i == 1){
              prevclust = cluster
            }else if (i == 2){
              posterior_trees2 = c(prevclust,cluster)
            }else{
              
              posterior_trees2 = c(posterior_trees2, cluster)
            }
          }
          
          allpp = prop.part(posterior_trees2)
          pp1 = prop.part(posterior_trees2[1:2000])
          pp2 = prop.part(posterior_trees2[2001:4000])
          pp3 = prop.part(posterior_trees2[4001:6000])
          pp4 = prop.part(posterior_trees2[6001:8000])
          chainlist = list(pp1,pp2,pp3,pp4)
          
          resdf = data.frame()
          for(i in 1:length(allpp)){
            part = allpp[[i]]
            for (chain in 1:4) {
              mychain = chainlist[[chain]]
              for(e in 1:length(mychain)){
                if(length(mychain[[e]]) == length(part)){
                  if(all(mychain[[e]] == part)){
                      resdf = rbind(resdf, data.frame('split' = i, 'chain' = chain, 'pp' = attr(mychain,"number")[e]/2000))
                  }
                }
              }
            }
          }
          
          resdf2 = resdf[resdf$split %in% as.numeric(names(table(resdf$split))[which(as.vector(table(resdf$split)) == 4)])
                         ,]
          sfreq = aggregate(pp ~ split, resdf2, sd)
          
          #only majority clades
          resdf = data.frame()
          for(i in 1:length(allpp)){
            part = allpp[[i]]
            for (chain in 1:4) {
              mychain = chainlist[[chain]]
              for(e in 1:length(mychain)){
                if(length(mychain[[e]]) == length(part)){
                  if(all(mychain[[e]] == part)){
                    if(attr(mychain,"number")[e]/2000 > 0.5){
                      resdf = rbind(resdf, data.frame('split' = i, 'chain' = chain, 'pp' = attr(mychain,"number")[e]/2000))
                    }
                  }
                }
              }
            }
          }
          
          resdf2 = resdf[resdf$split %in% as.numeric(names(table(resdf$split))[which(as.vector(table(resdf$split)) == 4)])
                         ,]
          sfreq_only_majority = aggregate(pp ~ split, resdf2, sd)
          
          tmpres = data.frame(
            'nlangs' = nlang,
            'nsites' = nsite,
            'covariance' = covariance,
            'rhat_mean' = mean(rhats),
            'rhat_phymat_mean' = mean(rhats_PhyMat),
            'splitfreq_mean' = mean(sfreq$pp),
            'splitfreq_majority_mean' = mean(sfreq_only_majority$pp),
            'majority_clades' = nrow(sfreq_only_majority)
          )
          write.table(tmpres, file = 'gridsearch_res.csv', sep = ",", 
                      append = TRUE, quote = FALSE, 
                      col.names = FALSE, row.names = FALSE) 
      }
    }
  }
}


