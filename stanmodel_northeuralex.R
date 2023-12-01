library(rstan)
library(ape)

#uncomment if files are not already in working directory
#download.file('http://www.sfs.uni-tuebingen.de/~jdellert/northeuralex/0.9/northeuralex-0.9-language-data.tsv', 'northeuralex-0.9-language-data.tsv', mode = "wb") 
#download.file('http://www.sfs.uni-tuebingen.de/~jdellert/pubs/jdellert-diss-supplements.tar.gz', 'jdellert-diss-supplements.tar.gz', mode = "wb") 
#untar("jdellert-diss-supplements.tar.gz",files="inferred-cognates.tsv")

df = read.csv('inferred-cognates.tsv', sep = '\t', header = F)

newdf = data.frame()

for(i in 1:length(unique(df$V1))){
  print(i)
  subdf = df[df$V1 == unique(df$V1)[i],]
  tmpdf = data.frame(matrix(ncol = length(unique(df$V2)), nrow = length(unique(subdf$V5))))
  colnames(tmpdf) <- unique(df$V2)
  for (row in 1:nrow(tmpdf)) {
    yesin = unique(subdf$V2[subdf$V5 == row])
    noin = unique(subdf$V2[subdf$V5 != row])
    tmpdf[row,yesin] = 1
    tmpdf[row,setdiff(noin,yesin)] = 0
  }
  newdf = rbind(newdf, tmpdf)
  
}

newdf = newdf[,c('isl', 'eng', 'nor', 'nld',
                 'dan', 'deu', 'smj', 'sma',
                 'swe', 'sme', 'fin', 'ekk',
                 'liv', 'lav', 'lit', 'pol',
                 'bel', 'smn', 'sms', 'krl',
                 'sjd', 'olo', 'vep', 'rus'
  
  
)]
newdf = newdf[rowSums(newdf) != 0,]
newdf = newdf[rowSums(newdf) != 1,]
newdf = newdf[rowSums(newdf) != ncol(newdf),]
newdf = newdf[rowSums(is.na(newdf))==0,]

chars = c()
langs = c()
for(i in 1:nrow(newdf)){
  chars = c(chars,as.vector(newdf[i,], mode = "numeric"))
  langs = c(langs,colnames(newdf))
}

geodat = read.csv("northeuralex-0.9-language-data.tsv", sep = '\t', header = T)
geodat = geodat[,c(3,6,7)]
geodat = geodat[geodat$iso_code %in% colnames(newdf), ]
geodat$iso_code = as.character(geodat$iso_code)
geodat = geodat[match(colnames(newdf),geodat$iso_code),]

library(geosphere)
geoDists = distm(geodat[,c(3,2)], fun = distHaversine)/1000000

dat_phylo = list(
  C = chars,
  L = as.numeric(as.factor(langs)),
  Nsites = nrow(newdf),
  Nlangs = ncol(newdf),
  N = length(chars),
  Ncomb = choose(ncol(newdf),2),
  Dmat = geoDists
)

dat_phylo$C[is.na(chars)] = -1


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



mod <- stan(model_code = stan_program, data = dat_phylo, iter=4000, chains = 4, cores = 4, init = 0)

save(mod, file = '/data/archiv/SCC/fhartmann/readymodel_northeuralex_confounded.rda')