library(rstan)
library(ape)

#uncomment if files are not already in working directory
#download.file('https://osf.io/czaym/metadata/?format=datacite-json', 'Indo-European.cc_sc.phy', mode = "wb") 
#download.file('https://osf.io/h86um/metadata/?format=datacite-json', 'geoDists.csv', mode = "wb") 

df = readLines('Indo-European.cc_sc.phy')
nlangs = strsplit(df[1], split = ' ')[[1]][1]
nsites = strsplit(df[1], split = ' ')[[1]][2]
df = df[-1]


chars = c()
langs = c()

for(i in 1:nlangs){
  temp = gsub(' +', ' ', df[i])
  temp = strsplit(temp, split = ' ')[[1]]
  langs = c(langs, rep(temp[1], nsites))
  chars = c(chars, as.numeric(strsplit(temp[2], split = '')[[1]]))
  
}

incllist = c(
  "IE.ARMENIAN.EASTERN_ARMENIAN", "IE.ARMENIAN.WESTERN_ARMENIAN", "IE.INDIC.URDU", "IE.INDIC.HINDI", "IE.IRANIAN.PERSIAN",
  "IE.GREEK.GREEK","IE.CELTIC.IRISH_GAELIC","IE.BALTIC.LATVIAN","IE.ALBANIAN.ALBANIAN", "IE.ROMANCE.ITALIAN", "IE.ROMANCE.FRENCH",
  "IE.ROMANCE.CORSICAN", "IE.GERMANIC.DANISH", "IE.GERMANIC.DUTCH","IE.GERMANIC.ENGLISH", "IE.GERMANIC.ICELANDIC",
  "IE.SLAVIC.UKRAINIAN", "IE.SLAVIC.RUSSIAN", "IE.SLAVIC.CZECH", "IE.SLAVIC.POLISH", "IE.CELTIC.GAELIC_SCOTTISH", "IE.ROMANCE.ROMANIAN",
  "IE.SLAVIC.BULGARIAN", "IE.SLAVIC.MACEDONIAN"
)

langs2 = langs[langs %in% incllist]
chars2 = chars[langs %in% incllist]

langs2 = gsub('.*\\.', '', langs2)
incllist = gsub('.*\\.', '', incllist)


geoDists <- read.table('geoDists.csv',
                       sep=',',row.names=1,header=F)
rownames(geoDists) <- sub('-','_',rownames(geoDists))
colnames(geoDists) <- rownames(geoDists)

distm = matrix(NA, nrow = length(incllist), ncol = length(incllist))
rownames(distm) <- incllist
colnames(distm) <- incllist

for (i in 1:length(incllist)) {
  for (e in 1:length(incllist)) {
    if(e>i){
      distm[i,e] = geoDists[incllist[i], incllist[e]]
        
      distm[e,i] = distm[i,e]
    }
  }
  distm[i,i] = 0
}

chars3 = c()
langs3 = c()
langlevelled = as.numeric(factor(langs2, levels = incllist))
tmpdf = data.frame('ID' = 1:as.numeric(nsites))
for (i in 1:length(unique(langlevelled))){
  tmpdf = cbind(tmpdf, data.frame(i = chars2[which(langlevelled == i)]))
}
tmpdf = tmpdf[,-1]

tmpdf = tmpdf[rowSums(tmpdf, na.rm = T)!=0,]
tmpdf = tmpdf[rowSums(is.na(tmpdf))==0,]

nsites = nrow(tmpdf)
for(i in 1:nrow(tmpdf)){
  chars3 = c(chars3, as.numeric(tmpdf[i,]))
  langs3 = c(langs3, c(1:24))
}

dat_phylo = list(
  C = chars3,
  L = langs3,
  Nsites = as.numeric(nsites),
  Nlangs = length(incllist),
  N = length(chars3),
  Ncomb = choose(length(incllist),2),
  Dmat = distm/1000
)

dat_phylo$C[is.na(chars3)] = -1



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

save(mod, file = 'stanmodel_ASJP.rda')