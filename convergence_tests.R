library(rstan)
library(phangorn)
##########

load('matrixsamples_stanmodel_northeuralex.rda')
posterior_draws = as.data.frame(smpls)
posterior_trees = c()

langvec = c('isl', 'eng', 'nor', 'nld',
            'dan', 'deu', 'smj', 'sma',
            'swe', 'sme', 'fin', 'ekk',
            'liv', 'lav', 'lit', 'pol',
            'bel', 'smn', 'sms', 'krl',
            'sjd', 'olo', 'vep', 'rus')
geodat = read.csv("northeuralex-0.9-language-data.tsv", sep = '\t', header = T)
geodat = geodat[geodat$iso_code %in% langvec, ]
geodat$iso_code = as.character(geodat$iso_code)

geodat = geodat[match(langvec, geodat$iso_code),]
incllist = geodat$name
nlangs = length(incllist)


library(stringr)

for(r in 1:8000){
  print(r)
  mymatrix = matrix(0, nrow = nlangs, ncol = nlangs,dimnames = list(incllist,incllist))
  counter = 1
  for(i in 1:nlangs){
    for(e in 1:nlangs){
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
  print(i)
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

sfreq = aggregate(pp ~ split, resdf2, sd)
mean(sfreq$pp)# 0.069, #0.031

ch1 = posterior_trees2[1:2000]
ch2 = posterior_trees2[2001:4000]
ch3 = posterior_trees2[4001:6000]
ch4 = posterior_trees2[6001:8000]

plot(consensus(ch1, rooted = T))
plot(consensus(ch2, rooted = T))
plot(consensus(ch3, rooted = T))
plot(consensus(ch4, rooted = T))


######
load('matrixsamples_stanmodel_ASJP.rda')
posterior_draws = as.data.frame(smpls)


posterior_trees = c()

incllist = c(
  
  "IE.ARMENIAN.EASTERN_ARMENIAN", "IE.ARMENIAN.WESTERN_ARMENIAN", "IE.INDIC.URDU", "IE.INDIC.HINDI", "IE.IRANIAN.PERSIAN",
  "IE.GREEK.GREEK","IE.CELTIC.IRISH_GAELIC","IE.BALTIC.LATVIAN","IE.ALBANIAN.ALBANIAN", "IE.ROMANCE.ITALIAN", "IE.ROMANCE.FRENCH",
  "IE.ROMANCE.CORSICAN", "IE.GERMANIC.DANISH", "IE.GERMANIC.DUTCH","IE.GERMANIC.ENGLISH", "IE.GERMANIC.ICELANDIC",
  "IE.SLAVIC.UKRAINIAN", "IE.SLAVIC.RUSSIAN", "IE.SLAVIC.CZECH", "IE.SLAVIC.POLISH", "IE.CELTIC.GAELIC_SCOTTISH", "IE.ROMANCE.ROMANIAN",
  "IE.SLAVIC.BULGARIAN", "IE.SLAVIC.MACEDONIAN"
)
incllist = gsub('.*\\.', '', incllist)
nlangs = length(incllist)
incllist = gsub('_', ' ', incllist)
library(stringr)
incllist = str_to_title(incllist)

for(r in 1:8000){
  print(r)
  mymatrix = matrix(0, nrow = nlangs, ncol = nlangs,dimnames = list(incllist,incllist))
  counter = 1
  for(i in 1:nlangs){
    for(e in 1:nlangs){
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
  print(i)
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

sfreq = aggregate(pp ~ split, resdf2, sd)
mean(sfreq$pp)# 0.016 #0.077


ch1 = posterior_trees2[1:2000]
ch2 = posterior_trees2[2001:4000]
ch3 = posterior_trees2[4001:6000]
ch4 = posterior_trees2[6001:8000]

plot(consensus(ch1, rooted = T, p=0.5))
plot(consensus(ch2, rooted = T, p=0.5))
plot(consensus(ch3, rooted = T, p=0.5))
plot(consensus(ch4, rooted = T, p=0.5))
