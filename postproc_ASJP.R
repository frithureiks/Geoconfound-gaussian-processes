library(rstan)
library(rethinking)

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
  Dmat = distm/1000,
  miss = as.numeric(is.na(chars3))
)

dat_phylo$C[is.na(chars3)] = -1


load(file ='fullsummary_gerharddat.rda')
shortsumm = smry$summary
res = data.frame('true' = dat_phylo$C, 'pred' = NA)

res$r = shortsumm[grepl('^r\\[', rownames(shortsumm)),1]
res$k = shortsumm[grepl('^k\\[', rownames(shortsumm)),1]
p =  shortsumm[grepl('^r\\[', rownames(shortsumm)),1] + 
  shortsumm[grepl('^k\\[', rownames(shortsumm)),1]

p = inv_logit(p)
res$pred = p

res$d = res$true-res$pred
hist(res$d)
sum(abs(res$d) < 0.3)/nrow(res)

res$disc = 0
res$disc[res$pred >0.5] = 1

caret::confusionMatrix(data = as.factor(res$disc), reference = as.factor(res$true))


#####
langvec = c(
  "IE.ARMENIAN.EASTERN_ARMENIAN", "IE.ARMENIAN.WESTERN_ARMENIAN", "IE.INDIC.URDU", "IE.INDIC.HINDI", "IE.IRANIAN.PERSIAN",
  "IE.GREEK.GREEK","IE.CELTIC.IRISH_GAELIC","IE.BALTIC.LATVIAN","IE.ALBANIAN.ALBANIAN", "IE.ROMANCE.ITALIAN", "IE.ROMANCE.FRENCH",
  "IE.ROMANCE.CORSICAN", "IE.GERMANIC.DANISH", "IE.GERMANIC.DUTCH","IE.GERMANIC.ENGLISH", "IE.GERMANIC.ICELANDIC",
  "IE.SLAVIC.UKRAINIAN", "IE.SLAVIC.RUSSIAN", "IE.SLAVIC.CZECH", "IE.SLAVIC.POLISH", "IE.CELTIC.GAELIC_SCOTTISH", "IE.ROMANCE.ROMANIAN",
  "IE.SLAVIC.BULGARIAN", "IE.SLAVIC.MACEDONIAN"
)

langvec = gsub('.*\\.', '', incllist)
datset = read.csv(file = 'dataset.tab', sep = '\t', fill = T, header = TRUE)

latlon_df = data.frame()
for (i in 1:length(incllist)){
  
  latlon_df = rbind(latlon_df, data.frame('lang' = incllist[i], 'lat' = datset[which(datset$names == incllist[i]), 'lat'],
                                          'lon' = datset[which(datset$names == incllist[i]), 'lon']))
}

matrconv <- function(Dmat2_vec, Nlangs) {
  Dmat2 = matrix(NA, nrow = Nlangs, ncol = Nlangs)
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
  return (Dmat2);
}


c1_vals = smry$summary[grepl('Dmat2_vec', rownames(smry$summary)),1]
c1_m = matrconv(c1_vals, length(langvec))
dimnames(c1_m) <-list(langvec,langvec)
my_upgma <- phangorn::upgma(as.dist(c1_m))


c1_vals = smry$c_summary[grepl('Dmat2_vec', rownames(smry$summary)),,1]
c1_m = matrconv(c1_vals, length(langvec))
dimnames(c1_m) <-list(langvec,langvec)
my_upgma <- phangorn::upgma(as.dist(c1_m))


cairo_pdf("tree_asjp_chain1.pdf", width = 12, height = 8)

plot(my_upgma)

dev.off()

c1_vals = smry$c_summary[grepl('Dmat2_vec', rownames(smry$summary)),,2]
c1_m = matrconv(c1_vals, length(langvec))
dimnames(c1_m) <-list(langvec,langvec)
my_upgma <- phangorn::upgma(as.dist(c1_m))


cairo_pdf("tree_asjp_chain2.pdf", width = 12, height = 8)

plot(my_upgma)

dev.off()



####


Nlang = length(langvec)

K <- matrix(0,nrow=Nlang,ncol=Nlang)
for ( i in 1:Nlang ){
  print(i/Nlang)
  for ( j in 1:Nlang){
    K[i,j] <- smry$summary[grepl('etasq$', rownames(smry$summary)),1] * exp( -(dat_phylo$Dmat[i,j]^2 / smry$summary[grepl('rhosq', rownames(smry$summary)),1]) )
  }}
diag(K) <- smry$summary[grepl('etasq$', rownames(smry$summary)),1] + 0.01

Rho <- round(cov2cor(K), 2 )
colnames(Rho) <- langvec
rownames(Rho) <- colnames(Rho)

langvec = gsub('_', ' ', langvec)
library(stringr)
langvec = str_to_title(langvec)
library(ggplot2)
cairo_pdf("geoconfound_strength_asjp.pdf", width = 12, height = 8)


plot( NULL, xlab="Longitude" , ylab="Latitude" ,
      col=col.alpha("black",0.8) , cex=1 , pch=1, xlim=c(-30,90), ylim=c(20,70), lwd = 3  )
labels <- langvec

world <- map_data("world", xlim=c(-30,90), ylim=c(20,70))

for(g in unique(world$group)){
  world2 = world[world$group == g,]
  polygon(world2[,1:2],
          col = NULL, border = col.alpha("black",0.3))
}

for( i in 1:Nlang ){
  for ( j in 1:Nlang ){
    if ( i < j ){
      lines( c( latlon_df$lon[i],latlon_df$lon[j] ) , c( latlon_df$lat[i],latlon_df$lat[j] ) ,
             lwd=2 , col=col.alpha("blue",Rho[i,j]^5)) 
    }}}
points( latlon_df$lon , latlon_df$lat,
      col=col.alpha("black",1) , cex=1 , pch=1, lwd = 3  )



dev.off()



cairo_pdf("geoconfound_strength_asjp2.pdf", width = 12, height = 8)


plot( NULL, xlab="Longitude" , ylab="Latitude" ,
      col=col.alpha("black",0.8) , cex=1 , pch=1, xlim=c(-30,90), ylim=c(20,70), lwd = 3  )
labels <- langvec

for( i in 1:Nlang ){
  for ( j in 1:Nlang ){
    if ( i < j ){
      lines( c( latlon_df$lon[i],latlon_df$lon[j] ) , c( latlon_df$lat[i],latlon_df$lat[j] ) ,
             #lwd=2 , col=col.alpha("black",Rho[i,j]^8 /max(Rho[upper.tri(Rho^8)]^8))) 
             lwd=2 , col=col.alpha("blue",Rho[i,j]^5)) 
    }}}
points( latlon_df$lon , latlon_df$lat,
        col=col.alpha("black",1) , cex=1 , pch=1, lwd = 3  )
text( latlon_df$lon , latlon_df$lat, labels=labels , cex=0.8 , pos =c(3,1,2,4,4,1,2,3,2,2,2,2,3,4,2,2,4,4,4,4,4,4,4,4))

dev.off()


####


library(ape)
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


library(phangorn)
cns = consensus(posterior_trees2, p=.5, rooted = T)

boot <- prop.clades(cns, posterior_trees2, rooted = T)
cns$node.label = round(boot/8000, digits = 2)

cairo_pdf("tree_asjp.pdf", width = 8, height = 6)
ggtree(cns) + geom_tiplab() + xlim(0, 8) + geom_nodelab(hjust=-.1, size=3)
dev.off()


mcc = maxCladeCred(posterior_trees2)
load('matrixsamples_stanmodel_ASJP_confounded.rda')
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
posterior_trees2_cnf = c()

taxa = posterior_trees[[1]]$tip.label

for (i in 1:length(posterior_trees)){
  print(i)
  cluster = rotateConstr(posterior_trees[[i]], taxa)
  if(i == 1){
    prevclust = cluster
  }else if (i == 2){
    posterior_trees2_cnf = c(prevclust,cluster)
  }else{
    
    posterior_trees2_cnf = c(posterior_trees2_cnf, cluster)
  }
}



library(phangorn)

cns2 = consensus(posterior_trees2_cnf, p=.5, rooted = T)
boot <- prop.clades(cns, posterior_trees2_cnf, rooted = T)
cns$node.label = round(boot/8000, digits = 2)


cairo_pdf("tree_asjp_confounded.pdf", width = 8, height = 6)
ggtree(cns2) + geom_tiplab() + xlim(0, 8) + geom_nodelab(hjust=-.1, size=3)
dev.off()

obj<-cophylo(cns2,cns, rotate = T)


cairo_pdf("treecompare_asjp.pdf", width = 8, height = 6)
plot(obj,fsize=c(0.8,0.8),link.type="curved",link.lwd=2, lty = 'dashed', ylim = c(0,1), align = T)
dev.off()

pp = prop.part(posterior_trees2)
pp_c = prop.part(posterior_trees2_cnf)

pp_comp = data.frame()
for(e in 1:length(pp_c)){
  tmpdf = data.frame()
  for(i in 1:length(pp)){
    if(all(pp_c[[e]] == pp[[i]])){
      tmpdf = data.frame('i_c' = e, 'i' = i)
    }
  }
  if(nrow(tmpdf) == 0){
    tmpdf = data.frame('i_c' = e, 'i' = NA)
  }
  pp_comp = rbind(pp_comp, tmpdf)
}


pp_comp$supp_c = attr(pp_c, 'number')/8000
pp_comp$supp = c(attr(pp, 'number')/8000)[pp_comp$i]
pp_comp$supp[is.na(pp_comp$supp)] = 0

pp_comp$d =  pp_comp$supp-pp_comp$supp_c


pp_comp$labels = NA

for(e in 1:length(pp_c)){
  pp_comp$labels[e] = paste(attr(pp_c, 'labels')[pp_c[[e]]], collapse = ', ')
}


head(pp_comp[order(pp_comp$d),], n= 5)

xtable::xtable(head(pp_comp[order(pp_comp$d),], n= 5))


head(pp_comp[order(pp_comp$d, decreasing = T),], n= 5)
xtable::xtable(head(pp_comp[order(pp_comp$d, decreasing = T),], n= 5))

####matrices
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
load(file ='matrixsamples_stanmodel_ASJP_confounded.rda')
posterior_draws = as.data.frame(smpls)


vals = colMeans(posterior_draws)
nlangs = length(incllist)
incllist = gsub('.*\\.', '', incllist)
counter = 1
mymatrix_cnf = matrix(0, nrow = nlangs, ncol = nlangs,dimnames = list(incllist,incllist))
for(i in 1:nlangs){
  for(e in 1:nlangs){
    if(e>i){
      mymatrix_cnf[i,e] = posterior_draws[1,counter]
      mymatrix_cnf[e,i] = mymatrix_cnf[i,e]
      counter = counter+1
    }
  }
}

load(file ='matrixsamples_stanmodel_ASJP.rda')
posterior_draws = as.data.frame(smpls)

vals = colMeans(posterior_draws)
nlangs = length(incllist)
incllist = gsub('.*\\.', '', incllist)
counter = 1
mymatrix = matrix(0, nrow = nlangs, ncol = nlangs,dimnames = list(incllist,incllist))
for(i in 1:nlangs){
  for(e in 1:nlangs){
    if(e>i){
      mymatrix[i,e] = posterior_draws[1,counter]
      mymatrix[e,i] = mymatrix[i,e]
      counter = counter+1
    }
  }
}
mx = max(mymatrix_cnf)
dmat = mymatrix_cnf/max(mymatrix_cnf)-mymatrix/max(mymatrix)


dmat_dup = c(dmat[order(dmat)][!duplicated(dmat[order(dmat)])][1:5], dmat[order(dmat, decreasing = T)][!duplicated(dmat[order(dmat, decreasing = T)])][1:5])
dmat_dup2 = as.data.frame(dmat_dup)
dmat_dup2$pair = NA
for(i in 1:nrow(dmat_dup2)){
  dmat_dup2$pair[i] = paste(rownames(which(dmat == dmat_dup2$dmat_dup[i],arr.ind =T)), collapse = ', ')
}
xtable::xtable(dmat_dup2)








