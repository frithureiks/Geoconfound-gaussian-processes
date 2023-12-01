library(rstan)

load(file = 'cog.rda')
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
  Dmat = geoDists,
  miss = as.numeric(is.na(chars))
)

dat_phylo$C[is.na(chars)] = -1


load(file ='fullsummary_readymodel_northeuralex.rda')
shortsumm = smry$summary
res = data.frame('true' = dat_phylo$C, 'pred' = NA)

res$r = inv_logit(shortsumm[grepl('^r\\[', rownames(shortsumm)),1])
res$k = inv_logit(shortsumm[grepl('^k\\[', rownames(shortsumm)),1])
p =  inv_logit(shortsumm[grepl('^r\\[', rownames(shortsumm)),1] + 
  shortsumm[grepl('^k\\[', rownames(shortsumm)),1])
res$r_raw = shortsumm[grepl('^r\\[', rownames(shortsumm)),1]
res$k_raw = shortsumm[grepl('^k\\[', rownames(shortsumm)),1]
res$k_raw_lower = shortsumm[grepl('^k\\[', rownames(shortsumm)),4]
res$k_lower = inv_logit(shortsumm[grepl('^k\\[', rownames(shortsumm)),4])
res$drk = res$r-res$k
res$drk_raw = res$r_raw-res$k_raw

library(rethinking)
res$pred = p

res$d = res$true-res$pred
hist(res$d)
sum(abs(res$d) < 0.3)/nrow(res)

res$disc = 0
res$disc[res$pred >0.5] = 1

caret::confusionMatrix(data = as.factor(res$disc), reference = as.factor(res$true))

sres = res[res$true == 1,]



df = read.csv('inferred-cognates.tsv', sep = '\t', header = F)
borr = read.csv('wold-nelex-intersect.tsv', sep = '\t', header = F)
langs = c('isl', 'eng', 'nor', 'nld',
          'dan', 'deu', 'smj', 'sma',
          'swe', 'sme', 'fin', 'ekk',
          'liv', 'lav', 'lit', 'pol',
          'bel', 'smn', 'sms', 'krl',
          'sjd', 'olo', 'vep', 'rus'
)



df$b = NA
for(i in 1:nrow(borr)){
  print(i/nrow(borr))
  i2 = strsplit(borr$V3[i], split = ':')
  i3 = strsplit(borr$V4[i], split = ':')
  if (i2[[1]][1] %in% langs & i3[[1]][1] %in% langs){
    l = i3[[1]][1]
    w = i3[[1]][2]
    df[df$V2 == l & df$V3 == w,'b'] = 1
  }
}

save(df, file = 'dellert_cognates2.rda')

load(file = 'dellert_cognates2.rda')

newdf2 = data.frame()

for(i in 1:length(unique(df$V1))){
  print(i/length(unique(df$V1)))
  subdf = df[df$V1 == unique(df$V1)[i],]
  tmpdf = data.frame(matrix(ncol = length(unique(df$V2)), nrow = length(unique(subdf$V5))))
  colnames(tmpdf) <- unique(df$V2)
  for (row in 1:nrow(tmpdf)) {
    yesin = unique(subdf$V2[subdf$V5 == row])
    noin = unique(subdf$V2[subdf$V5 != row])
    tmpdf[row,yesin] = i
    tmpdf[row,setdiff(noin,yesin)] = NA
  }
  newdf2 = rbind(newdf2, tmpdf)
  
}
save(newdf2, file = 'conceptdf.rda')
load(file = 'conceptdf.rda')


newdf3 = data.frame()

for(i in 1:length(unique(df$V1))){
  print(i/length(unique(df$V1)))
  subdf = df[df$V1 == unique(df$V1)[i],]
  tmpdf = data.frame(matrix(ncol = length(unique(df$V2)), nrow = length(unique(subdf$V5))))
  colnames(tmpdf) <- unique(df$V2)
  for (row in 1:nrow(tmpdf)) {
    if(1 %in% subdf[subdf$V5==row,'b']){
      l = subdf[which(subdf$V5==row & subdf$b == 1),'V2']
      tmpdf[row,l] = 1
    }

  }
  newdf3 = rbind(newdf3, tmpdf)
  
}
save(newdf3, file = 'borrdf.rda')
load(file = 'borrdf.rda')


load(file = 'cog.rda')
newdf = newdf[,c('isl', 'eng', 'nor', 'nld',
                 'dan', 'deu', 'smj', 'sma',
                 'swe', 'sme', 'fin', 'ekk',
                 'liv', 'lav', 'lit', 'pol',
                 'bel', 'smn', 'sms', 'krl',
                 'sjd', 'olo', 'vep', 'rus'
                 
                 
)]
d1 = which(rowSums(newdf) != 0)
newdf = newdf[d1,]
d2 = which(rowSums(newdf) != 1)
newdf = newdf[d2,]
d3 = which(rowSums(newdf) != ncol(newdf))
newdf = newdf[d3,]
d4 = which(rowSums(newdf) != ncol(newdf))
newdf = newdf[d4,]


newdf3 = newdf3[,c('isl', 'eng', 'nor', 'nld',
                 'dan', 'deu', 'smj', 'sma',
                 'swe', 'sme', 'fin', 'ekk',
                 'liv', 'lav', 'lit', 'pol',
                 'bel', 'smn', 'sms', 'krl',
                 'sjd', 'olo', 'vep', 'rus'
)]
newdf3 = newdf3[d1,]
newdf3 = newdf3[d2,]
newdf3 = newdf3[d3,]
newdf3 = newdf3[d4,]

dmtrx = matrix(res$k_lower, nrow = nrow(newdf3), ncol = ncol(newdf3), byrow = T)
highgeo = res[which(res$k_lower > 0.4 & res$true == 1), 'k']

dmtrx2 = matrix(0, nrow = nrow(newdf3), ncol = ncol(newdf3), byrow = T)
for(i in 1:length(highgeo)){
  dmtrx2[which(dmtrx == highgeo[i], arr.ind = T)] = 1
}

library(mltools)
library('MLmetrics')
newdf3[is.na(newdf3)] = 0
mltools::auc_roc(c(dmtrx2), c(as.matrix(newdf3)))
F1_Score(c(as.matrix(newdf3)),c(dmtrx2))
ConfusionMatrix(c(dmtrx2), c(as.matrix(newdf3)))
mcc(c(dmtrx2), c(as.matrix(newdf3)))

F1_Score(c(as.matrix(newdf3)),c(dmtrx2), positive = '1')*(sum(c(as.matrix(newdf3)))/length(c(as.matrix(newdf3)))) +
F1_Score(c(as.matrix(newdf3)),c(dmtrx2), positive = '0')*((length(c(as.matrix(newdf3)))-sum(c(as.matrix(newdf3))))/length(c(as.matrix(newdf3))))

###2ndapp
dmtrx = matrix(res$k_raw, nrow = nrow(newdf3), ncol = ncol(newdf3), byrow = T)
ddf = data.frame('lik' = c(dmtrx), 'true' = c(as.matrix(newdf3)))

ddf$true[ddf$true == 0] = 'no'
ddf$true[ddf$true == 1] = 'yes'
ddf$true = as.factor(ddf$true)

g <- ggplot(ddf, aes(x=true, y=lik, fill = true)) + geom_boxplot(notch = T,outlier.size = 0.2,outlier.alpha = 0.2) + 
  xlab('Loan?')+ ylab('Log-likelihood') +
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

cairo_pdf("loanboxplot.pdf", width = 6, height = 4)

g

dev.off()

library(brms)

mod = brm(lik ~ true, data= ddf, prior = set_prior('normal(0,1)', class = 'b')+set_prior('exponential(1)', class = 'sigma'))

mod

ggplot(ddf, aes(x=true, y=inv_logit(lik), fill = true)) + geom_boxplot(notch = T) + 
  xlab('Loan?')+ ylab('Log-likelihood') +
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )

#####Confoundstrength


load(file = 'cog.rda')
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
  Dmat = geoDists,
  miss = as.numeric(is.na(chars))
)

dat_phylo$C[is.na(chars)] = -1

load(file ='fullsummary_readymodel_northeuralex.rda')
shortsumm = smry$summary

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


#c1_vals = smry$summary[grepl('Dmat2_vec', rownames(smry$summary)),1]
#c1_m = matrconv(c1_vals, length(langvec))
#dimnames(c1_m) <-list(langvec,langvec)
#my_upgma <- phangorn::upgma(as.dist(c1_m))


#cairo_pdf("tree_northeuralex.pdf", width = 12, height = 8)
#
#plot(my_upgma)
#
#dev.off()



Nlang = length(langvec)

K <- matrix(0,nrow=Nlang,ncol=Nlang)
for ( i in 1:Nlang ){
  print(i/Nlang)
  for ( j in 1:Nlang){
    K[i,j] <- smry$summary[grepl('etasq$', rownames(smry$summary)),1] * exp( -(dat_phylo$Dmat[i,j]^2 / smry$summary[grepl('rhosq', rownames(smry$summary)),1]) )
  }}
diag(K) <- smry$summary[grepl('etasq$', rownames(smry$summary)),1] + 0.01

# convert to correlation matrix
Rho <- cov2cor(K)
# add row/col names for convenience
colnames(Rho) <- langvec
rownames(Rho) <- colnames(Rho)

crstr = Rho
crstr[!upper.tri(Rho)] = 0
highest = order(crstr, decreasing = T)[1:10]
values = sort(crstr, decreasing = T)[1:10]

crstrdf = data.frame(
  'strength' = values,
  'L1' = NA,
  'L2' = NA
)

for (i in 1:length(values)){
  crstrdf$L1[i] = colnames(crstr)[which(crstr == values[i], arr.ind = T)[1]]
  crstrdf$L2[i] = colnames(crstr)[which(crstr == values[i], arr.ind = T)[2]]
}

# plot raw data and labels

cairo_pdf("geoconfound_strength_northeuralex.pdf", width = 12, height = 8)

plot( NULL, xlab="Longitude" , ylab="Latitude" ,
      col=col.alpha("black",0.8) , cex=1 , pch=1, xlim=c(-30,60), ylim=c(50,70), lwd = 3  )
labels <- langvec

world <- map_data("world", xlim=c(-30,60), ylim=c(50,70))

for(g in unique(world$group)){
  world2 = world[world$group == g,]
  polygon(world2[,1:2],
          col = NULL, border = col.alpha("black",0.3))
}

for( i in 1:Nlang ){
  for ( j in 1:Nlang ){
    if ( i < j ){
      lines( c( geodat$lon[i],geodat$lon[j] ) , c( geodat$lat[i],geodat$lat[j] ) ,
             #lwd=2 , col=col.alpha("black",Rho[i,j]^8 /max(Rho[upper.tri(Rho^8)]^8))) 
             lwd=2 , col=col.alpha("blue",Rho[i,j]^12)) 
    }}}
points( geodat$lon , geodat$lat,
        col=col.alpha("black",1) , cex=1 , pch=1, lwd = 3  )
#text( geodat$lon , geodat$lat, labels=labels , cex=0.8 , pos=1)




dev.off()



cairo_pdf("geoconfound_strength_northeuralex2.pdf", width = 12, height = 8)

plot( NULL, xlab="Longitude" , ylab="Latitude" ,
      col=col.alpha("black",0.8) , cex=1 , pch=1, xlim=c(-30,60), ylim=c(50,70), lwd = 3  )
labels <- langvec

#world <- map_data("world", xlim=c(-30,60), ylim=c(50,70))

#for(g in unique(world$group)){
#  world2 = world[world$group == g,]
#  polygon(world2[,1:2],
#          col = NULL, border = col.alpha("black",0.3))
#}

for( i in 1:Nlang ){
  for ( j in 1:Nlang ){
    if ( i < j ){
      lines( c( geodat$lon[i],geodat$lon[j] ) , c( geodat$lat[i],geodat$lat[j] ) ,
             #lwd=2 , col=col.alpha("black",Rho[i,j]^8 /max(Rho[upper.tri(Rho^8)]^8))) 
             lwd=2 , col=col.alpha("blue",Rho[i,j]^12)) 
    }}}
points( geodat$lon , geodat$lat,
        col=col.alpha("black",1) , cex=1 , pch=1, lwd = 3  )
text( geodat$lon , geodat$lat, labels=labels , cex=0.8 , 
      pos=c(1,1,2,1,2,1,2,2,3,3,3,4,2,4,4,1,1,3,3,4,4,4,4,1))




dev.off()


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

#langvec = geodat$name[order(geodat$iso_code)]
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


#library(phangorn)
#cns1 = consensus(posterior_trees, p=.5, rooted = T)
#plot(cns1)
#nodelabels(round(cns1$node.label, digits = 2),adj=c(1,-0.2),frame="none")
#cns2 = ls.consensus(posterior_trees2)
#plot(cns2)
#nodelabels(round(cns2$node.label, digits = 2),adj=c(1,-0.2),frame="none")

cns = consensus(posterior_trees2, p=.5, rooted = T)
boot <- prop.clades(cns, posterior_trees2, rooted = T)
cns$node.label = round(boot/8000, digits = 2)

cairo_pdf("tree_northeu.pdf", width = 8, height = 6)
ggtree(cns) + geom_tiplab() + xlim(0, 8) + geom_nodelab(hjust=-.1, size=3)
dev.off()


mcc = maxCladeCred(posterior_trees2)
#plot(mcc)
#nodelabels(round(mcc$node.label, digits = 2),adj=c(1,-0.2),frame="none")

#mcc2 = maxCladeCred(posterior_trees)
#plot(mcc2)
#nodelabels(round(mcc2$node.label, digits = 2),adj=c(1,-0.2),frame="none")


load('matrixsamples_stanmodel_northeuralex_confounded.rda')
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

#langvec = geodat$name[order(geodat$iso_code)]
incllist = geodat$name
nlangs = length(incllist)

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



#plot(cns)
#nodelabels(round(cns$node.label, digits = 2),adj=c(1,-0.2),frame="none")
cns2 = consensus(posterior_trees2_cnf, p=.5, rooted = T)
#plot(cns2)
#nodelabels(round(cns2$node.label, digits = 2),adj=c(1,-0.2),frame="none")

boot <- prop.clades(cns2, posterior_trees2_cnf, rooted = T)
cns2$node.label = round(boot/8000, digits = 2)

cairo_pdf("tree_northeu_confounded.pdf", width = 8, height = 6)
ggtree(cns2) + geom_tiplab() + xlim(0, 8) + geom_nodelab(hjust=-.1, size=3)
dev.off()



#plot(mcc)
#nodelabels(round(mcc$node.label, digits = 2),adj=c(1,-0.2),frame="none")
#mcc2 = maxCladeCred(posterior_trees2_cnf)
#plot(mcc2)
#nodelabels(round(mcc2$node.label, digits = 2),adj=c(1,-0.2),frame="none")

#mcc2 = maxCladeCred(posterior_trees)
#plot(mcc2)
#nodelabels(round(mcc2$node.label, digits = 2),adj=c(1,-0.2),frame="none")

obj<-cophylo(cns2,cns, rotate = T)


cairo_pdf("treecompare_northeu.pdf", width = 8, height = 6)
plot(obj,fsize=c(0.8,0.8),link.type="curved",link.lwd=2, lty = 'dashed', ylim = c(0,1), align = T)
dev.off()


#partcheck

prevcount = 0
counter = 0
for(tree in 1:length(posterior_trees2)){
  pp_test = prop.part(posterior_trees2[[tree]])
  for (part1 in 1:23) {
    if (all(pp_test[[part1]] == c(1,  2,  3,  4,  5,  6,  9, 14, 15, 16, 17, 24))){
      for (part2 in 1:23) {
        if (all(pp_test[[part2]] == c(1, 2, 3, 4, 5, 6, 9))){
          prevcount = prevcount +1
          for (part3 in 1:23) {
            if (all(pp_test[[part3]] == c(2, 4, 6))){
              counter = counter +1
              break
            }
            
          }
          
        }
      
      }
    }

    
    
  }
}

#pp = prop.clades(cns, part = prop.part(posterior_trees2), rooted = T)
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


head(pp_comp[order(pp_comp$d),], n= 20)

xtable::xtable(head(pp_comp[order(pp_comp$d),], n= 5))


head(pp_comp[order(pp_comp$d, decreasing = T),], n= 20)
xtable::xtable(head(pp_comp[order(pp_comp$d, decreasing = T),], n= 5))

####matrices

load(file ='matrixsamples_stanmodel_northeuralex_confounded.rda')
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

load(file ='matrixsamples_stanmodel_northeuralex.rda')
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
dmat = mymatrix_cnf/sd(mymatrix_cnf)-mymatrix/sd(mymatrix)

dmat_dup = c(dmat[order(dmat)][!duplicated(dmat[order(dmat)])][1:5], dmat[order(dmat, decreasing = T)][!duplicated(dmat[order(dmat, decreasing = T)])][1:5])
dmat_dup2 = as.data.frame(dmat_dup)
dmat_dup2$pair = NA
for(i in 1:nrow(dmat_dup2)){
  dmat_dup2$pair[i] = paste(rownames(which(dmat == dmat_dup2$dmat_dup[i],arr.ind =T)), collapse = ', ')
}
xtable::xtable(dmat_dup2)




