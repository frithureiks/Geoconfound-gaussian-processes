gridsearch_res = read.csv('gridsearch_res.csv', header = F)

colnames(gridsearch_res) <- c('nlangs',
                              'nsites',
                              'covariance',
                              'rhat_mean',
                              'rhat_phymat_mean',
                              'splitfreq_mean',
                              'splitfreq_majority_mean',
                              'majority_clades')

gridsearch_res$splitfreq_majority_mean[gridsearch_res$majority_clades == 1] = NA

cor(gridsearch_res$rhat_mean, gridsearch_res$rhat_phymat_mean)

mean(gridsearch_res$splitfreq_mean[gridsearch_res$majority_clades != 1] - gridsearch_res$splitfreq_majority_mean[gridsearch_res$majority_clades != 1])

mean(gridsearch_res$splitfreq_mean[gridsearch_res$majority_clades != 1] - 
       gridsearch_res$splitfreq_majority_mean[gridsearch_res$majority_clades != 1])/mean(gridsearch_res$splitfreq_mean[gridsearch_res$majority_clades != 1])


mod_rhat = lm(rhat_phymat_mean ~ (nlangs + nsites + covariance)^2, data = gridsearch_res)
smry_rhat = summary(mod_rhat)
mod_splitfreq = lm(splitfreq_mean ~ (nlangs + nsites + covariance)^2, data = gridsearch_res)
smry_splitfreq = summary(mod_splitfreq)
mod_splitfreq2 = lm(splitfreq_majority_mean ~ (nlangs + nsites)^2, data = gridsearch_res)
smry_splitfreq2 = summary(mod_splitfreq2)

