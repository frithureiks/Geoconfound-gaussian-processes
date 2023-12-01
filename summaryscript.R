library(rstan)


load('stanmodel_ASJP.rda')
smry = summary(mod, pars = c('r', 'k', 'Dmat2_vec', 'rhosq', 'etasq', 'etasq2'))
save(smry, file ='fullsummary_stanmodel_ASJP.rda')
smpls = extract(mod, pars = c('Dmat2_vec'))
save(smpls, file ='matrixsamples_stanmodel_ASJP.rda')

load('stanmodel_ASJP_confounded.rda')
smry = summary(mod, pars = c('r', 'Dmat2_vec', 'etasq2'))
save(smry, file ='fullsummary_stanmodel_ASJP_confounded.rda')
smpls = extract(mod, pars = c('Dmat2_vec'))
save(smpls, file ='matrixsamples_stanmodel_ASJP_confounded.rda')

load('readymodel_northeuralex.rda')
smry = summary(mod, pars = c('r', 'k', 'Dmat2_vec', 'rhosq', 'etasq', 'etasq2'))
save(smry, file ='fullsummary_readymodel_northeuralex.rda')
smpls = extract(mod, pars = c('Dmat2_vec'))
save(smpls, file ='matrixsamples_stanmodel_northeuralex.rda')


load('stanmodel_northeuralex_confounded.rda')
smry = summary(mod, pars = c('r', 'Dmat2_vec', 'etasq2'))
save(smry, file ='fullsummary_stanmodel_northeuralex_confounded.rda')
smpls = extract(mod, pars = c('Dmat2_vec'))
save(smpls, file ='matrixsamples_stanmodel_northeuralex_confounded.rda')



