data("PoliticalDemocracy")
dados <- PoliticalDemocracy[, c(9:11, 1:8)]

model <- '
   # latent variables
     
     y1  ~  x1 + x2 + x3
     y5 ~ y1 + y7 + x3
   # 
'
fit <- sem(model,
           data=PoliticalDemocracy, meanstructure = T, fixed.x=F)
stanFit <- stan('lavaanPath.stan', data=buildDataList(fit, dados, rep(1, dim(dados)[1])),
                iter = 1000, warmup=500, chains=4, thin=2, control = list(adapt_delta=0.85))
                #pars=c('Alpha', 'Nu', 'Lambda', 'Psi', 'PHI', 'PPP', 'phi'),
                #init=lapply(1:4, function(x) initf(fit)))
