# lavaanStan

`lavaanStan` is a package designed to fit simple Bayesian SEM (*Structural Equation Modeling*) and CFA (*Confirmatory Factor Analysis*) written with `lavaan` syntax using precompiled `Stan` models.

## Usage

`lavaanStan` works just like `lavaan`. Define a model using `lavaan` intuitive syntax and call `lavaanStan` function identifying the model type (CFA or SEM).

```r
source(lavaanStan.R)

data("PoliticalDemocracy")
data <- PoliticalDemocracy
model <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
   # residual covariances
     y1 ~~ y5
     y2 ~~ y4 + y6
     y3 ~~ y7
     y4 ~~ y8
     y6 ~~ y8
'
fit <- lavaanStan(model=model, data=data, model.type='sem')
modelMatrix(fit, 'est')
```

## Features

* Multiple groups SEM and CFA;
* User-defined equality constraints (including across groups);
* Posterior predictive p-values based on simulated data $\chi^2$;
* Model comparison available using LOO and WAIC (via `loo` package).

## Not working yet

* Latent growth models;
* Models with observed exogenous variables;
* Categorical indicator variables;
* Inequality constraints;
* Custom priors (`lavaanStan` uses "default" weakly informative priors);
* ...

## Links

[`lavaan`](http://lavaan.ugent.be/) (News, tutorials and support for the `lavaan` package.)

[`blavaan`](https://github.com/ecmerkle/blavaan) (Another (much better and flexible) implementation of Bayesian SEM using `lavaan` syntax)

[`Stan`](http://mc-stan.org/) (If you are here, `Stan` requires no introductions)
