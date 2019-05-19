# Ecological Model
Birch implementation of the ecological model.
```
  x[t] = x[t-1] + b0 + b2 exp(b3 x[t-1]) + w[t]
  y[t] = x[t] + v[t]
```
Where `w[t] ~ N(0, Q)` and `v[t] ~ N(0, R)`.

We use conjugate priors for `b0`, `b2`, `Q`, and `R`. The parameter `b3` has a
normal prior:
```
  R ~ InverseGamma(α_R, β_R)
  Q ~ InverseGamma(α_Q, β_Q)
  [b0, b2] ~ Gaussian(μ, S*Q)
  b3 ~ Gaussian(μ3, s2)
  x[t] ~ Gaussian(x[t-1] + b0 + b2 exp(b3 x[t-1]), Q)
  y[t] ~ Gaussian(x[t], R)
```

Build and install
```
birch build
birch install
```

The file `input/observation_data.json` contains data from observations of song
birds between 1976 and 1998 (extracted from Saether et al., "Estimating the time to Extinction
in an island population of song sparrows", Proc. Biol. Sci., Vol. 267, No. 1443
(2000), pp. 621-626).

## Particle Marginal Metropolis Hastings
The code in `sampler/ParticleMarginalMetropolisHastings.bi` is experimental and is currently not working.
To compile it, add `model/ForwardModel.bi` to the standard library and recompile it. There is a problem with the sampler, however, and the proposal function is not called correctly (see [this issue][1]).

[1]:https://github.com/uu-sml/birch-ecological-model/issues/1
