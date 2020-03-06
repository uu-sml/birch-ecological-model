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
birds between 1976 and 1998 (extracted from [Saether et al., "Estimating the time to Extinction
in an island population of song sparrows", Proc. Biol. Sci., Vol. 267, No. 1443
(2000), pp. 621-626][2]).

## Particle Marginal Metropolis Hastings
The model has a `proposal` function that can be use together with the
`ParticleMarginalMetropolisHasting` sampler in the `metropolis-hastings` branch
of the [Birch standard library][1]. Note that you need to use [commit 55a047c][3] of the
Birch language to run the simulation.


[1]:https://github.com/lawmurray/Birch.Standard/tree/metropolis-hastings
[2]:https://www.jstor.org/stable/2665985
[3]:https://github.com/lawmurray/Birch/commit/55a047c6670c750c50b6bd55d832e7c04d4867b5
