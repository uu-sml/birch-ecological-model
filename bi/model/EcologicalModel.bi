class EcoState {
  L:Random<Real>;

  function write(buffer:Buffer) {
    buffer.set("L", L);
  }
}

class EcoParameter {
  b:Random<Real[_]>;
  c:Random<Real>;
  Q:Random<Real>;
  R:Random<Real>;

  function write(buffer:Buffer) {
    buffer.set("b", b);
    buffer.set("c", c);
    buffer.set("Q", Q);
    buffer.set("R", R);
  }
}

class EcologicalModel < HMMWithProposal<EcoParameter, EcoState, Random<Real>> {

  μ:Real <- 0.01;

  σ2:Real <- 1.01;

  α_Q:Real <- 2.01;
  β_Q:Real <- 10.01;

  α_R:Real <- 2.01;
  β_R:Real <- 10.01;

  fiber parameter(θ:EcoParameter) -> Event {
    θ.c ~ Normal(μ, σ2);
  }

  fiber initial(x:EcoState, θ:EcoParameter) -> Event {
    x.L ~ Gaussian(40.0, 10.0);
    θ.Q ~ InverseGamma(α_Q, β_Q);

    θ.b ~ Gaussian(vector(μ, 2), identity(2)*θ.Q);
    θ.R ~ InverseGamma(α_R, β_R);
  }

  fiber transition(x':EcoState, x:EcoState, θ:EcoParameter) -> Event {
    auto e <- exp(x.L*θ.c);
    auto A <- [1.0, e.value()]; 
    x'.L ~ Gaussian(dot(A, θ.b) + x.L , θ.Q);
  }

  fiber observation(y:Random<Real>, x:EcoState, θ:EcoParameter) -> Event {
    y ~ Gaussian(x.L.value(), θ.R); // force realization to prevent pruning
                                    // of the delayed graph
  }

  function propose(x:ForwardModel) -> (Real, Real) {

    auto x_old <- (HMMWithProposal<EcoParameter,EcoState,Random<Real>>?(x))!;

    auto θ_old <- x_old.θ; // Parameter from previous model
    
    auto σ2 <- 5.5;
    auto Q_old <- Normal(θ_old.c, σ2); // q(θ' | θ)

    θ.c <- Q_old.simulate(); // Draw new parameter for this model
    auto q <- Q_old.observe(θ.c);  // log q(θ | θ') (new given old)

    auto Q_new <- Normal(θ.c, σ2); // q(θ' | θ) 
    auto q_old <- Q_new.observe(θ_old.c); // log q(θ | θ') (old given new)

    return (q, q_old);
  }
}
