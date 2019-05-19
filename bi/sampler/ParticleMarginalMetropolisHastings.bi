class ParticleMarginalMetropolisHastings < ParticleFilter {

  m:ForwardModel?;
  py:Real;
  q:Real;
  m':ForwardModel;
  py':Real;
  q':Real;
  A:Real <- 0;
  R:Real <- 0;


  function initialize() {
    if !m? {
      // No pervious model, run as normal
      super.initialize();
    } else {
      /* There is a previous model, create particles from  proposal model */
      py <- sum(Z.walk()); // Save previous evidence
      
      m' <- clone<ForwardModel>(archetype!);
      (q', q) <- m'.propose(m!);

      Z.clear();
      ess.clear();
      memory.clear();
      elapsed.clear();

      w <- vector(0.0, N);
      a <- iota(1, N);
      o <- vector(1, N);

      x1:Vector<ForwardModel>;
      x1.enlarge(N, clone<ForwardModel>(m'));
      x <- x1.toArray();

      parallel for auto n in 1..N {
        x[n] <- clone<ForwardModel>(x[n]);
      }

      tic();
    }
  }
  
  function finalize() {
    // The filter has drawn x[b] and the associated 
    if m? {
      m' <- x[b];
      py' <- sum(Z.walk());
      if (log(simulate_uniform(0.0, 1.0)) < py' + q - py - q') {
        // accept and save for next iteration
        m <- clone<ForwardModel>(m');
        py <- py';
        A <- A + 1;
      } else {
        // reject and replace drawn model with previous model
        R <- R + 1;
        x[b] <- clone<ForwardModel>(m!);
      }
    } else {
      m <- clone<ForwardModel>(x[b]);
    }
    if verbose {
      stderr.print("acceptance rate: " + A/(A+R) + "\n");
    }
  }

}