class ParticleMarginalMetropolisHastings < ConditionalParticleFilter {
  w':Real;


  function start() {
    if !h'? {
      /* no previous particle, run normal particle filter */
      super.start();
      w' <- -inf; // Always accept first sample
    } else {
      /* There is a reference particle. We create a proposal particle
       * by cloning the reference particle and the event handler
       */


       auto x' <- this.x[b];
       auto h <- x[b].h;
       h.setMode(PLAY_DELAY);
       auto forwardCount <- h.trace.forwardCount;
       auto forward <- h.trace.takeForward();

       auto x <- clone<ForwardModel>(this.x[b]);
       auto f <- clone<StackNode<Event>>(forward!);
       x.h.setMode(REPLAY_DELAY);
       x.h.trace.putForward(f, forwardCount);

      q':Real;
      q:Real;
       (q', q) <- x.propose(x');

       auto w <- x.play();

       if (log(simulate_uniform(0.0, 1.0)) < w + q' - w' - q) {
         x' <- clone<ForwardModel>(x);
         w' <- w;
       }

       this.x[b] <- x';
       h <- this.x[b].h;
       h.setMode(REPLAY_DELAY);
       f <- clone<StackNode<Event>>(forward!);
       cpp{{
         f.finish();
       }}
       h.trace.putForward(f, forwardCount);
    }
  }
}


 //      auto x' <- this.x[b]; // Previous model
 //      auto forwardCount <- x'.h.trace.forwardCount;
 //      auto forward <- x'.h.trace.takeForward();

 //      q:Real; // log q(θ | θ') (new given old)
 //      q':Real; // log q(θ' | θ) (old given new)

 //      auto x <- clone<ForwardModel>(x'); // Proposal model
 //      x.h.setMode(REPLAY_DELAY);
 //      // copy the trace from the reference into the proposal particle
 //      x.h.trace.forward <- clone<StackNode<Event>>(forward!);

 //      (q', q) <- x.propose(x'); //  Propose parameters from reference

 //      auto w <- x.play(); // Compute p(y|x,θ)p(x|θ)p(θ)

 //      if (log(simulate_uniform(0.0, 1.0)) < w + q' - w' - q) {
 //        x' <- clone<ForwardModel>(x);
 //        w' <- w;
 //      }

 //      this.x[b] <- x';
 //      auto h <- this.x[b].h;
 //      h.setMode(REPLAY_DELAY);
 //      h.trace.putForward(forward, forwardCount);
