/**
 * Forward model. To simulate the model, first use `start()`, followed by
 * `step()` any number of times.
 */
class ForwardModel < Model {   
  /**
   * Start.
   */
  function start() -> Real;

  /**
   * Step.
   */
  function step() -> Real;
  
  /**
   * Rewind to start.
   */
  function rewind();

  /**
   * Proposal function.
   */
  function propose(x:ForwardModel) -> (Real, Real);

  /**
   * Number of steps.
   */
  function size() -> Integer;
}
