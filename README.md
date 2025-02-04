# nowcasting


# Notes about modifications


- Function: `generate_exponential_q`
    * It is not taking into account an intercept.
    * It is accumulating the probabilities of reporting two times!
- Model: fixed_q
    * alpha and beta can not be estimated, one value should be fixed
    * set a uniform distribution on the simplex vector
- Gamma distribution:
    * Is useful for theoretical results, but is very unstable for sampling.
- Model: fixed_b
    * the priors should respect the support of the parameter
- Exponential curve:
    * Values greater than 4 are not identifiable. Place strong prior with low chance
      of values greater than 4.
- Random walk (Use centered random walks):
    * Do not set-up initial values, it will lead to higher variance through time.
    Impose a condition instead such as the mean is 0.
- Scenario random walk:
    * Inadequate, we end-up with full zeros, meaning that cases are not reported. Can
      not estimate or predict on those cases.
