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
    * Impose a condition instead such as the mean is 0.
    * Under log and logistic transformation, only small variances are identifiable.
    * Log: Standard deviation higher than 0.2 explodes.
    * Logistic: Standard deviation higher than 0.2 cover all possible ranges already.
- Scenario random walk:
    * Inadequate, we end-up with full zeros, meaning that cases are not reported. Can
      not estimate or predict on those cases.
- Priors:
    * Lambda: should be a sensible value the allows high number of cases, but can
    quickly move toward 0:
    * Moderate: rlnorm(100000, meanlog = 0, sdlog = 2)
    * Allows very high numbers: rlnorm(100000, meanlog = 0, sdlog = 3)
