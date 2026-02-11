# CDNP
Conditional Discrete Non-Parametric tools for characterising and simulating streamflow errors.

This package contains tools for estimating daily streamflow errors for use during stochastic simulation. The probability of errors are assumed dependent on the previous days error (autocorrelation assumption), the current days flow (heteroscedasticity assumption), and optionally another upstream site's error (spatial dependence assumption).
    The main approach categorises previous day errors (e_t-1) and current day streamflow values (q_t) using k-means clustering. Each cluster represents a separate empirical distribution. Thus, for any pair of e_t-1 and q_t values, we can determine the relevant empirical distribution by computing which cluster centroid is closest to the pair of values.
