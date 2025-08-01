# CDNP
Conditional Discrete Non-Parametric tools for characterising and simulating streamflow errors.

This package contains tools for estimating daily streamflow errors for use during stochastic simulation. The probability of errors are assumed dependent on the previous days error (autocorrelation assumption) and the current days flow (heteroscedasticity assumption).
    The approach determines conditional probabilities by discretising the current day streamflow errors (e_t) into separate ranges or bins. Also, previous day errors (e_t-1) and current day streamflow values (q_t) are categorised using k-means clustering. The probability of each streamflow error bin is computed for each cluster. Thus, for any pair of e_t-1 and q_t values, we can determine the cluster it belongs to (by computing the shortest distance to cluster centroids). Given this, we can determine the probability of each current day streamflow error bin.
