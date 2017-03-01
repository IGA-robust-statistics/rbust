q_scale = function(data, c=2.21914){
  # Return robust estimate of scale proposed by
  # Rosseeuw, Croux: Time-Efficient Algorithms for
  # Two Highly robust Estimators of Scale
  #
  # Input:
  #   data - array, sample observations or time series
  #   c - float, consistency factor, default for Gaussian Distribution
  #
  # Return:
  #   Q_n - float, estimated scale

  n = length(data)
  k = as.integer(((n * (n - 1) / 2) + 2) / 4 ) + 1
  diffs = sort(apply(combn(data,2), 2, function(x) abs(x[1] - x[2])))
  Q_n = c * diffs[k]
  return(Q_n)
}

q_acf = function(ts, lag.max=5, cor=F){
  # Return robust estimate of autocovariance function proposed by
  # Ma, Genton: Highly Robust Estimation of the Autocovariance function
  #
  # Input:
  #   ts - array, time series
  #   max.lag - int, maximum time lag for which acf is calculated
  #   cor - bool, it true return autocorrelation, else autocovariance (default)
  #
  # Return:
  #   acf - float, estimated robust autocovariance/autocorrelation

  if (!is.numeric(lag.max) || lag.max < 0)
    stop("'lag.max' must be integer > 0")
  acf = rep(0, lag.max)
  n = length(ts)
  for(lag in 1:lag.max){
    u = ts[1:(n-lag)]
    v = ts[(lag + 1):n]
    Q_plus = q_scale(u+v)^2
    Q_minus = q_scale(u-v)^2
    if(cor){
      acf[lag] = (Q_plus - Q_minus) / (Q_plus + Q_minus)
    } else {
      acf[lag] = (Q_plus - Q_minus) / 4
    }
  }
  return(acf)
}

# Compare with std estimate
ts = arima.sim(model=list(ar=c(.9, -0.1)),n=100)
robust_acf = q_acf(ts, lag.max=20, cor=T)
std_acf = acf(ts)$acf[2:21]

plot(robust_acf, std_acf)
abline(0, 1)
