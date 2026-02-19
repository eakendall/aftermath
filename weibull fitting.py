import numpy as np
import scipy
from scipy.stats import weibull_min
from scipy.optimize import minimize

# Parameters for truncation
a, b = 0, 540 - 2 * 17    # truncated domain [0, 506], corresponds to diagnosis by day 540 with 17*2 days of symptoms

# Target mean and median for the truncated distribution of time to symptom onset (after disaese course analysis.R, line 92)
target_mean = 184
target_median = 146
proportion_before_90 = 0.35
proportion_before_360 = 0.85

# Function to compute truncated mean and truncated median from weibull parameters
def truncated_statistics(params):
    # parameters for weibull distribution
    k, theta = params
    dist = weibull_min(c=k, scale=theta)

    # Compute the mean of the truncated distribution
    # mean_truncated = dist.mean() * (dist.cdf(b) - dist.cdf(a)) / (dist.cdf(b) - dist.cdf(a))
    mean_truncated = dist.expect(lb = a, ub = b, conditional= True)

    # Compute the median of the truncated distribution    
    p50 = 0.5 * (dist.cdf(b) - dist.cdf(a))  # Adjust for truncation
    median_truncated = dist.ppf(p50 + dist.cdf(a))

    # Compute the proportion before day 90
    p90 = dist.cdf(90)
    proportion_90 = (p90 - dist.cdf(a)) / (dist.cdf(b) - dist.cdf(a))

    # Compute the proportion before day 360
    p360 = dist.cdf(360)
    proportion_360 = (p360 - dist.cdf(a)) / (dist.cdf(b) - dist.cdf(a))
    
    return mean_truncated, median_truncated, proportion_90, proportion_360

# Objective function to minimize the difference between target and estimated mean/median
def objective(params):
    mean_est, median_est, p90_est, p360_est = truncated_statistics(params)
    
    # Compute squared errors for mean and median
    error_mean = (mean_est - target_mean) ** 2
    error_median = (median_est - target_median) ** 2
    error_proportion90 = ((p90_est - proportion_before_90)*1000) ** 2 # weight proportion error more
    error_proportion360 = ((p360_est - proportion_before_360)*1000) ** 2 # weight proportion error more
    
    return error_mean + error_median + error_proportion90 + error_proportion360

# Initial guess for shape (k) and scale (theta)
initial_guess = [0.75, 320]

# Fit the weifull distribution parameters
result = minimize(objective, initial_guess, bounds=[(0.1, None), (0.1, None)])

# Optimal shape and scale for weibull distribution
optimal_k, optimal_theta  = result.x

print("Optimal shape (k, c):", optimal_k)
print("Optimal scale (theta, lambda):", optimal_theta)

# Validate the fit
mean_truncated, median_truncated, p90, p360 = truncated_statistics((optimal_k, optimal_theta))
print("Fitted truncated mean:", mean_truncated)
print("Fitted truncated median:", median_truncated)
print("Fitted proportion before day 90:", p90)
print("Fitted proportion before day 360:", p360)