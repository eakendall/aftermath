import numpy as np
import scipy
from scipy.stats import gamma
from scipy.optimize import minimize

# Parameters for truncation
a, b = 0, 540 - 2 * 17   # truncated domain [0, 506]

# Target mean and median for the truncated distribution
target_mean = 192
target_median = 149

# Function to compute truncated mean and truncated median from gamma parameters
def truncated_statistics(params):
    k, theta = params  # shape (k) and scale (theta)
    dist = gamma(k, scale=theta)
    
    # Compute the mean of the truncated distribution
    mean_truncated = dist.mean() * (dist.cdf(b) - dist.cdf(a)) / (dist.cdf(b) - dist.cdf(a))
    
    # Compute the median of the truncated distribution
    p50 = 0.5 * (dist.cdf(b) - dist.cdf(a))  # Adjust for truncation
    median_truncated = dist.ppf(p50 + dist.cdf(a))
    
    return mean_truncated, median_truncated

# Objective function to minimize the difference between target and estimated mean/median
def objective(params):
    mean_est, median_est = truncated_statistics(params)
    
    # Compute squared errors for mean and median
    error_mean = (mean_est - target_mean) ** 2
    error_median = (median_est - target_median) ** 2
    
    return error_mean + error_median

# Initial guess for shape (k) and scale (theta)
initial_guess = [2, 2]

# Fit the gamma distribution parameters
result = minimize(objective, initial_guess, bounds=[(0.1, None), (0.1, None)])

# Optimal shape and scale for gamma distribution
optimal_k, optimal_theta = result.x

print("Optimal shape (k):", optimal_k)
print("Optimal scale (theta):", optimal_theta)

# Validate the fit
mean_truncated, median_truncated = truncated_statistics((optimal_k, optimal_theta))
print("Fitted truncated mean:", mean_truncated)
print("Fitted truncated median:", median_truncated)