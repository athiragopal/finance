# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 21:52:19 2014

@author: josh
"""
import numpy
import math
import scipy.stats


def trapezoidal_rule(function, a, b, n):
    # The trapezoidal rule is a technique for approximating a definite integral by estimating the area under the graph
    # of the integrand as a trapezoid.  For more information see http://en.wikipedia.org/wiki/Trapezoidal_rule. a and
    # b are the integration bounds, and n is an integer number of strips to use.
    answer = 0
    # calculate the strip size 
    h = float(b - a) / n
    assert type(n) == int
    # evaluate the endpoints
    answer += function(a) + function(b)
    # evaluate the midpoints
    answer += sum([2 * function(a + k * h) for k in xrange(1, n)])
    answer *= (h / 2)
    return answer


def simpsons_rule(function, a, b, n):
    # Simpson's rule is a method for approximating a definite integral by replacing the integrand with a polynomial
    # that takes the same value as the integrand at the endpoints and the midpoint. This function is an implementation
    # of the composite Simpson's rule. See http://en.wikipedia.org/wiki/Simpson's_rule for more information. a and b are
    # the integration bounds, and n is the integer number of strips to use in the approximation.
    answer = 0
    # calculate the strip size
    h = float(b - a) / (2 * n)
    assert type(n) == int
    # evaluate the endpoints
    answer += function(a) + function(b)
    # evaluate the odd terms
    answer += sum([2 * function(a + (2 * k) * h) for k in xrange(1, n)])
    # evaluate the even terms
    answer += sum([4 * function(a + (2 * k - 1) * h) for k in xrange(1, n + 1)])
    # multiply by the final factor
    answer *= h / 3
    return answer


def crude_monte_carlo(function, a, b, n):
    # The crude Monte Carlo approximation uses a discretisation of the average 
    # value of a function formula from first year calculus.  Pseudo-random
    # numbers are selected from a uniform distribution. For more information see
    # http://en.wikipedia.org/wiki/Monte_Carlo_integration. a and b are the bounds of the integral and n is the integer
    # number of samples.
    assert type(n) == int
    total = sum([function(scipy.random.uniform(a, b)) for k in xrange(n)])
    answer = ( total / n) * (b - a)
    return answer


def black_scholes_european_call_price(spot, strike, r, sigma, maturity):
    # The Black-Scholes option pricing formula provides an exact formula for the fair value of a european vanilla
    # call option.  For a review of the Black-Scholes model see
    # http://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model.  spot is the spot price of the underlying asset, strike
    # is the strike price, r is the risk free interest rate, sigma is the volatility of the underlying asset, and
    # maturity is the time to maturity.
    if spot < 0:  # No arbitrage condition
        return None
    else:
        d_1 = ((math.log(spot / strike) + (r + sigma ** 2 / 2) * maturity)
               / (sigma * math.sqrt(maturity)))
        d_2 = d_1 - sigma * math.sqrt(maturity)
        answer = (spot * scipy.stats.norm.cdf(d_1) - strike * math.exp(-r * maturity)
                  * scipy.stats.norm.cdf(d_2))
        return answer


def black_scholes_european_put_price(spot, strike, r, sigma, maturity):
    # The Put price can be calculated from the call price using put-call parity. More information may be found at
    # http://en.wikipedia.org/wiki/Put%E2%80%93call_parity.
    if spot < 0:  # No arbitrage condition
        return None
    else:
        d_1 = ((math.log(spot / strike) + (r + sigma ** 2 / 2) * maturity)
               / (sigma * math.sqrt(maturity)))
        d_2 = d_1 - sigma * math.sqrt(maturity)
        answer = (strike * math.exp(-r * maturity) * scipy.stats.norm.cdf(-d_2) - spot
                  * scipy.stats.norm.cdf(-d_1))
        return answer


def monte_carlo_european_call_price(spot, strike, r, sigma, maturity, n):
    # This is a monte carlo simulation of a european call by estimating the
    # discounted expected payoff of the option under a risk neutral 
    # probability. See http://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model for more information.
    if spot < 0:  # No arbitrage condition
        return None
    else:
        s_factor = spot * math.exp(maturity * (r - 0.05 * sigma ** 2))
        s_current = [s_factor * math.exp(math.sqrt(sigma ** 2 * maturity)
                                         * numpy.random.normal(0, 1)) for k in xrange(n)]
        total = [max(Asset - strike, 0.0) for Asset in s_current]
        answer = scipy.mean(total) * math.exp(-r * maturity)
        return answer


def monte_carlo_european_put_price(spot, strike, r, sigma, maturity, n):
    # This is a monte carlo simulation of a european put by estimating the
    # discounted expected payoff of the option under a risk neutral 
    # probability. See http://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model for more information.
    if spot < 0:  # No arbitrage condition
        return None
    else:
        s_factor = spot * math.exp(maturity * (r - 0.05 * sigma ** 2))
        s_current = [s_factor * math.exp(math.sqrt(sigma ** 2 * maturity)
                                         * numpy.random.normal(0, 1)) for k in xrange(n)]
        total = [max(strike - Asset, 0.0) for Asset in s_current]
        answer = scipy.mean(total) * math.exp(-r * maturity)
        return answer


def black_scholes_european_call_3d(sigma, r, strike, maturity, asset_steps):
    ds = 2 * strike / asset_steps
    dt = 0.9 / (sigma ** 2 * asset_steps ** 2)
    time_steps = int(maturity / dt) + 1
    dt = maturity / time_steps
    v = numpy.zeros((asset_steps, time_steps))
    s = numpy.zeros(asset_steps)

    for i in range(asset_steps):
        s[i] = i * ds
        v[i, 0] = max(s[i] - strike, 0)

    for k in range(1, time_steps):
        for i in range(1, len(s) - 1):
            delta = (v[i + 1, k - 1] - v[i - 1, k - 1]) / (2 * ds)
            gamma = (v[i + 1, k - 1] - 2 * v[i, k - 1] + v[i - 1, k - 1]) / ds ** 2
            theta = -0.5 * (sigma ** 2) * s[i] ** 2 * gamma - r * s[i] * delta + r * v[i, k - 1]
            v[i, k] = v[i, k - 1] - dt * theta

        v[0, k] = v[0, k - 1] * (1 - r * dt)
        v[len(s) - 1, k] = 2 * v[len(s) - 2, k] - v[len(s) - 3, k]

    return v


# noinspection PyNoneFunctionAssignment
def monte_carlo_asset_price_path(spot, mu, sigma, time_horizon, asset_paths,
                                 time_steps):
    # This function simulates the price of an asset undergoing a geometric brownian motion using a vectorised euler
    # discretisation. See http://en.wikipedia.org/wiki/Geometric_Brownian_motion for a description of this process.
    # spot is the spot price of the asset, mu is the drift, sigma is the volatility, T is the time horizon, asset_paths
    # is the number of asset asset paths, and time_steps is the number of time steps.
    s = numpy.zeros((asset_paths, time_steps + 1))
    dt = time_horizon / time_steps
    s[:, 0] = spot
    epsilon = numpy.random.normal(0, 1, (asset_paths, time_steps))
    s[:, 1:] = numpy.exp((mu - 0.5 * sigma ** 2) * dt + epsilon * sigma * numpy.sqrt(dt))
    s = numpy.cumprod(s, axis=1)
    return s

print black_scholes_european_call_price(100, 100, 0.05, 0.02, 1)


