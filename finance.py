# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 21:52:19 2014

@author: josh
"""
import numpy
import math
import scipy.stats

def trapezoidal_rule(function, a, b, n):
    # The trapezoidal rule can be found in nay first year calculus textbook. 
    # See for example Stewart.
    answer = 0
    # calculate the strip size 
    h = float(b-a) / n
    # evaluate the endpoints
    answer += function(a) + function(b)
    # evaluate the midpoints
    answer += sum([2*function(a + k*h) for k in xrange(1,n)])
    answer *= (h/2) 
   
    return answer


def simpsons_rule(function, a, b, n):
    # Simpsons Rule can be found in any first year calculus textbook.  See for 
    # example Stewart.
    answer = 0
    # calculate the strip size
    h = (b-a) / (2*n)
    # evaluate the endpoints
    answer += function(a) + function(b)
    # evaluate the odd terms
    answer += sum([2*function(a + (2*k)*h) for k in xrange(1,n) ])
    # evaluate the even terms
    answer += sum([4*function(a + (2*k-1)*h) for k in xrange(1,n+1)])
    # multiply by the final factor
    answer *= h/3
    return answer


def crude_monte_carlo(function, a, b, n):
    # The crude Monte Carlo approximation uses a discretisation of the average 
    # value of a function formula from first year calculus.  Psuedo-random 
    # numbers are selected from the numpy uniform distribution special 
    # function.  a and b are the bounds of the integral and n is the number
    # of samples.
    total =[function(scipy.random.uniform(a, b)) for k in xrange(n)]
    return scipy.mean(total)


def stratified_monte_carlo(function, a, b, n):
    # The stratified monte carlo approximation is an improvement on the crude
    # algorithm by dividing the integration region into subintervals, which 
    # reduces the error in the approximation.
    midpoint = (b-a) / 2 
    total_1 = scipy.mean([function(numpy.random.uniform(a, midpoint)) 
                        for k in xrange(n)])
    total_2 = scipy.mean([function(numpy.random.uniform(midpoint, b)) 
                        for k in xrange(n)])
    
    return total_1 + total_2 

def black_scholes_european_call_price(spot, strike, r, sigma, maturity):
    # The Black-Scholes formalism can be found in any textbook on mathematical 
    # finance, for example Hull.  spot is the current value of the underlying 
    # asset. strike is the price at which one can buy the stock.
    # r is the risk free interest rate. sigma is the volatility.  The maturity 
    # is the time until the option expires.
    d_1 = ((math.log(spot / strike) + (r + sigma**2 / 2) * maturity) 
            / (sigma * math.sqrt(maturity)))
    d_2 = d_1 - sigma * math.sqrt(maturity)
    return (spot*scipy.stats.norm.cdf(d_1) - strike*math.exp(-r * maturity) 
            *scipy.stats.norm.cdf(d_2)) 
    

def black_scholes_european_put_price(spot, strike, r, sigma, maturity):
    # The Put price can be calulated from the call price using put-call parity. 
    # See for example Hull.
    d_1 = ((math.log(spot/strike) + (r + sigma**2 / 2)*maturity) 
            / (sigma * math.sqrt(maturity)))
    d_2 = d_1 - sigma*math.sqrt(maturity)
    return (strike*math.exp(-r * maturity)*scipy.stats.norm.cdf(-d_2) - spot 
            *scipy.stats.norm.cdf(-d_1))


def monte_carlo_european_call_price(spot, strike, r, sigma, maturity, n):
    # This is a monte carlo simulation of a european call by estimating the
    # discounted expected payoff of the option under a risk neutral 
    # probability. 
    S_factor = spot*math.exp(maturity*(r - 0.05*sigma**2))
    S_current = [S_factor*math.exp(math.sqrt(sigma**2*maturity)
                *numpy.random.normal(0,1)) for k in xrange(n)]
    total = [max(Asset - strike,0.0) for Asset in S_current]

    return scipy.mean(total) *  math.exp(-r*maturity)
    

def monte_carlo_european_put_price(spot, strike, r, sigma, maturity, n):
    # This is a monte carlo simulation of a european put by estimating the
    # discounted expected payoff of the option under a risk neutral 
    # probability. 
    S_factor = spot*math.exp(maturity*(r-0.05*sigma**2))
    S_current = [S_factor*math.exp(math.sqrt(sigma**2*maturity)
                *numpy.random.normal(0,1)) for k in xrange(n)]
    total = [max(strike - Asset,0.0) for Asset in S_current]
        
    return scipy.mean(total) *  math.exp(-r*maturity)

 
def black_scholes_european_call_3d(sigma, r, strike, maturity, asset_steps):
    ds = 2*strike / asset_steps
    dt = 0.9 / ( sigma**2 * asset_steps**2)
    time_steps = int(maturity/dt) + 1
    dt = maturity / time_steps
    v = numpy.zeros((asset_steps, time_steps))
    s = numpy.zeros(asset_steps)
    
    for i in range(asset_steps):
        s[i] = i * ds
        v[i, 0] = max(s[i] - strike, 0) 
    
    for k in range(1, time_steps):
        for i in range(1, len(s)-1):
            delta = (v[i+1,k-1]-v[i-1,k-1])/(2*ds)
            gamma = (v[i+1,k-1]-2*v[i,k-1]+v[i-1,k-1])/ds**2
            theta = -0.5*(sigma**2)*s[i]**2*gamma - r*s[i]*delta + r*v[i,k-1]
            v[i,k] = v[i,k-1] - dt * theta
            
    v[0,k] = v[0,k-1] * (1-r*dt)
    v[len(s)-1,k] = 2*v[len(s)-2, k] - v[len(s)-3, k]
    
    return v


def monte_carlo_asset_price_path(spot, mu, sigma, maturity,asset_paths, 
                                 time_steps):
    # This is a vectorised euler method discretisation of a GBM asset price 
    # path. See for example Hull.
    S = numpy.zeros((asset_paths, time_steps+1))
    dt = maturity/time_steps
    S[:,0] = spot
    
    epsilon = numpy.random.normal(0, 1, (asset_paths,time_steps))
    S[:,1:] = numpy.exp((mu-0.5*sigma**2)*dt + epsilon*sigma*numpy.sqrt(dt))
    S = numpy.cumprod(S, axis = 1)
    
    return S
    

    




    



    
    




