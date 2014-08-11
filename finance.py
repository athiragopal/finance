# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 21:52:19 2014

@author: josh
"""
import numpy as np
import math
import scipy.stats

def Trapezoidal_Rule(function, a, b,n):
    answer = 0
    """calculate the strip size """
    h = float((b-a)) / n
    """ evaluate the endpoints """
    answer += function(a) + function(b)
    """evaluate the midpoints"""
    for k in range(1,n):
        answer += 2*function(a + k*h)
    """multiply by the final factor"""    
    answer *= h/2
    return answer

def Simpsons_Rule(function,a,b,n):
    answer = 0
    """calculate the strip size"""
    h = float((b-a)) / (2*n)
    """ evaluate the endpoints """
    answer += function(a) + function(b)
    """evaluate the odd terms"""
    for k in range(1,n):
        answer += 2*function(a+(2*k)*h)
    """evaluate the even terms"""
    for k in range(1,n+1):
        answer += 4*function(a+(2*k-1)*h)
    """multiply by the final factor"""
    answer *= h/3
    return answer

def Crude_Monte_Carlo(function, a,b,n):
    total = 0 
    range = b-a 
    counter = 1 
    
    while (counter <= n):
        random =np.random.uniform (a, b) 
        answer = function(random)
        total += answer 
        counter = counter+1 
        
    return (total / n) * range

def MC_stratified(function, a, b, n):
    """ find midpoint """
    midpoint = (b-a) / 2 
    
    total_1 = 0 
    total_2 = 0 
    
    range_1 = midpoint - a 
    range_2 = b - midpoint 
    
    counter = 1 ;
    
    while (counter <= n):
        random =np.random.uniform (a, b) 
        answer = function(random)
        total_1 += answer
        counter = counter + 1 
        
    result_1 =  (total_1 / n) * range_1;
    counter = 1

    while (counter <= n):
        random =np.random.uniform (a, b) 
        answer = function(random)
        total_2 += answer 
        counter = counter + 1
    result_2 =  (total_2 / n) * range_2
    return result_1 + result_2 

def gaussian_twopoint(function, a,b):
    x_1 = 0.57735502691896257 
    x_2 = -0.5773502691896257 
    f_1 = float(b-a )
    f_2 = float(b+a )
    answer = float(function((f_1*x_1+f_2)/2) + function((f_1*x_2+f_2)/2) )
    result = (f_1/2) * answer 
    return result 

def Black_Scholes_European_Call_Price(spot,strike,r,sigma,maturity ):
    d_1 = (math.log(spot / strike) + (r + sigma**2 / 2) * maturity) / (sigma * math.sqrt(maturity))
    d_2 = d_1 - sigma*math.sqrt(maturity)
    return spot * scipy.stats.norm.cdf(d_1) - strike * math.exp(-r * maturity) * scipy.stats.norm.cdf(d_2)   
    
def Black_Scholes_European_Put_Price(spot,strike,r,sigma,maturity ):
     d_1 = (math.log(spot / strike) + (r + sigma**2 / 2) * maturity) / (sigma * math.sqrt(maturity))
     d_2 = d_1 - sigma*math.sqrt(maturity)
     return spot * scipy.stats.norm.cdf(-d_1) - strike * math.exp(-r * maturity) * scipy.stats.norm.cdf(-d_2)
     
   
    
S0 = 100.0
K = 110.0
r = 0.03
sigma = 0.1
T = 0.5

result  = Black_Scholes_European_Call_Price(S0,K,r,sigma,T )
print result

def Monte_Carlo_European_Call_Price(spot,strike,r,sigma,maturity,n):
    S_factor = spot*math.exp(maturity*(r-0.05*sigma**2))
    S_current = 0.0
    total = 0.0
    counter = 0
    a,b = 0, 1
    for counter in range (0, n):
        random =np.random.normal(a,b) 
        S_current = S_factor * math.exp(math.sqrt(sigma*sigma*maturity)*random)
        total += max(S_current-strike,0.0)

    return ((total / n)* math.exp(-r*maturity))


result = Monte_Carlo_European_Call_Price(S0,K,r,sigma,T,5000000) 
print result   
    
