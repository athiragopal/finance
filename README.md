# Finance


Finance is a python module for calculations in mathematical finance. It includes
* Numerical integration algorithms
* An implementation of the Black-Scholes european option pricing results and numerical approximation algorithms
* A geometric brownian motion simulator

# Using Finance

1. [Getting Finance](https://github.com/joshuadebellis/finance)
2. [Installation](# Installation)
3. [Finance Tutorial](# Finance Tutorial)
  1. [numerical integration](## numerical integration)
  2. [Black-Scholes option pricing](## Black-Scholes option pricing)
  3. [geometric brownian motion simulator](## geometric brownian motion simulator) 

# Installation
Finance can be installed by placing finance.py in the same directory as the main python script.  This module depends on
on both the [Numpy](http://www.numpy.org/) and [Scipy](http://www.scipy.org/) packages.

## numerical integration

Finance includes four one-dimensional numerical integration algorithms. It includes the trapezoidal rule, Simpson's Rule, and
monte-carlo integration algorithm. These are implemented as `trapezoidal_rule(function, a, b, n)`, `simpsons_rule(function, a, b, n)`,
and`crude_monte_carlo(function, a, b, n)`, 

Simpson's Rule and the trapezoidal rule may be used in the following way.
```
>>> print finance.trapezoidal_rule(lambda x:x**2, 0, 1, 20)
>>> 0.33375
```
where:
* `function` is the integrand of the definite integral
* `a` and `b` are the bounds of integration
* `n` is an integer number of strips to use in the approximations.

It returns a `float` which approximates the area under the integrand.

The monte-carlo integration algorithm is invoked in a similar way.
```
>>> print finance.crude_monte_carlo(lambda x:x**2, 0, 1, 1000)
>>> 0.332749640773
```

where: 
* `function` is the integrand of the definite integral
* `a` and `b` are the bounds of integration
* `n` is an integer number of pseudo-random samples to take

It returns a `float` which approximates the area under the integrand.

## Black-Scholes option pricing

Finance includes several functions related to european vanilla call and put option pricing. The first is an exact 
calculation of the vanilla call and put option prices. These is implemented in Finance as 
`black_scholes_european_call_price(spot, strike, r, sigma, maturity)` and `black_scholes_european_put_price(spot, 
strike, r, sigma, maturity)`.  

The call and put pricing functions are called in an identical way.
```
>>> print finance.black_scholes_european_call_price(100, 100, 0.05, 0.02, 1)
>>> 4.88096669701
```
where:
* `spot` is the spot price of the underlying asset
* `strike` is the price at which the holder of an option may exercise
* `r` is the risk-free interest rate
* `sigma` is the volatility of the underlying asset
* `maturity`, is the time at which the contract expires.

It returns a `float` which is the analytical solution to the Black-Scholes partial differential equation. `None` is
returned in the case of arbitrage opportunities.


Finance provides monte-carlo estimations of the put and call prices by estimating the discounted expected payoff at
maturity under a risk-neutral probability. In Finance these functions are `monte_carlo_european_call_price(spot, strike,
 r, sigma, maturity, n)` and `monte_carlo_european_put_price(spot, strike, r, sigma, maturity, n)`. For example, the european call estimator is invoked as
```
>>> print finance.monte_carlo_european_call_price(100, 100, 0.05, 0.02, 1, 1000)
4.75598195465
```

where:
* `n` is the integer number of simulations to perform

It returns a `float` which estimates the analytical result of the Black-Scholes pricing formula.




A finite difference approximation scheme to the Black-Scholes equation for a european call is included in finance as 
`black_scholes_european_call_3d(sigma, r, strike, maturity, asset_steps)`. 



It returns an `array` containing spot price,time step data, and the corresponding option value data. 

The time steps are hardcoded to be as large as possible whilst ensuring stability. This function can be used to study the option surface. For example using the 
[matplotlib](http://matplotlib.org/) package we can create the call option surface with the following code:
```
import finance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



result = finance.black_scholes_european_call_3d(0.2, 0.05, 100, 1.0, 20)
x = [k * 0.1 for k in xrange(18)]
y = [k * 10 for k in xrange(20)]
fig = plt.figure()
ax = Axes3D(fig)

x, y = numpy.meshgrid(x, y)
ax.plot_surface(x, y, result, rstride=1, cstride=1, cmap=plt.cm.hot)

ax.set_xlabel('Time')
ax.set_ylabel('Asset')
ax.set_zlabel('Option Value')

plt.show()

```
This results in
![Screenshot](http://imgur.com/lZ3TJxj.png)




## geometric brownian motion simulator

A vectorised euler algorithm for simulating asset prices using geometric brownian motion is included in Finance as
 `monte_carlo_asset_price_path(spot, mu, sigma, time_horizon, asset_paths, time_steps)`.

where:
* `asset_paths` is the number of asset paths to simulated
* `time_steps` is the number of time steps in the discretisation
* `time_horizon` is the time duration of the simulation

It returns an `array` containing time step data, and asset price data.





which when plotted using the following code
```
import finance
import pylab

simulation = finance.monte_carlo_asset_price_path(spot=100, mu=0.05, sigma=0.3, time_horizon=1.0, asset_paths=1000, time_steps=50)


pylab.plot(numpy.linspace(0,1,50+1), numpy.percentile(simulation, 95, axis = 0))
pylab.title('asset price')
pylab.xlabel('time')
pylab.show()

```
produces
![Screenshot](http://imgur.com/IegmzUc.png)



