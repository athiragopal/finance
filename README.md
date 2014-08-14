# Finance


Finance is a python module for calculations in quantitative finance. It includes
* Numerical integration algorithms
* An implementation of the Black-Scholes European option pricing results and numerical approximations
* A geometric brownian motion simulator

# Using Finance

1. [Getting Finance](https://github.com/joshuadebellis/finance)
2. [Installation](# Installation)
3. [Finance Tutorial](# Finance Tutorial)
  1. [numerical integration](## numerical integration)
  2. [Black-Scholes option pricing](## Black-Scholes option pricing)
  3. [geometric brownian motion simulator](## geometric brownian motion simulator) 

# Installation
Finance can be installed by placing finance.py in the same directory as the main python script.

## numerical integration

Finance includes four one-dimensional numerical integration algorithms. It includes the trapezium rule, Simpson's Rule, a
monte-carlo integration rule, and a stratified monte-carlo algorithm such that the error in the approximation is reduced.These are implemented as `trapezoidal_rule(function, a, b, n)`, `simpsons_rule(function, a, b, n)`,`crude_monte_carlo(function, a, b, n)`, and `stratified_monte_carlo(function, a, b, n)`.

Simpson's Rule and the trapezium rule may be used in the following way.
```
print trapezoidal_rule(lambda x: x**2, a, b, n)
>>> 0.33375
```
where a and b are integer bounds of integration, and n is the number of midpoints and trapezoids, respectively.

The monte-carlo integration algorithms are invoked in a similar way.
```
print crude_monte_carlo(lambda x: x**2, a, b, n)
>>> 0.372776750963
```
```
print stratified_monte_carlo(lambda x: x**2, a, b, n)
>>> 0.382230814658
```
where a and b are integer bounds of integration, and n is the number of pseudo-random samples to take.







## Black-Scholes option pricing

Finance includes several functions related to european call and put option pricing. The first is an exact calculation of the vanilla call and put option prices. These is implemented in Finance as `black_scholes_european_call_price(spot, strike, r, sigma, maturity)` and `black_scholes_european_put_price(spot, strike, r, sigma, maturity)`. `spot` is the
spot price of the option, `strike` is the price at which the holder of an option may exercise, `r` is the risk-free 
interest rate, `sigma` is the volatility, and `maturity`, is the time at which the contract expires.

The call and put pricing functions are called in an identical way.
```
print black_scholes_european_call_price(spot=100, strike=100, r=0.05, sigma0.02, maturity=1)
>>> 4.88096669701
```

Finance provides monte-carlo estimations of the put and call prices by estimating the discounted expected payoff at
maturity under a risk-neutral probability. In Finance these functions are `monte_carlo_european_call_price(spot, strike, r, sigma, maturity, n)` and `monte_carlo_european_put_price(spot, strike, r, sigma, maturity, n)`. Here `n` is the 
number of simulations to perform. For example, the european call estimator is invoked as
```
print monte_carlo_european_call_price(spot=100, strike=100, r=0.05, sigma0.02, maturity=1, n=1000)
>>> 4.91577454629
```
A finite difference approximator to the Black-Scholes equation for a european call is included in finance as `black_scholes_european_call_3d(sigma, r, strike, maturity, asset_steps)`. It returns an array containing option prices and time steps. The time steps are hardcoded to be as large as possible whilst ensuring stability.



## Geometric Brownian Motion Simulator

An vectorised euler algorithm for simulating option prices using geometric brownian motion is included in Finance as `monte_carlo_asset_price_path(spot, mu, sigma, maturity,asset_paths, time_steps)`. `spot` is the spot price of the asset, `mu` is the drift, `sigma` is the volatility, `maturity`, is the time at which the contract expires, `asset_paths` is the number of underlying asset paths to simulated, and `time_steps` is the number of time steps in the discretisation. It returns an array `S` consisitng of spot price and time step data. For example,

```simulation = finance.monte_carlo_asset_price_path(spot=100, mu=0.05, sigma=0.3, maturity=1.0, asset_paths=1000, time_steps=50)```

Produces an array `S` which when plotted produces a graph
![Screenshot](http://imgur.com/AAGb8DH.png)



