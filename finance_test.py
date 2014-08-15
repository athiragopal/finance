__author__ = 'josh'

import unittest
import finance


class TestIntegrationAlgorithms(unittest.TestCase):
    def setUp(self):
        self.trapezoid_example = finance.trapezoidal_rule(lambda x: x, 1, 2, 1)
        self.simpsons_example = finance.simpsons_rule(lambda x: x ** 3, 1, 2, 1)
        self.black_scholes_call = finance.black_scholes_european_call_price(-5, 100, 0.5, 0.2, 1)
        self.black_scholes_put = finance.black_scholes_european_put_price(-5, 100, 0.5, 0.2, 1)

    def test_integration_exact(self):
        # The trapezoidal rule is exact for linear integrands with one strip and Simpson's Rule is exact for polynomials
        # up to cubic order using a single strip.
        self.assertEqual(self.trapezoid_example, 1.5)
        self.assertEqual(self.simpsons_example, 3.75)

    def test_no_arbitrage(self):
        # The Black-Scholes results assume no arbitrage opportunities. If the spot price of the underlying is negative
        # the value of the option should be zero.
        self.assertEqual(self.black_scholes_call, None)
        self.assertEqual(self.black_scholes_put, None)


if __name__ == '__main__':
    unittest.main()
