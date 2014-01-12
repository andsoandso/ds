import unittest
import numpy as np


class DiscreteTests(unittest.TestCase):
    self.fns = [lambda x: (2.5*x)*(1-x),
        lambda x: x**3 + 1,
        lambda x, c1, c2: np.sqrt(c1/(x+c2))]
    
    self.seeds = [-10, -1, 0.5, 0, 0.5, 1, 10]

    # TODO

if __name__ == "__main__":
    unittest.main()
