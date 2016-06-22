import numpy as np
import logging

logger = logging.getLogger(__name__)

class ScorePolynomial():
    def __init__(self, params, dist):
        self.cutoff, self.dim = params
        self.dist = dist

    def __str__(self):
        return ("Score parameters: normal, params: " +
                str((self.cutoff, self.dim)))

    def calculate(self, multiplier, a, b):
        return multiplier / (max(self.dist(a, b), self.cutoff)**self.dim)

class ScoreExp():
    def __init__(self, params, dist):
        self.cutoff, self.scale = params
        self.dist = dist

    def __str__(self):
        return ("Score parameters: exp, params: " +
                str((self.cutoff, self.scale)))

    def calculate(self, multiplier, a, b):
        return (multiplier / 
                max(np.exp(self.dist(a, b) / self.scale), self.cutoff))

class ScoreGauss():
    def __init__(self, params, dist):
        self.sigma = params
        self.dist = dist
        self.constant = 1/((2 * np.pi * self.sigma**2) ** 0.5)

    def __str__(self):
        return ("Score parameters: Gaussian, params: " +
                str(self.sigma))

    def calculate(self, multiplier, a, b):
        return  (self.constant * multiplier *
                 np.exp(-self.dist(a, b)**2 / (2 * self.sigma**2)))

def set_score_func(args, dist):
        
    type = args.score_type
    if type == 'normal':
        return ScorePolynomial((args.cutoff, args.dim), dist)

    elif type == 'exp':
        return ScoreExp((args.cutoff, args.scale), dist)

    elif type == 'Gaussian':
        return ScoreGauss(args.sigma, dist)

    else:
        logger.error("Invalid score_type in set_score_params!")
        quit()

