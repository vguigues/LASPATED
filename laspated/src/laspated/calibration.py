import numpy as np


class Param(object):
    def __init__(self) -> None:
        self.EPS = 1e-5
        self.sigma = 0.5
        self.gap = 0.01
        self.max_iter = 30
        self.lower_lambda = 1e-6
        self.upper_lambda = 1.0
        self.beta_bar = 2.0
        self.cv_proportion = 0.2


class RegularizedModel(object):
    def __init__(
        self,
        N,
        M,
        durations,
        groups,
        weights,
        alphas,
        distance,
        type_region,
        neighbors,
        param,
    ):
        if len(N.shape) != 3:
            raise ValueError(f"N dimension must be 3: C,R,T, but got {len(N.shape)}")
        if len(M.shape) != 3:
            raise ValueError(f"M dimension must be 3: C,R,T, but got {len(M.shape)}")

        self.C, self.R, self.T = N.shape
        self.nb_observations = N
        self.nb_arrivals = M
        self.durations = durations
        self.alpha = alphas
        self.weights = weights
        self.distance = distance
        self.neighbors = neighbors
        self.type_region = type_region
        self.groups = groups
        self.nb_groups = len(groups)
        self.param = param

        self.which_group = np.zeros((self.T))
        for g in range(self.nb_groups):
            for i in groups[g]:
                self.which_group[i] = g
                
        self.emp_rates_by_class = np.zeros((self.C))
        for c in range(self.C):
            sum_c = 0.0
            for r in range(self.R):
                for t in range(self.T):
                    sum_c += self.nb_arrivals[c,r,t] / self.nb_observations
                    
        

    def f(self, x):
        obj = 0
        nb_observations = self.nb_observations
        nb_arrivals = self.nb_arrivals
        durations = self.durations
        alpha = self.alpha
        distance = self.distance

        for c in range(self.C):
            for r in range(self.R):
                for t in range(self.T):
                    obj += nb_observations[c, r, t] * x[c, r, t] * durations[
                        t
                    ] - nb_arrivals[t] * np.log(x[c, r, t])
                    for s in self.neighbors[r]:
                        obj += (
                            0.5
                            * alpha[r, s]
                            * nb_observations[c, r, t]
                            * nb_observations[c, s, t]
                            * ((x[c, r, t] - x[c, s, t]) ** 2)
                            / distance[r, s]
                        )

        for c in range(self.C):
            for r in range(self.R):
                for grindex in range(self.nb_groups):
                    group = self.groups[grindex]
                    for t in group:
                        for t1 in group:
                            if t != t1:
                                obj += (
                                    0.5
                                    * self.weights[grindex]
                                    * nb_observations[c, r, t]
                                    * nb_observations[c, r, t1]
                                    * (x[c, r, t] - x[c, r, t1]) ** 2
                                )
        return obj

    def gradient(self, x):
        grad = np.zeros(x.shape)

        return grad

    def projection(self, x):
        z = np.zeros(x.shape)
        for c in range(self.C):
                for r in range(self.R):
                    for t in range(self.T):
                        z = max(self.param.lower_lambda, x[c,r,t])
        if not self.param.relax_emp_fix:
            for c in range(self.C):
                sum_c = 0.0
                for r in range(self.R):
                    for t in range(self.T):
                        sum_c += z[c,r,t]*self.durations[t]
                diff = self.emp_rates_by_class[c] - sum_c
                if s
        return z
