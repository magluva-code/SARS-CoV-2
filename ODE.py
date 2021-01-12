import numpy as np


class State(object):
    def __init__(self, init_pop, n, init_state="S"):
        self.pop = init_pop
        self.state_array = np.zeros(n)
        self.state_array[0] = self.pop
        self.state_history = {0 : init_state}

    def __call__(self):
        return self.get_state()

    def __str__(self):
        return f'symptomatic: {self.symptomatic}, infectious: {self.infectious}, immunized: {self.immunized}'

    def __repr__(self):
        return f'Susceptible(symptomatic={self.symptomatic}, infectious={self.infectious}, immunized={self.immunized})'

    def advance_pop(self, new_pop):
        self.pop = new_pop

    def get_pop(self):
        return self.pop

    def advance_state(self, t):
        self.state_array[t] = self.pop

    def get_state(self):
        return self.state_array

    def update_state_history(self, state):
        self.state_history[t] = state

    def get_state_history(self):
        return self.state_history


class ODERK4(object):
    def __init__(self, f):
        self.f = f

    def __str__(self):
        pass

    def __repr__(self):
        pass

    def __call__(self):
        pass

    def set_init_cnd(self, u0):
        u0 = np.asarray(u0)
        self.size_u = np.shape(u0)[0]
        self.u0 = u0

    def solve(self, tp):
        self.t = np.asarray(tp)
        n = self.t.size
        self.u = np.empty((n, self.size_u), dtype=object)
        self.u[0] = self.u0
        for k in range(n-1):
            self.k = k
            self.u[k+1] = self.advance()
        return self.u[:k+2], self.t[:k+2]

    def advance(self):
        u, f, k ,t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        K1 = dt*f(u[k], t[k])
        K2 = dt*f(u[k] + (1/2)*K1, t[k] + (1/2)*dt)
        K3 = dt*f(u[k] + (1/2)*K2, t[k] + (1/2)*dt)
        K4 = dt*f(u[k] + K3, t[k] + dt)
        _u = u[k] + (1/6)*(K1 + 2*K2 + 2*K3 + K4)
        return _u


class ProblemSIDRV(object):
    def __init__(self, S, I, D, R, V):
        self.S = S
        self.I = I
        self.D = D
        self.R = R
        self.V = V

    def init_params(self, alpha, beta, gamma, delta):
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta

    def __call__(self, u, t):
        S, I, D, R, V = u
        _s = -self.beta(t)*self.S.pop*self.I.pop - self.gamma(t)*self.S.pop
        _i = self.beta(t)*self.S.pop*self.I.pop - self.delta(t)*self.I.pop - self.alpha(t)*self.I.pop
        _d = self.delta(t)*self.I.pop
        _r = self.alpha(t)*self.I.pop
        _v = self.gamma(t)*self.S.pop
        return np.array([_s, _i, _d, _r, _v])



class SolverSIDRV(object):
    def __init__(self, problem, T, dt):
        self.problem = problem
        self.T = T
        self.dt = dt

    def solve(self, method=ODERK4):
        n = int(self.T/self.dt)
        t = np.linspace(0, self.T, n, True)
        solver = method(self.problem)
        solver.set_init_cnd([self.problem.S, self.problem.I, self.problem.D, self.problem.R, self.problem.V])
        u, self.t = solver.solve(t)
        self.S, self.I, self.D, self.R, self.V = u[:, 0], u[:, 1], u[:, 2], u[:, 3], u[:, 4]


class SARSCoV2(object):
    def __init__(self):
        pass


def main():
    S0 = 1500
    I0 = 1
    D0 = 0
    R0 = 0
    V0 = 0
    T = 60
    dt = 0.5
    n = int(T/0.5)

    S = State(S0, n)
    I = State(I0, n, init_state="I")
    D = State(D0, n, init_state="D")
    R = State(R0, n, init_state="R")
    V = State(V0, n, init_state="V")
    problem = ProblemSIDRV(S, I, D, R, V)
    solver = SolverSIDRV(problem, T, dt)
    solver.solve()


if __name__ == "__main__":
    main()
