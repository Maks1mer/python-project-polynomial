
from math import sqrt

class Polynomial:

    @staticmethod
    def Trim(coefs):
        while len(coefs) > 1 and coefs[-1] == 0:
            coefs.pop()
        return coefs

    def __init__(self, *args):
        if len(args) == 0:
            self.coefs = [0]
            return
        if len(args) == 1:
            if isinstance(args[0], list):
                self.coefs = args[0].copy()
                return
            if isinstance(args[0], dict):
                self.coefs = [0 for _ in range(max(args[0].keys()) + 1)]
                for deg, coef in args[0].items():
                    self.coefs[deg] = coef
                return
            if isinstance(args[0], Polynomial):
                self.coefs = args[0].coefs.copy()
                return
        self.coefs = Polynomial.Trim(list(args))

    def __repr__(self):
        return "Polynomial "+ self.coefs.__str__()

    def __iter__(self):
        self.iter_range = zip(range(len(self.coefs)), self.coefs)
        return self.iter_range.__iter__()

    def __next__(self):
        return self.iter_range.__next__()

    @staticmethod
    def TermToString(coef, deg, var='x'):
        if coef == 0:
            return ""
        output = "+ {0}" if coef > 0 else "- {0}"
        coef = abs(coef)
        if deg > 0:
            coef = "" if coef == 1 else coef
            output += var
        if deg > 1:
            output += "^{1}"
        return output.format(coef, deg) + " "

    def __str__(self):
        self.coefs = Polynomial.Trim(self.coefs)
        if len(self.coefs) == 0:
            return ""
        output = Polynomial.TermToString(self.coefs[-1], len(self.coefs) -1)
        output = output[0] + output[2:]
        if output.startswith("+"):
            output = output[1:]
        for deg in range(len(self.coefs) - 2, -1, -1):
            output += Polynomial.TermToString(self.coefs[deg], deg)
        return output.strip()

    def __eq__(self, other):
        return isinstance(other, Polynomial) and Polynomial.Trim(self.coefs) == Polynomial.Trim(other.coefs)

    def __add__(self, other):
        if isinstance(other, Polynomial):
            new_coefs = self.coefs.copy() if len(self.coefs) > len(other.coefs) else other.coefs.copy()
            for i in range(min(len(self.coefs), len(other.coefs))):
                new_coefs[i] = self.coefs[i] + other.coefs[i]
            res = Polynomial(Polynomial.Trim(new_coefs))
        else:
            res = Polynomial(self)
            res.coefs[0] += other
        return res

    __radd__ = __add__

    def __neg__(self):
        return Polynomial([-coef for coef in self.coefs])

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __call__(self, x):
        res = 0
        for deg in range(len(self.coefs)):
            res += self.coefs[deg] * x ** deg
        return res

    def degree(self):
        return len(Polynomial.Trim(self.coefs)) - 1

    def der(self, d=1):
        old_coefs = self.coefs
        for _ in range(d):
            new_coefs = list()
            for deg in range(1, len(old_coefs)):
                new_coefs.append(deg * old_coefs[deg])
            old_coefs = new_coefs
        return Polynomial(old_coefs)

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            new_coefs = [0 for _ in range(self.degree() + other.degree() + 1)]
            for deg1, coef1 in self:
                for deg2, coef2 in other:
                    new_coefs[deg1 + deg2] += coef1 * coef2
            return Polynomial(new_coefs)
        else:
            return Polynomial([coef * other for coef in self.coefs])

    __rmul__ = __mul__

    def __mod__(self, other):
        if isinstance(other, Polynomial):
            result = Polynomial(self)
            while result.degree() >= other.degree():
                result = result - Polynomial([0 for _ in range(result.degree() - other.degree())] + other.coefs) * (result.coefs[-1] / other.coefs[-1])
            return result
        else:
            return Polynomial([coef % other for coef in self.coefs])

    def __rmod__(self, other):
        if self.degree() == 0:
            return other % self.coefs[0]

    def gcd(self, other):
        while self != Polynomial() and other != Polynomial():
            if self.degree() > other.degree():
                self = self % other
            else:
                other = other % self
        return self + other

class RealPolynomial(Polynomial):
    def find_root(self):
        self.coefs = Polynomial.Trim(self.coefs)
        neg = -1
        pos = 1
        mid = 0
        if  self.coefs[-1] < 0:
            pos, neg = neg, pos
        while not self(neg) < 0 < self(pos):
            neg *= 2
            pos *= 2
        while abs(self(mid)) > 1e-6:
            if self(mid) > 0:
                pos = mid
            elif self(mid) < 0:
                neg = mid
            else:
                break
            mid = (pos + neg) / 2
        return mid

class QuadraticPolynomial(Polynomial):
    def solve(self):
        D = self.coefs[1] ** 2 - 4 * self.coefs[0] * self.coefs[2]
        if D > 0:
            return [(-self.coefs[1] - sqrt(D)) / (2 * self.coefs[2]), (-self.coefs[1] + sqrt(D)) / (2 * self.coefs[2])]
        if D == 0:
            return [-self.coefs[1] / (2 * self.coefs[2])]
        return []
