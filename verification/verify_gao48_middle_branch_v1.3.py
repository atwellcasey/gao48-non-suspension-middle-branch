#!/usr/bin/env python3
"""Exact symbolic verification for Atwell, Gao middle-branch non-suspension witness v1.3."""
import sympy as sp

x, y, t, gamma = sp.symbols('x y t gamma')
i = sp.I

# Witness
G2 = y*t**2
G3 = 2*x*y*t**2 - t
G4 = sp.Rational(4, 3)*i*y**3*t**3
H = sp.expand(G3 - 2*x*G2)

# Gao characteristic field for c=1, phi=0 is W = d/dy.
WG2 = sp.diff(G2, y)
assert sp.expand(WG2 - t**2) == 0

# Gao coefficients (equation (12) setup)
D2 = -sp.diff(H, t)
D1 = sp.expand(
    sp.diff(H, y)*sp.diff(G4, t)
    - sp.diff(H, t)*sp.diff(G4, y)
    - (sp.diff(G2, x)*sp.diff(G3, t) - sp.diff(G3, x)*sp.diff(G2, t))
)
D0 = sp.expand(sp.Matrix([
    [sp.diff(F, v) for v in (x, y, t)] for F in (G2, G3, G4)
]).det())

assert sp.expand(D2 - 1) == 0
assert sp.expand(D1 - 4*(1+i)*y**2*t**3) == 0
assert sp.expand(D0 - 8*i*y**4*t**6) == 0
assert sp.expand(4*D0 - D1**2) == 0

# Reconstructed threefold
X1 = sp.expand(-D1/2)
X2 = sp.expand(G2 + x*X1)
X3 = sp.expand(G3 + x**2*X1)
X4 = sp.expand(G4 + y*X1)
X = (X1, X2, X3, X4)

Delta = (1, x, x**2, y)
S = tuple(sp.expand(Xj + gamma*dj) for Xj, dj in zip(X, Delta))
JS = sp.Matrix([[sp.diff(F, v) for v in (gamma, x, y, t)] for F in S])
assert sp.expand(JS.det() - gamma**2) == 0

JX123 = sp.Matrix([[sp.diff(F, v) for v in (x, y, t)] for F in (X1, X2, X3)])
minor = sp.factor(JX123.det())
assert sp.expand(minor - (-4+12*i)*y**3*t**6) == 0

# Elimination certificate
R = sp.expand(X1*X3 - X2**2)
kappa = 2*(1+i)
rho = 1+2*i
d = -2-sp.Rational(2,3)*i
lam = sp.simplify(rho**3*d**2/kappa**6)
expected_lam = sp.Rational(41,576) - sp.Rational(19,288)*i
assert sp.simplify(lam - expected_lam) == 0
assert sp.expand(R - rho*y**2*t**4) == 0
assert sp.expand(R**3*X4**2 - lam*X1**6) == 0

# Singular-locus factor display for normalized lambda=1 model
z1,z2,z3,z4 = sp.symbols('z1 z2 z3 z4')
Rz = z1*z3-z2**2
F = sp.expand(Rz**3*z4**2-z1**6)
grad = [sp.factor(sp.diff(F,z)) for z in (z1,z2,z3,z4)]
assert grad[3] == 2*z4*Rz**3

# Check P1 and P2 are contained in the singular locus exactly by substitution.
for expr in [F] + grad:
    assert sp.expand(expr.subs({z1:0,z2:0})) == 0  # P1
    assert sp.expand(expr.subs({z1:0,z4:0})) == 0  # P2

print('PASS: exact witness verification')
print('W(G2) =', sp.factor(WG2))
print('D2 =', sp.factor(D2))
print('D1 =', sp.factor(D1))
print('D0 =', sp.factor(D0))
print('det J(S) =', sp.factor(JS.det()))
print('rank-3 minor =', minor)
print('R =', sp.factor(R))
print('lambda =', lam)
print('normalized singular gradient:')
for z, g in zip((z1,z2,z3,z4), grad):
    print(f'  dF/d{z} = {g}')
