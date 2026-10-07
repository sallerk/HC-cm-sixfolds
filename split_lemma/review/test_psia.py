# test_psia.py -- checks of Lemma lem:psi-a (note.tex, Section 3) with own code.
#  (i)  psi_a(x,y) = Tr_{E/Q}(sqrt(-d) a x ybar) is alternating, E-compatible; it is a Riemann form
#       (psi(Jx,Jy)=psi(x,y), psi(x,Jx)>0, Milne's convention) iff sign v(a) = eps(v) for all v,
#       where J is the complex structure defined by the CM type Phi with signs eps.  (numerical)
#  (ii) Milne's phi(x,y) = psi(x,beta y) + beta psi(x,y) equals 2d * conj(h_a(x,y))  (exact),
#       det(h_a on the power basis) = N(a) * disc(power basis)  (exact), signature of h_a.
import random, sys, time
from fractions import Fraction
import numpy as np
import mpmath
from alg import Alg, det_frac
from inst import random_instance, random_pattern, SHAPES

rng = random.Random(20261006)
mpmath.mp.dps = 50


def E_mul(Fa, x, y, d):
    # x = (xr, xi) meaning xr + sqrt(-d) xi, with xr, xi in F
    xr, xi = x
    yr, yi = y
    r = [[p - d * q for p, q in zip(A, B)] for A, B in zip(Fa.mul(xr, yr), Fa.mul(xi, yi))]
    i = [[p + q for p, q in zip(A, B)] for A, B in zip(Fa.mul(xr, yi), Fa.mul(xi, yr))]
    return (r, i)


def conjE(x):
    return (x[0], [[-c for c in comp] for comp in x[1]])


def TrF(Fa, z):
    return sum(F.trace(zi) for F, zi in zip(Fa.factors, z))


def psi(Fa, a, x, y, d):
    # Tr_{E/Q}(sqrt(-d) a x ybar);  Tr_{E/Q}(u + sqrt(-d) v) = 2 Tr_{F/Q}(u)
    z = E_mul(Fa, x, conjE(y), d)
    az = (Fa.mul(a, z[0]), Fa.mul(a, z[1]))
    w = ([[Fraction(0)] * F.m for F in Fa.factors], Fa.one())
    t = E_mul(Fa, w, az, d)
    return 2 * TrF(Fa, t[0])


def h(Fa, a, x, y, d):
    # Tr_{E/K}(a x ybar) = Tr_{F/Q}(real part) + sqrt(-d) Tr_{F/Q}(imag part); returned as (re, im)
    z = E_mul(Fa, x, conjE(y), d)
    return (TrF(Fa, Fa.mul(a, z[0])), TrF(Fa, Fa.mul(a, z[1])))


def rand_F(Fa, B=4):
    return [[Fraction(rng.randint(-B, B)) for _ in range(F.m)] for F in Fa.factors]


def rand_E(Fa):
    return (rand_F(Fa), rand_F(Fa))


def basisE(Fa):
    """Q-basis of E: (f_k, 0) then (0, f_k), f_k running over power bases of the factors."""
    out = []
    for part in (0, 1):
        for i, F in enumerate(Fa.factors):
            for k in range(F.m):
                comp = [[Fraction(0)] * G.m for G in Fa.factors]
                comp[i][k] = Fraction(1)
                zero = [[Fraction(0)] * G.m for G in Fa.factors]
                out.append((comp, zero) if part == 0 else (zero, comp))
    return out


stats = dict(inst=0, phi_ok=0, det_ok=0, riem_ok=0, riem_pos=0, riem_neg=0, sig_ok=0, alt_ok=0)
bad = []
t0 = time.time()
shapes = ['field2', 'field4', 'field6', '1+1', '2+2', '1+3', '3+3', '2+4', '2+2+2', '1+1+4']
for it in range(120):
    shape = shapes[it % len(shapes)]
    polys, d = random_instance(shape, rng)
    Fa = Alg(polys)
    n = Fa.n
    N = 2 * n
    stats['inst'] += 1
    # random a, random pattern; half the time force a to have the pattern's signs
    eps = random_pattern(N, rng.randint(0, N), rng)
    for trial in range(4):
        a = rand_F(Fa, 6)
        if any(all(c == 0 for c in comp) for comp in a):
            continue
        try:
            sg = Fa.signs(a)
        except ValueError:
            continue
        # (ii) exact: phi = 2d conj(h)
        for _ in range(3):
            x, y = rand_E(Fa), rand_E(Fa)
            dd = Fraction(d)
            # beta y = sqrt(-d) * y = (-d yi, yr)
            by = ([[-d * c for c in comp] for comp in y[1]], y[0])
            p1 = psi(Fa, a, x, by, d)
            p2 = psi(Fa, a, x, y, d)
            phi = (p1, p2)                      # p1 + sqrt(-d) p2
            hr, hi = h(Fa, a, x, y, d)
            if phi == (2 * d * hr, -2 * d * hi):
                stats['phi_ok'] += 1
            else:
                bad.append(('phi', shape, polys, d, a))
            if psi(Fa, a, x, x, d) == 0 and psi(Fa, a, x, y, d) == -psi(Fa, a, y, x, d):
                stats['alt_ok'] += 1
            else:
                bad.append(('alt', shape))
        # (ii) exact: det of h_a on the power basis
        G = Fa.gram_ha(a)
        if det_frac(G) == Fa.norm(a) * Fa.disc():
            stats['det_ok'] += 1
        else:
            bad.append(('det', shape, polys, d, a))
        # signature of h_a (real symmetric G) = (#pos signs, #neg signs) of a
        ev = np.linalg.eigvalsh(np.array([[float(c) for c in row] for row in G]))
        if (sum(ev > 0), sum(ev < 0)) == (sg.count(1), sg.count(-1)):
            stats['sig_ok'] += 1
        else:
            bad.append(('sig', shape))
        # (i) numerical Riemann-form test with J from the CM type with signs eps
        B = basisE(Fa)
        Psi = np.array([[float(psi(Fa, a, u, v, d)) for v in B] for u in B])
        # Theta: V_R -> C^N, x -> (phi_v(x))_v ; phi_v(xr + w xi) = v(xr) + eps_v i sqrt(d) v(xi)
        places = Fa.places()
        Th = np.zeros((2 * N, 2 * N))
        for col, (xr, xi) in enumerate(B):
            for pi, (fi, ri) in enumerate(places):
                F = Fa.factors[fi]
                r = F.roots[ri]
                vr = float(F.evalf(xr[fi], r))
                vi = float(F.evalf(xi[fi], r))
                Th[2 * pi, col] = vr
                Th[2 * pi + 1, col] = eps[pi] * np.sqrt(d) * vi
        I2 = np.zeros((2 * N, 2 * N))
        for pi in range(N):
            I2[2 * pi, 2 * pi + 1] = -1.0   # multiplication by i on (re, im)
            I2[2 * pi + 1, 2 * pi] = 1.0
        J = np.linalg.solve(Th, I2 @ Th)
        S = Psi @ J                     # S[u,v] = psi(u, J v)
        inv_ok = np.allclose(J.T @ Psi @ J, Psi, atol=1e-6 * np.abs(Psi).max())
        sym_ok = np.allclose(S, S.T, atol=1e-6 * np.abs(S).max())
        evS = np.linalg.eigvalsh((S + S.T) / 2)
        posdef = bool(evS.min() > 0)
        predicted = all(s == e for s, e in zip(sg, eps))
        if inv_ok and sym_ok and posdef == predicted:
            stats['riem_ok'] += 1
            stats['riem_pos' if predicted else 'riem_neg'] += 1
        else:
            bad.append(('riemann', shape, polys, d, a, eps, inv_ok, sym_ok, posdef, predicted))
        # also test a with exactly the pattern's signs (positive case) when cheap
        if trial == 0 and N <= 6:
            for _ in range(200):
                a2 = rand_F(Fa, 6)
                try:
                    if Fa.signs(a2) == eps:
                        break
                except ValueError:
                    pass
            else:
                continue
            Psi2 = np.array([[float(psi(Fa, a2, u, v, d)) for v in B] for u in B])
            S2 = Psi2 @ J
            ev2 = np.linalg.eigvalsh((S2 + S2.T) / 2)
            if ev2.min() > 0 and np.allclose(S2, S2.T, atol=1e-6 * np.abs(S2).max()):
                stats['riem_ok'] += 1
                stats['riem_pos'] += 1
            else:
                bad.append(('riemann+', shape, polys, d, a2, eps))
print(stats)
print("failures:", len(bad))
for b in bad[:10]:
    print(b)
print(f"time {time.time()-t0:.1f}s")
