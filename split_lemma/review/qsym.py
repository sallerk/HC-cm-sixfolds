# qsym.py -- own implementation (review of Theorem A) of Hilbert symbols over Q_p and R,
# square classes and small helpers.  Pure Python integers / Fractions; no PARI.
# Formulas: Serre, A Course in Arithmetic, Ch. III, Thm. 1.
from fractions import Fraction


def vp(n, p):
    """p-adic valuation of a nonzero integer."""
    n = abs(int(n))
    if n == 0:
        raise ValueError("vp(0)")
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def to_int_class(x):
    """Return a nonzero integer in the same square class as the nonzero rational x."""
    x = Fraction(x)
    if x == 0:
        raise ValueError("zero")
    return x.numerator * x.denominator


def legendre(u, p):
    u %= p
    if u == 0:
        return 0
    return 1 if pow(u, (p - 1) // 2, p) == 1 else -1


def hilbert(a, b, p):
    """Hilbert symbol (a,b)_p for nonzero rationals a,b; p prime, or p == 0 for the real place."""
    a = to_int_class(a)
    b = to_int_class(b)
    if p == 0:
        return -1 if (a < 0 and b < 0) else 1
    al = vp(a, p)
    be = vp(b, p)
    u = a // p ** al
    w = b // p ** be
    if p != 2:
        s = 1
        if (al * be * ((p - 1) // 2)) % 2 == 1:
            s = -s
        if be % 2 == 1:
            s *= legendre(u, p)
        if al % 2 == 1:
            s *= legendre(w, p)
        return s
    # p == 2
    def eps(x):
        return ((x - 1) // 2) % 2

    def omg(x):
        return ((x * x - 1) // 8) % 2
    e = (eps(u) * eps(w) + al * omg(w) + be * omg(u)) % 2
    return -1 if e else 1


def factor_over(n, primes):
    """Factor the nonzero integer n completely over the given primes; return (dict, cofactor)."""
    n = abs(int(n))
    f = {}
    for p in primes:
        if n == 1:
            break
        k = 0
        while n % p == 0:
            n //= p
            k += 1
        if k:
            f[p] = k
    return f, n


def small_primes(B):
    s = bytearray([1]) * (B + 1)
    s[0] = s[1] = 0
    for i in range(2, int(B ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(B + 1) if s[i]]


def is_norm_from_K(x, d, primes_hint=None):
    """Decide whether the nonzero rational x is a norm from K = Q(sqrt(-d)) (d>0 squarefree),
    by the Hasse norm theorem: (x,-d)_p = 1 at all places.  Needs the primes dividing x and 2d;
    these are found with sympy.factorint unless primes_hint covers them."""
    xi = to_int_class(x)
    if xi < 0:
        return False, ('inf',)
    ps = set([2]) | set(prime_factors(d))
    if primes_hint is not None:
        f, cof = factor_over(xi, sorted(primes_hint))
        if cof != 1:
            from sympy import factorint
            ps |= set(factorint(cof).keys())
        ps |= set(f.keys())
    else:
        ps |= set(prime_factors(xi))
    bad = [p for p in sorted(ps) if hilbert(xi, -d, p) == -1]
    return (len(bad) == 0), tuple(bad)


def prime_factors(n):
    from sympy import factorint
    n = abs(int(n))
    if n <= 1:
        return []
    return sorted(factorint(n).keys())


def squarefree_part(n, primes=None):
    """Signed squarefree part of a nonzero integer."""
    n = int(n)
    sgn = -1 if n < 0 else 1
    n = abs(n)
    if primes is None:
        from sympy import factorint
        f = factorint(n)
    else:
        f, cof = factor_over(n, primes)
        if cof != 1:
            from sympy import factorint
            f.update(factorint(cof))
    r = 1
    for p, k in f.items():
        if k % 2:
            r *= p
    return sgn * r


# ---------------------------------------------------------------------------------------------
# brute-force (theory-free) local solvability, used to cross-check `hilbert`
def has_primitive_zero_mod(coeffs, p, k):
    """True iff sum c_i x_i^2 == 0 mod p^k has a solution with some x_i not divisible by p.
    coeffs: list of integers (3 or 4 entries).  Exhaustive, vectorised with numpy."""
    import numpy as np
    m = p ** k
    r = np.arange(m, dtype=np.int64)
    sq = (r * r) % m
    c = [int(ci) % m for ci in coeffs]
    nvar = len(c)
    # values c_i x_i^2 mod m for each variable, plus a mask of units (x not divisible by p)
    vals = [(ci * sq) % m for ci in c]
    unit = (r % p) != 0
    # primitive: at least one coordinate is a unit.  Count solutions with total sum == 0 mod m
    # and subtract those with all coordinates divisible by p.
    def count(vlist):
        # convolution of value distributions mod m
        dist = np.bincount(vlist[0], minlength=m).astype(np.int64)
        for v in vlist[1:]:
            h = np.bincount(v, minlength=m).astype(np.int64)
            new = np.zeros(m, dtype=np.int64)
            nz = np.nonzero(h)[0]
            for t in nz:
                new += np.roll(dist, t) * h[t]
            dist = new
        return dist[0]
    total = count(vals)
    nonunit_vals = [v[~unit] for v in vals]
    nonprim = count(nonunit_vals)
    return (total - nonprim) > 0


def hilbert_bruteforce(a, b, p, k=None):
    """(a,b)_p via exhaustive search for primitive zeros of z^2 - a x^2 - b y^2 mod p^k."""
    a = squarefree_part(to_int_class(a))
    b = squarefree_part(to_int_class(b))
    if k is None:
        k = 2 if p != 2 else 5
    ok = has_primitive_zero_mod([1, -a, -b], p, k)
    return 1 if ok else -1
