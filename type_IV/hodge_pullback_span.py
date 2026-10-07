# hodge_pullback_span.py -- produces hodge_pullback_span.out: numerical rank of the span of the pullbacks
# f^* Omega (f in Hom(A^n, A) (x) Q = L^n) of the Weil class, for n = 1..4, using part1_pullback_span of
# hodge_ring_gen.py with enough random trials (the default of 60 trials caps the rank at 60, which is why
# hodge_ring_gen_n3.out reports 60 of 100 for n = 3). Each output line: n (rank, (C(n+2,3))^2).
# First run as an inline command; this file contains the same code.
import hodge_ring_gen as h

for n, t in [(1, 20), (2, 60), (3, 250), (4, 600)]:
    print(n, h.part1_pullback_span(n, trials=t))
