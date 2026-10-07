\\ produce_lib.gp -- PARI/GP library for the main computation (run through cypari by produce.py).
\\ E = F.K with F totally real sextic (polF) and K imaginary quadratic (polK), absolute polynomial R(x),
\\ theta = a + k*b (polcompositum flag 1), a = root of polF, b = root of polK.
\\ Complex conjugation tau: a -> a, b -> tr(b) - b.

buildE(polF, polK) =
{
  my(C = polcompositum(polF, polK, 1));
  if (#C != 1, error("compositum not a field"));
  C = C[1];
  my(R = C[1], a = C[2], b = C[3], k = C[4]);
  my(trb = -polcoeff(polK, 1) / polcoeff(polK, 2));
  my(tauimg = lift(a + k * (trb - b)));
  [R, lift(a), lift(b), k, tauimg];
}

\\ apply tau to an nf element given as polynomial/polmod in x mod R
tauapply(R, tauimg, z) = Mod(subst(lift(Mod(z, R)), x, tauimg), R);

\\ identify the 12 degree-1 primes of E above p with pairs (j, eta):
\\ j = index of the root of polF mod p (sorted), eta = index (0/1) of the root of polK mod p (sorted)
primepairs(nf, polF, polK, a, b, p) =
{
  my(dec = idealprimedec(nf, p), rF, rK, tab = matrix(6, 2), cnt = 0);
  if (#dec != 12, return(0));
  rF = vecsort(apply(z -> lift(z), polrootsmod(polF, p)));
  rK = vecsort(apply(z -> lift(z), polrootsmod(polK, p)));
  if (#rF != 6 || #rK != 2, return(0));
  for (i = 1, 12,
    my(P = dec[i], jj = 0, ee = 0);
    if (P.f != 1, return(0));
    for (j = 1, 6, if (nfeltval(nf, Mod(a, nf.pol) - rF[j], P) > 0, if (jj, error("double F"), jj = j)));
    for (e = 1, 2, if (nfeltval(nf, Mod(b, nf.pol) - rK[e], P) > 0, if (ee, error("double K"), ee = e)));
    if (!jj || !ee, error("unidentified prime"));
    tab[jj, ee] = i; cnt++);
  [dec, tab, rF, rK];
}

\\ order of an ideal class: bnfisprincipal flag 0 gives the class vector
classorder(bnf, I) =
{
  my(e = bnfisprincipal(bnf, I, 0), cyc = bnf.cyc, h = 1);
  for (i = 1, #cyc, h = lcm(h, cyc[i] / gcd(e[i], cyc[i])));
  h;
}

\\ find unit v of E with v * tau(v) = u, or 0
normunit(bnf, R, tauimg, u) =
{
  my(fu = bnf.fu, r = #fu, t = varhigher("t"));
  forvec(e = vector(r, i, [0, 1]),
    my(eps = Mod(1, R));
    for (i = 1, r, if (e[i], eps *= Mod(lift(fu[i]), R)));
    my(w = Mod(lift(u), R) / (eps * tauapply(R, tauimg, eps)));
    my(rts = nfroots(bnf, t^2 - lift(w)));
    for (i = 1, #rts,
      my(y0 = Mod(lift(rts[i]), R));
      if (tauapply(R, tauimg, y0) == y0, return(eps * y0))));
  0;
}

\\ Weil number for the CM type given by the choice vector S01 (length 6, entries 1/2 = eta index).
\\ returns [pi (polmod), q, h, usedsquare]
weilnumber(bnf, R, tauimg, dec, tab, p, S01, qmax = 0) =
{
  my(nf = bnf.nf, I = idealhnf(nf, 1), h, Ih, al, u, v, piw, q, sq = 0);
  for (j = 1, 6, I = idealmul(nf, I, dec[tab[j, S01[j]]]));
  h = classorder(bnf, I);
  if (qmax && p^h > qmax, return([0, p^h, h, 0]));
  Ih = idealpow(nf, I, h);
  al = bnfisprincipal(bnf, Ih, 5);
  if (al[1] != 0 * al[1], error("not principal"));
  al = al[2]; if (type(al) == "t_MAT", al = nffactorback(nf, al));
  al = Mod(nfbasistoalg(nf, al), R);
  if (idealhnf(nf, lift(al)) != Ih, error("generator check failed"));
  u = al * tauapply(R, tauimg, al) / p^h;
  if (abs(norm(u)) != 1 || denominator(content(lift(u))) != 1 && 0, error("u not unit"));
  v = normunit(bnf, R, tauimg, u);
  if (v != 0, piw = al / v; q = p^h, piw = al^2 / u; q = p^(2*h); sq = 1);
  if (piw * tauapply(R, tauimg, piw) != q, error("pi*taupi != q"));
  [piw, q, h, sq];
}

\\ lcm of all n with eulerphi(n) | m  (multiple of #roots of unity of any field of degree m)
rootsofunitybound(m) =
{
  my(N = 1);
  for (n = 1, 2 * m^2 + 10, if (m % eulerphi(n) == 0, N = lcm(N, n)));
  N;
}

\\ squarefreeness of charpoly(pi^N) mod ell, computed from the Weil polynomial P
geomsimple_modl(P, N, ell) =
{
  my(Pl = P * Mod(1, ell), z = Mod(x, Pl)^N, Q = charpoly(z));
  poldegree(gcd(Q, deriv(Q))) == 0;
}
