# locate.py -- helpers (Sage/GAP) to locate a given (G,Phi) inside the stored enumeration enum_g{g}.json
import json, os
from sage.all import libgap
from cmenum import W_group

HERE = os.path.dirname(os.path.abspath(__file__))


class Stored:
    def __init__(self, g):
        self.g = g
        self.D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
        self.W, self.rho = W_group(g)
        self.groups = {G['gid']: G for G in self.D['groups']}
        self.cases = {}
        for c in self.D['cases']:
            self.cases.setdefault(c['gid'], []).append(c)
        self._gap = {}

    def gapgroup(self, gid):
        if gid not in self._gap:
            gens = [libgap.PermList([x + 1 for x in p]) for p in self.groups[gid]['gens']]
            self._gap[gid] = libgap.Group(gens)
        return self._gap[gid]

    def locate(self, G, Phi):
        """G: GAP group <= W containing rho; Phi: list of 0-based points. Returns (gid, stored case, conj elt)."""
        n = 2 * self.g
        order = int(libgap.Size(G))
        orbs = sorted(len(o) for o in libgap.Orbits(G, list(range(1, n + 1))))
        hits = []
        for gid, Gd in self.groups.items():
            if Gd['order'] != order or sorted(len(o) for o in Gd['orbits']) != orbs:
                continue
            H = self.gapgroup(gid)
            x = libgap.RepresentativeAction(self.W, G, H)
            if x == libgap.eval('fail'):
                continue
            img = libgap.OnSets(libgap.Set([p + 1 for p in Phi]), x)
            orb = libgap.Orbit(H, img, libgap.OnSets)
            mn = min(tuple(sorted(int(y) - 1 for y in s)) for s in orb)
            case = [c for c in self.cases[gid] if tuple(c['Phi']) == mn]
            assert len(case) == 1
            hits.append((gid, case[0]))
        assert len(hits) <= 1, "group conjugate to several stored reps?!"
        return hits[0] if hits else (None, None)
