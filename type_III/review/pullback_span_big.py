"""pullback_span_big.py: k = 4 (all five test cases) and k = 5 (two cases) of pullback_span.run_case."""
import json, sys, time
from pullback_span import run_case
cases = [(-1, -1, [0, 1, 0, 0]), (-1, -1, [0, 1, 1, 1]), (-2, -5, [0, 0, 1, 0]),
         (-3, -7, [0, 1, 2, 0]), (-11, -13, [0, 1, 1, 0])]
res = []
for (a, c, q) in cases:
    t0 = time.time(); r = run_case(a, c, q, 4, seed=4 * 17 + a); r["secs"] = round(time.time() - t0, 1)
    print(json.dumps(r)); sys.stdout.flush(); res.append(r)
for (a, c, q) in cases[1:3]:
    t0 = time.time(); r = run_case(a, c, q, 5, seed=5 * 17 + a); r["secs"] = round(time.time() - t0, 1)
    print(json.dumps(r)); sys.stdout.flush(); res.append(r)
json.dump(res, open("pullback_span_big_out.json", "w"), indent=0)
