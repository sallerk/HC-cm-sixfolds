# fetch_apfv_data.py -- download the two data files of Arango-Pineros, Frengley, Vemulapalli [APFV] that are read by
# controls.py and verify_controls.py, and check them against the SHA-256 of the copies used for the stored outputs.
# Source: public GitHub repository SamFrengley/exceptional-tate-classes, folder 1-4-theorem-data.
# The files are saved in apfv_data/ next to this script (apfv_data/ is not part of this repository).
#   usage:  python fetch_apfv_data.py          (pure Python 3, standard library only)
import hashlib, json, os, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = 'SamFrengley/exceptional-tate-classes'
FOLDER = '1-4-theorem-data'
FILES = [
    ('4-0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00.m',
     'b677bff503e8bd79c11046a2e6b7b4aa1ace2c5b8573a5aee9286b21089ba059'),
    ('6-0.000_0.000_0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00_1.00_1.00.m',
     '29b40025944cb7bd7f3f95ce4af5d8ee83c804f59017f070ce6f90275f909de6'),
]


def default_branch():
    try:
        with urllib.request.urlopen('https://api.github.com/repos/%s' % REPO, timeout=30) as r:
            return json.load(r)['default_branch']
    except Exception as e:
        print('GitHub API not reachable (%s); trying branch "main"' % e)
        return 'main'


def main():
    branch = default_branch()
    print('repository %s, branch %s, folder %s' % (REPO, branch, FOLDER))
    outdir = os.path.join(HERE, 'apfv_data')
    os.makedirs(outdir, exist_ok=True)
    allok = True
    for name, expected in FILES:
        url = 'https://raw.githubusercontent.com/%s/%s/%s/%s' % (REPO, branch, FOLDER, name)
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read()
        with open(os.path.join(outdir, name), 'wb') as f:
            f.write(data)
        h = hashlib.sha256(data).hexdigest()
        ok = (h == expected)
        allok &= ok
        print('%s  %d bytes  sha256 %s  %s' % (name, len(data), h, 'OK' if ok else 'DIFFERENT (expected %s)' % expected))
    print('all files identical to the copies used for the stored outputs' if allok else
          'WARNING: the downloaded data differ from the copies used for the stored outputs')
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
