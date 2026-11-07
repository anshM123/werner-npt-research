"""Tables of NOTES.md, Section 8, from the logs and the meta data of the runs (not part of any proof).
python tables.py sd LOG [LOG ...]       (tsdet.py D logs: family, vectors, parts, B, terms of N~, degrees, time)
python tables.py d2 DIR [DIR ...]       (tdet2.py runs: family, Q2g size, bounds, L_v, grid, terms of N~, time)
python tables.py d3 DIR [DIR ...]       (tdet3.py runs)"""
import sys, os, re, pickle, glob


def sd_row(fn):
    s = open(fn).read()
    fam = re.search(r'^(M[TU] sym \[[^\]]*\])$', s, re.M)
    nv = re.search(r'form: (\d+) vectors', s)
    parts = re.search(r'parts (\d+) \+ (\d+)', s)
    q2 = re.search(r'Q2g (\d+) x', s)
    B = re.search(r'deg_i T <= B = (\d+)', s)
    nt = re.search(r'N~ with nu = .*?: (\d+) terms, (\d+) negative', s)
    dg = re.search(r'N~\(z, nu, s\): (\d+) terms, degrees \(([^)]*)\)', s)
    res = re.search(r'\] (PD on complement.*|FAILED.*) \(([\d.]+)s\)$', s, re.M)
    lab = fam.group(1).replace('MU sym ', '').replace('MT sym ', '') if fam else fn
    return (f'| {lab} | {nv.group(1) if nv else "?"} | {parts.group(1) + " + " + parts.group(2) if parts else "?"} | '
            f'{q2.group(1) if q2 else "?"} | {B.group(1) if B else "?"} | {dg.group(1) if dg else "?"} | '
            f'({dg.group(2) if dg else "?"}) | {nt.group(1) if nt else "?"} | {nt.group(2) if nt else "?"} | '
            f'{res.group(2) if res else "?"} | {"PD" if res and res.group(1).startswith("PD") else "not proved"} |')


def run_time(dirpath, prefix):
    """the sum over the worker logs of their reported total time for this family"""
    tag = os.path.basename(dirpath.rstrip('/\\'))
    tot = 0
    for fn in glob.glob(os.path.join('logs', f'{prefix}_w*_{tag}.log')):
        for line in open(fn):
            m = re.search(r'worker \d+(/\d+)?: (\d+) slices, (\d+)s', line)
            if m:
                tot += int(m.group(3))
    return tot


def final_info(prefix, tag):
    # a final step computed before the run script reran it (the rerun only reports the marker final_ok) is kept in
    # logs/PREFIX_final_TAG.computed.log
    fn = os.path.join('logs', f'{prefix}_final_{tag}.computed.log')
    if not os.path.exists(fn):
        fn = os.path.join('logs', f'{prefix}_final_{tag}.log')
    if not os.path.exists(fn):
        return None
    s = open(fn).read()
    m = re.search(r'(\d+) positive and (\d+) negative coefficients, degrees \(y, z, w, s\) <= \(([^)]*)\)', s) or \
        re.search(r'(\d+) positive and (\d+) negative coefficients, degrees <= \(([^)]*)\)', s)
    ok = s.strip().endswith('PASS')
    return (m.group(1), m.group(2), m.group(3)) if m else None, ok


def d_rows(kind, dirs):
    for d in dirs:
        with open(os.path.join(d, 'meta.pkl'), 'rb') as fh:
            meta = pickle.load(fh)
        tag = os.path.basename(d.rstrip('/\\'))
        prefix = 'd2' if kind == 'd2' else 'd3'
        fi = final_info(prefix, tag)
        nsl = len(glob.glob(os.path.join(d, 's_*.pkl')))
        if kind == 'd2':
            grid = (meta['Bv'] + 1) * (meta['Bi'] + 1)
            b = f"{meta['v']}: {meta['Bv']}, i: {meta['Bi']}"
            lv = f"{meta['v']}: {meta['Lv']}"
        else:
            grid = (meta['B1'] + 1) * (meta['B2'] + 1)
            b = f"{meta['v1']}: {meta['B1']}, {meta['v2']}: {meta['B2']} (free {meta['u']}: {meta['bounds'][meta['u']]})"
            lv = f"{meta['v1']}: {meta['L1']}, {meta['v2']}: {meta['L2']}"
        if fi and fi[0]:
            terms, neg, deg = fi[0]
            st = 'PD' if fi[1] else 'not proved'
        else:
            terms = neg = deg = '-'
            st = 'slices %d of %d' % (nsl, grid)
        q2 = len(meta['keep2'])
        print(f"| {meta['spec']} | {q2} + {len(meta['three'])} | {b} | {lv} | {grid} | {terms} | {neg} | ({deg}) | "
              f"{run_time(d, prefix)} | {st} |")


if __name__ == '__main__':
    kind = sys.argv[1]
    if kind == 'sd':
        for fn in sys.argv[2:]:
            print(sd_row(fn))
    else:
        d_rows(kind, sys.argv[2:])
