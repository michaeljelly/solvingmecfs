"""Parse Metabolomics Workbench mwTab text files and separate results files."""
import re, io
import numpy as np, pandas as pd

RAW = '/mnt/project-files/mecfs-research/downloads/mwb/raw/'


def ssf(path):
    """SUBJECT_SAMPLE_FACTORS -> DataFrame (sample, subject, factors..., additional data...)."""
    rows = []
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('SUBJECT_SAMPLE_FACTORS'):
            p = line.rstrip('\n').split('\t')
            p = p[1:] + [''] * 4
            subj, samp, fac, add = p[0].strip(), p[1].strip(), p[2], p[3]
            r = dict(sample=samp, subject=subj)
            for kv in fac.split('|'):
                if ':' in kv:
                    k, v = kv.split(':', 1); r[k.strip()] = v.strip()
            for kv in re.split(r';\s*', add):
                if '=' in kv:
                    k, v = kv.split('=', 1); r['x_' + k.strip()] = v.strip()
            rows.append(r)
    return pd.DataFrame(rows)


def block(path, start, end):
    lines, on = [], False
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith(start): on = True; continue
        if line.startswith(end): break
        if on: lines.append(line.rstrip('\n'))
    return lines


def data(path):
    """MS_METABOLITE_DATA -> (X samples x metabolites, float), plus METABOLITES annotation table."""
    L = block(path, 'MS_METABOLITE_DATA_START', 'MS_METABOLITE_DATA_END')
    if not L: L = block(path, 'NMR_BINNED_DATA_START', 'NMR_BINNED_DATA_END')
    hdr = L[0].split('\t')
    body = [l.split('\t') for l in L[1:] if not l.startswith('Factors')]
    n = len(hdr)
    body = [(b + [''] * n)[:n] for b in body]
    D = pd.DataFrame(body, columns=['name'] + [h.strip() for h in hdr[1:]])
    D = D[D.name.str.strip() != '']
    D = D.set_index('name')
    D = D.loc[:, [c for c in D.columns if c]]
    X = D.apply(lambda s: pd.to_numeric(s.str.replace(',', ''), errors='coerce')).T
    A = block(path, 'METABOLITES_START', 'METABOLITES_END')
    ann = None
    if A:
        h = A[0].split('\t'); rr = [(a.split('\t') + [''] * len(h))[:len(h)] for a in A[1:]]
        ann = pd.DataFrame(rr, columns=h)
    return X, ann
