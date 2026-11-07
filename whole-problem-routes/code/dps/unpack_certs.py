"""Decompress certs/*.txt.xz to certs/*.txt (the scripts exact_verify*.py read plain text; verify_independent.py
reads either).  usage: python unpack_certs.py"""
import glob, lzma, os

for p in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'certs', '*.txt.xz'))):
    out = p[:-3]
    if not os.path.exists(out):
        with lzma.open(p, 'rb') as f, open(out, 'wb') as g:
            g.write(f.read())
    print(out)
