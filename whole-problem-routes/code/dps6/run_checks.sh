#!/bin/sh
# Exact checks of the six-copy certificate (Theorem 7.32) and of the four-copy test certificate.
# Needs Python with numpy and python-flint; about four minutes in all.
cd "$(dirname "$0")"
python verify_l2.py certs/cert_n6_a3_4_t1_16384.json.xz | tee logs/verify_n6.log
python cert_l2.py 6 3 4 check certs/cert_n6_a3_4_t1_16384.json.xz | tee logs/check_n6.log
python verify_l2.py certs/cert_n4_a3_4_t1_512.json.xz | tee logs/verify_n4.log
python cert_l2.py 4 3 4 check certs/cert_n4_a3_4_t1_512.json.xz | tee logs/check_n4.log
