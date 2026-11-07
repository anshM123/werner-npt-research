#!/usr/bin/env bash
# Exact verification of every certificate in certs/, and negative controls (run from this folder; about 10 minutes,
# under 2.5 GB).
# Each run goes through guard25.sh, which stops a run above 2.5 GB or past its time limit.
set -u
for c in certs/cert_n3_a3_4_t1_64_exact.txt certs/cert_n3_a3_4_t1_64_exactB.txt certs/cert_n3_a3_4_t1953_125000.txt \
         certs/cert_n3_a3_4_t1_128.txt certs/cert_n3_a81_125_t0.txt certs/cert_n4_a3_4_t1_256_exact.txt.xz \
         certs/cert_n4_a3_4_t1_512.txt.xz certs/cert_n5_a3_4_t1_1024_exact.txt.xz certs/cert_n5_a3_4_t1_2048.txt.xz; do
  b=$(basename "$c"); b=${b%.xz}; b=${b%.txt}
  bash guard25.sh "logs/independent_${b#cert_}.log" 3600 verify_independent.py "$c" 27
  echo "$b: $(tail -n 2 "logs/independent_${b#cert_}.log" | head -n 1)"
done
for c in certs/cert_n3_a3_4_t1_64_exact.txt certs/cert_n3_a3_4_t1_64_exactB.txt; do
  b=$(basename "$c" .txt)
  bash guard25.sh "logs/crosscheck3_${b#cert_}.log" 1800 crosscheck3.py "$c" 6
  echo "$b: $(tail -n 1 "logs/crosscheck3_${b#cert_}.log")"
done
python unpack_certs.py > /dev/null
bash guard25.sh logs/verify_n3_a3_4_t1_64_exact.log 600 exact_verify.py certs/cert_n3_a3_4_t1_64_exact.txt
echo "exact_verify n3: $(tail -n 1 logs/verify_n3_a3_4_t1_64_exact.log)"
bash guard25.sh logs/modverify_n5_a3_4_t1_1024_exact.log 1800 exact_verify_mod.py certs/cert_n5_a3_4_t1_1024_exact.txt
echo "exact_verify_mod n5: $(tail -n 1 logs/modverify_n5_a3_4_t1_1024_exact.log)"
# negative controls: the exact certificates with t raised just above the product value 4^-n must fail
for spec in "certs/cert_n3_a3_4_t1_64_exact.txt t=1/63" "certs/cert_n4_a3_4_t1_256_exact.txt.xz t=1/255"             "certs/cert_n5_a3_4_t1_1024_exact.txt.xz t=1/1023"; do
  set -- $spec
  b=$(basename "$1"); b=${b%.xz}; b=${b%.txt}; tt=${2#t=}; tt=${tt/\//_}
  bash guard25.sh "logs/independent_negative_${b#cert_}_t${tt}.log" 3600 verify_independent.py "$1" 27 "$2"
  echo "negative control $b with $2: $(tail -n 2 "logs/independent_negative_${b#cert_}_t${tt}.log" | head -n 1)"
done
rm -f logs/*.pid
