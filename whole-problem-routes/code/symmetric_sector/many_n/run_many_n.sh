#!/usr/bin/env bash
# Exact verification of the certificates of Proposition 7.16 for n = 9..21 (coupled-basis verifier verify_exact.py),
# plus the consistency checks of the coupled basis (check_cpl.py) and of the word blocks against the monomial-basis code
# of the parent folder (check_words.py).  About ten minutes; run from this folder: bash run_many_n.sh > logs/run_many_n.out
cd "$(dirname "$0")"
export PYTHONPATH="$(pwd)/..:$(pwd)"
python check_cpl.py > logs/check_cpl.log 2>&1; echo "check_cpl: $(tail -n 1 logs/check_cpl.log)"
python check_words.py > logs/check_words.log 2>&1; echo "check_words: $(tail -n 1 logs/check_words.log)"
for n in 9 10 11 12 13 14 15 16 17 18 19 20 21; do
  case $n in 9) c=certs/cert_K1_n9.pkl;; 10) c=certs/restr_K1_n10.pkl;; 21) c=certs/cert_K2_n21.pkl;; *) c=certs/cert_K1_n$n.pkl;; esac
  python verify_exact.py $n $c 200 > logs/verify_exact_n$n.log 2>&1
  echo "n = $n: $(grep RESULT logs/verify_exact_n$n.log | sed 's/^ *//')"
done
