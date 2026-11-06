#!/usr/bin/env bash
cd "$(dirname "$0")/../.."
bash guard.sh logs/dev/pe_small.log 600 prove_gen.py logs/dev/pe_small.pkl -1,n,n 0,n,n 1,n,n --imin=6; echo "exit $?" >> logs/dev/pe_small.log
bash guard.sh logs/dev/pe_smallt.log 600 prove_gen.py logs/dev/pe_smallt.pkl n,-2,n n,-1,n n,0,n --imin=5; echo "exit $?" >> logs/dev/pe_smallt.log
bash guard.sh logs/dev/pe_q2.log 900 prove_gen.py logs/dev/pe_q2.pkl 2,n,n --imin=6; echo "exit $?" >> logs/dev/pe_q2.log
bash guard.sh logs/dev/pe_t1.log 900 prove_gen.py logs/dev/pe_t1.pkl n,1,n --imin=5; echo "exit $?" >> logs/dev/pe_t1.log
echo done > logs/dev/pe_done.txt
