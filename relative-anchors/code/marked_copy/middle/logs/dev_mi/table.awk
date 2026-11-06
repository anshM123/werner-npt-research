# summarises a prove_gen.py log: one line per family (vectors, split blocks, pivots num/den, time)
/^\[.*\] parts/ { fam=$0; sub(/ parts.*/, "", fam); piv=""; nv=""; sp="" }
/form computed/ { match($0, /[0-9]+ vectors/); nv=substr($0, RSTART, RLENGTH) }
/congruence done/ { s=$0; sub(/.*split blocks /, "", s); sub(/; parts.*/, "", s); sp=s }
/part[23] pivot / { n=$0; sub(/.*pivot /, "", n); nm=n; sub(/:.*/, "", nm); a=n; sub(/.*num /, "", a); sub(/ terms.*/, "", a); b=n; sub(/.*den /, "", b); sub(/ terms.*/, "", b); piv = piv (piv=="" ? "" : ", ") nm " " a "/" b }
/^\[.*\] (PD|PART|empty|PIVOT|TIME)/ { r=$0; t=r; sub(/.*\(/, "", t); sub(/s\).*/, "", t); res=r; sub(/^\[[^]]*\] /, "", res); sub(/  \(.*/, "", res); print "- " fam ": " nv ", split blocks " sp "; pivots " piv "; " t " s (" res ")" }
