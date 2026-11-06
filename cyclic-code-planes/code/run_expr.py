"""Run one of the JavaScript certificate programs of this folder with Node (BigInt arithmetic, no packages).
Each program is a single expression that evaluates to its verdict; it throws if a check fails.
python run_expr.py FILE.js"""
import sys, subprocess, os, tempfile
code = open(sys.argv[1], encoding='utf-8').read().strip()
rep = "(k,v)=>typeof v==='bigint'?v.toString():v"
fd, path = tempfile.mkstemp(suffix='.js')
with os.fdopen(fd, 'w', encoding='utf-8') as f:
    f.write("const __r=" + code + ";\nconsole.log(JSON.stringify(__r," + rep + ",1));\n")
try:
    r = subprocess.run(['node', '--stack-size=65500', path], capture_output=True, text=True)
    sys.stdout.write(r.stdout); sys.stderr.write(r.stderr)
    sys.exit(r.returncode)
finally:
    os.remove(path)
