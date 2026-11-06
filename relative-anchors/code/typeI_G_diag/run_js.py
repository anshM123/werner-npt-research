"""Run one of the JavaScript programs of this folder with Node (BigInt arithmetic, no packages).
python run_js.py FILE.js [RUNNER]
A file without RUNNER is a function body that returns its verdict; with RUNNER the file defines functions and RUNNER()
is called.  The internal time budgets of the certificates (45 s per step on the machine they were written for) are raised to
ten hours; nothing else is changed.  The verdict is printed as JSON."""
import sys, re, subprocess, os, tempfile
src = sys.argv[1]
call = sys.argv[2] if len(sys.argv) > 2 else None
code = open(src, encoding='utf-8').read()
code = re.sub(r'rfDeadline=rfStart\+\d+', 'rfDeadline=rfStart+36000000', code)
code = re.sub(r'mpDeadline=mpStart\+\d+', 'mpDeadline=mpStart+36000000', code)
rep = "(k,v)=>typeof v==='bigint'?v.toString():v"
if call:
    body = code + f"\n;const __r={call}();\nconsole.log(JSON.stringify(__r,{rep},1));\n"
else:
    body = "const __r=(function(){\n" + code + f"\n}})();\nconsole.log(JSON.stringify(__r,{rep},1));\n"
fd, path = tempfile.mkstemp(suffix='.js')
with os.fdopen(fd, 'w', encoding='utf-8') as f:
    f.write(body)
try:
    r = subprocess.run(['node', '--stack-size=65500', path], capture_output=True, text=True)
    sys.stdout.write(r.stdout); sys.stderr.write(r.stderr)
    sys.exit(r.returncode)
finally:
    os.remove(path)
