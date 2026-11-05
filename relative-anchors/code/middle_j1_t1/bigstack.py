"""Run a script in a thread with a 128 MB stack.  Two runs of bareiss_t1.py 5 ended with a segmentation fault inside
the polynomial library; a stack overflow on the default 1 MB Windows stack is the likely cause, and the run through
this wrapper passed that point.  Usage: python bigstack.py script.py args..."""
import sys, threading, runpy

script = sys.argv[1]
sys.argv = sys.argv[1:]
threading.stack_size(128 * 1024 * 1024)
th = threading.Thread(target=lambda: runpy.run_path(script, run_name='__main__'))
th.start()
th.join()
