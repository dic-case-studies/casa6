import os
import sys
print("you cannot run the casatasks module",file=sys.stderr)
sys.stderr.flush()
os._exit(1)
