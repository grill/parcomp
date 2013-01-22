#!/usr/bin/python
import sys
import os

usage = \
"usage: jupiterrun.py [number of threads] [first node] [programm name]"

if(len(sys.argv) < 4):
    print(usage)
    sys.exit(1)

threads = int(sys.argv[1])
node = int(sys.argv[2])
prog = sys.argv[3:]

#print(threads)
#print(node)
#print(prog)

command = ["mpirun"]

while(threads > 0):
    if(threads >= 8):
        command += ["-host", "jupiter%d" % node, "-np", "8"]
        node += 1
        threads -= 8
    else:
        command += ["-host", "jupiter%d" % node, "-np" , str(threads)]
        threads = 0

command += prog

#print(command)

os.execvp("mpirun", command)

