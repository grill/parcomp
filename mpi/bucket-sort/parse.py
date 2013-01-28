import re
import os
import sys
import math

path = sys.argv[1]
pathout = sys.argv[2]
dirs = os.listdir(path)
mis = open("missing.txt", "a")
d = {}

def isFloat(string):
     try:
         float(string)
         return True
     except ValueError:
         return False

for filename in dirs:
    inf = open(path + "/" + filename)
    line = inf.readline()
    if (not line or not isFloat(line)):
       mis.write(filename + "\n")
       line = "1"
    pars = filename.split("_")
    if (not d.has_key(pars[1] + "_" + pars[2])):
       d[pars[1] + "_" + pars[2]] = float(line)
    outf = open(pathout + "/" + pars[1] + "_" + pars[2], "a")
    outf.write(str(float(line)/float(d[pars[1] + "_" + pars[2]])) + "\n")