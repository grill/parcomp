import re
import os
import sys
import math

path = sys.argv[1]
pathout = sys.argv[2]
dirs = os.listdir(path)

def isFloat(string):
     try:
         float(string)
         return True
     except ValueError:
         return False

for filename in dirs:
   if "time" in filename or "reduce" in filename:
      inf = open(path + "/" + filename)
      outf = open(pathout + "/" + filename, "w")
      line = inf.readline()
      alle = line.split(';')
      first = float(alle[1])
      for x in alle[1:len(alle)-1]:
         outf.write(str(float(x)/first) + "\n")
   if "work" in filename:
      inf = open(path + "/" + filename)
      outf = open(pathout + "/" + filename, "w")
      outf2 = open(pathout + "/" + filename + "_time", "w")
      lines = inf.readlines()
      for line in lines[:len(lines)-1]:
         alld = line.split(';')
         floats = map(float, filter(isFloat, alld[1:]))
         outf.write(str(math.fsum(floats)) + "\n")
         if(line != lines[0]):
            outf2.write(alld[0] + "\n")
