import re
import os
import sys
import math

path = sys.argv[1]
pathout = sys.argv[2]
dirs = os.listdir(path)
mis = open("missing.txt", "w")

def isFloat(string):
     try:
         float(string)
         return True
     except ValueError:
         return False

for filename in dirs:
   if "iterative" in filename:
      inf = open(path + "/" + filename)
      line = inf.readline()
      if(not line):
         mis.write(filename + "\n")
      else:
         pars = filename.split("_")
         outf = open(pathout + "/" + "taskpara_"+pars[1]+"_"+pars[2], "a")
         outf.write(line + "\n")
         #os.close(outf)
   if "prefixsum" in filename:
      inf = open(path + "/" + filename)
      line = inf.readline()
      if(not line):
         mis.write(filename + "\n")
      else:
         pars = filename.split("_")
         outf = open(pathout + "/" + "iterative_"+pars[1]+"_"+pars[2], "a")
         outf.write(line + "\n")
         #os.close(outf)
