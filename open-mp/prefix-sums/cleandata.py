import re
import os
import sys

path = sys.argv[1]
pathout = sys.argv[2]
dirs = os.listdir(path)
for filename in dirs:
   inf = open(path + "/" + filename)
   outf = open(pathout + "/" + filename, "w+")
   line = inf.readline()
   all = line.split(';')
   for x in all[1:len(all)-1]:
      outf.write(x + "\n")
