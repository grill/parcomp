import sys
import re

f = open(sys.argv[1])
line = " "
bla = float(sys.argv[2])

while(line != ""):
    line = f.readline()
    r = re.compile("^(\d+) (\d+\.\d+)$")
    m = r.match(line)
    if(m != None):
        print(m.group(1), float(m.group(2))/bla)
