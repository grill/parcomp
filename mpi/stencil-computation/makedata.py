import sys
from sys import argv
import re

f = open(sys.argv[1])
line = " "
data = {}

def verify(k, r, c):
    return ((k == "1" and c == 1) or (k == "2" and r == 1) or (k == "3" and r == c) or \
            (k == "0" and c != 1 and r != 1 and r != c) or
            (k == "4" and c != 1 and c < r) or (k == "5" and r != 1 and r < c))

while(line != ""):
    line = f.readline()
    r = None
    r = re.compile("stencil computation, n={0}, m={1}, r=(\d+), c=(\d+), i=10, time=(\d+\.\d+)".format(argv[3], argv[4]))
    m = r.match(line)
    if(m != None):
        r = int(m.group(1))
        c = int(m.group(2))
        t = float(m.group(3))
        if(verify(argv[2], r, c)):
            if((r,c) in data):
                data[(r,c)] = data[(r,c)] + [t]
            else:
                data[(r,c)] = [t]

for i in data:
    print(i[0] * i[1], min(data[i]))

