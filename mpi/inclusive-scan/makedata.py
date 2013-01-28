import sys
import re

f = open(sys.argv[1])
line = " "
data = {}

while(line != ""):
    line = f.readline()
    r = re.compile("^inclusive-scan np=(\d+) s=\d+ time=(\d+\.\d+)$")
    m = r.match(line)
    if(m != None):
        i = int(m.group(1))
        if(i in data):
            data[i] = data[i] + [float(m.group(2))]
        else:
            data[i] = [float(m.group(2))]

for i in data:
    print(i, min(data[i]))
