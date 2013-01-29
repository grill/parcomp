import sys
import re

f = open(sys.argv[1])
line = " "

while(line != ""):
    line = f.readline()
    r = re.compile("^(\d+) threads, (\d+.\d+), (\d+.\d+), (\d+.\d+), (\d+.\d+), (\d+.\d+)$")
    m = r.match(line)
    if(m != None):
        print("{0} {1}".format(m.group(1), min(map(float, m.groups()[1:]))))
