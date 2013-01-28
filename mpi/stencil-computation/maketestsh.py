#!/usr/bin/python

i = 1
print("./jupiterrun.py " + str(i) + " 29 stencil-computation -n 1770 -m 1770 -r " + str(i) + " -c 1 -i 10 >> test.log")
for i in range(2, 49, 2):
   print("./jupiterrun.py " + str(i) + " 29 stencil-computation -n 1770 -m 1770 -r " + str(i) + " -c 1 -i 10 >> test.log")

for i in range(2, 49, 2):
   print("./jupiterrun.py " + str(i) + " 29 stencil-computation -n 1770 -m 1770 -r 1 -c " + str(i) + " -i 10 >> test.log")

for i in range(2, 9, 2):
    for j in range(2, 9, 2):
        print("./jupiterrun.py " + str(i*j) + " 29 stencil-computation -n 1770 -m 1770 -r "+str(i)+" -c " + str(j) + " -i 10 >> test.log")
