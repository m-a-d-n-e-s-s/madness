
import math

for i in [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,71,77,83,89,95,101,113,119,125,131]:
    filename = "lebedev_%03d.txt" % i;
    data = open(filename).readlines()
    filename = "leb_%03d.txt" % i;
    f = open(filename,"w")
    f.write("%d\n" % len(data))
    for line in data:
        f.write(line)
    f.close()
