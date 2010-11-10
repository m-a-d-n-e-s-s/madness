#!/usr/bin/env python
import sys, re
try:
    sys.argv[1]
except:
    print "Which file?"
    sys.exit()

#Parse MADNESS output 
inputFile = sys.argv[1]
print inputFile
f = open(inputFile, 'r')
exFile = open("ex.dat", 'w')
ionFile = open("ion.dat", 'w')
RlFile = open("Rl.dat", 'w')
lines = f.readlines()
while 1:
    if(lines):
        line = lines.pop(0)
        exState = re.match( "^[1-9]", line)
        kState = re.match( "^0", line)
        Rl = re.match( "^Y(\d+)0", line)
        if exState:
            words = line.split()
            for word in words:
                exFile.write(word + "\t")
                exFile.write("\n")
        if kState:
            words = line.split()
            words.pop(0)
            for word in words:
                if( word == 'took'):
                    break
                ionFile.write(word + "\t")
            ionFile.write("\n")
        if Rl:
            words = line.split()
            l = Rl.group(1)
            Pl = lines.pop(0)[0:-1]
            RlFile.write( l + "\t" + Pl + "\t" + " ".join(words[1:]) + "\n" )
    else:
        break
f.close()

#Get last line of wf.num
f = open("wf.num", "r")
lines = f.readlines()
maxStep = int(lines.pop())
f.close()
print "maxStep = ", maxStep

#Sort STEP, atomic time, and walltime
f = open("time.dat", 'r')
lines = f.readlines()
legacyTime = 0
lastTime = 0
lastStep = 0
fout = open("t.dat", "w")
for line in lines:
    time = line.split()
    step  = float(time[0])
    aTime = time[1]
    wTime = float(time[2]) # wall time
    if step <= lastStep:   # no repeat time steps
        continue
    if step == maxStep:
        tFile = open("tMAX.dat", 'w')
        tFile.write( aTime )
        tFile.close()
    if wTime < lastTime:   # add to legacy time
        legacyTime += lastTime
    totalTime = wTime + legacyTime
    fout.write(time[0] + "\t" + time[1] + "\t" + str(totalTime) + "\n")
    lastStep = step
    lastTime = wTime
f.close()

#make dr.dat
f = open("input", 'r')
lines = f.readlines()
for line in lines:
    if line:
        word = line.split()
        if word[0] == 'L':
            L = float(word[1])
f.close()
f = open("input2", 'r')
lines = f.readlines()
for line in lines:
    if line:
        word = line.split()
        if word[0] == 'n':
            n = int(word[1])
f.close()
f = open("dr.dat", 'w')
f.write( str(L/n) )
f.close()
