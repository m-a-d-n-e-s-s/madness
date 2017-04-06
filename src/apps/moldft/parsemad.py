import sys
import collections
import csv

def parsemad(filename):
    f = open(filename,'r')

    energies = []
    conviters = []
    timers = collections.defaultdict(float)
    nnodes = -1
    nthreads = -1
    time = -1
    iter = -1
    energy = -1
    nucrep = -1
    nalpha = nbeta = -1
    norbalpha = norbbeta = -1
    spin_restricted = ''
    aevals = []
    bevals = []
    dipole = -1
    grad = []

    while 1:
        line = f.readline()
        if line=='': break
        data = line.strip().split()

        if data == []:
            pass
        elif data[0] == 'number' and data[1] == 'of' and data[2] == 'electrons':
            nalpha = int(data[3])
            nbeta = int(data[4])
        elif data[0] == 'number' and data[1] == 'of' and data[2] == 'orbitals':
            norbalpha = int(data[3])
            norbbeta = int(data[4])
        elif data[0] == 'Total' and data[1] == 'Dipole' and data[2] == 'Moment:':
            dipole = float(data[3])
        elif data[0] == 'spin' and data[1] == 'restricted':
            spin_restricted = (data[2] == 'true')
        elif data[0] == 'Iteration':
            iter = int(data[1])
        elif data[0] == 'Converged!':
            conviters.append(iter)
            energies.append(energy)
        elif line[0:22] == '                total ':
            energy = float(data[1])
        elif data[0] == 'nuclear-repulsion':
            nucrep = float(data[1])
        elif data[0] == 'alpha' and data[1] == 'eigenvalues':
            data = map(float,f.readline()[3:].strip().split())
            aevals = data
        elif data[0] == 'beta' and data[1] == 'eigenvalues':
            data = map(float,f.readline()[3:].strip().split())
            bevals = data
        elif data[0] == '#nodes':
            nnodes = int(data[1])
        elif data[0] == '#total' and data[1] == 'threads':
            nthreads = int(data[2])
        elif line[0:24] == '         Total wall time':
            time = float(data[3][0:-2])
        elif data[0] == 'timer:':
            t = float(data[-1][0:-1])
            what = data[1]
            for s in data[2:-2]:
                what = what + ' ' + s
            timers[what] += t
        elif line[0:19] == ' Derivatives (a.u.)':
            grad = []
            f.readline(); # skip -----
            f.readline(); # skip blank
            f.readline(); # skip headers
            f.readline(); # skip -----
            while 1:
                line = f.readline()
                data = line.strip().split()
                if len(data) != 7: break
                grad.append(map(float,[data[4], data[5], data[6]]))

    # print "na, nb", nalpha, nbeta
    # print "norba, norbb", norbalpha, norbbeta
    # print "spin restricted", spin_restricted
    # print "nucrep", nucrep
    # print "energies", energies
    # print "iterations", conviters
    # print "nodes", nnodes
    # print "threads", nthreads
    # print "wall time", time
    # print "dipole", dipole
    # print "gradients", grad
    # print "alpha evals", aevals
    # if not spin_restricted:
    #     print " beta evals", bevals
    # print timers

    results = {}
    results["nalpha"] = nalpha
    results["nbeta"] = nbeta
    results["norba"] = norbalpha
    results["norbb"] = norbbeta
    results["spin restricted"] = spin_restricted
    results["nucrep"] = nucrep
    results["energies"] = energies
    results["iterations"] = conviters
    results["nodes"] = nnodes
    results["threads"] = nthreads
    results["wall time"] = time
    results["dipole"] = dipole
    results["gradients"] = grad
    results["alpha evals"] = aevals
    if not spin_restricted:
        results[" beta evals"] = bevals
    results["timers"] = timers

    return results

results = {}
for filename in sys.argv[1:]:
   results[filename] = parsemad(filename)

outf = open('out.csv', 'w')
out = csv.writer(outf)
F0 = results.keys()[0]
for tag in results[F0]['timers']:
    for filename in results.keys():
        nt = results[filename]['threads']
        TBB = (filename[0] == "T")
        tt = "Old"
        if TBB:
            nt = nt - 1
            tt = "TBB"
        out.writerow([tt,tag,nt,results[filename]['timers'][tag]])

for filename in results.keys():
    nt = results[filename]['threads']
    TBB = (filename[0] == "T")
    tt = "Old"
    if TBB:
        nt = nt - 1
        tt = "TBB"
    out.writerow([tt,'total',nt,results[filename]['wall time']])
    

   
