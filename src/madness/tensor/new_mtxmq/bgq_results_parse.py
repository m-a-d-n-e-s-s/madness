#!/usr/bin/env python3

import fileinput

class COLORS:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

# [ (variant, pthreads, (shape, i, j, k), gflops), ]

def get_data():
    results = []
    variants = set()
    tests = set()
    pthreads = set()
    for line in fileinput.input():
        line = line.strip().replace('[', '').replace('pthreads]:', '').replace('(avg', 'avg').replace(':', '')[:-1].split()[:-6]
        test = tuple([line[2]] + [int(x) for x in line[3:-1]])
        gflops = float(line[-1])
        results.append((line[0], int(line[1]), test, gflops))
        variants.add(line[0])
        tests.add(test)
        pthreads.add(int(line[1]))
    return (results, sorted(variants), sorted(pthreads), sorted(tests))

def print_best_per_benchmark(results):
    best_avg, blas_avg, imp_avg = 0, 0, 0
    for test in sorted(results):
        best = max(results[test], key=lambda x: x[1])
        best_avg += best[1]
        color = COLORS.OKGREEN
        print(color + '{:40}'.format(' '.join([str(x) for x in test])) + COLORS.ENDC, '{:40}{:6.2f}'.format(*best))
    best_avg /= len(results)
    print('{:40} {:40}{:6.2f}'.format('', 'average', best_avg))

def print_top_n_variants(results, n=10):
    variants = {}
    for test in results:
        for variant, gflops in results[test]:
            if variant not in variants:
                variants[variant] = 0
            variants[variant] += gflops / len(results)

    topn = sorted(variants, key=lambda x: variants[x], reverse=True)[:n]
    for x in topn:
        print('{:40}{:6.2f}'.format(x, variants[x]))


if __name__ == '__main__':
    results, variants, pthreads, tests = get_data()
    for pc in pthreads:
        pcfilter = list(filter(lambda x: x[1] == pc, results))
        for test in tests:
            tfilter = sorted(filter(lambda x: x[2] == test, pcfilter), key=lambda x: x[3])
            if tfilter:
                print(pc, test, tfilter[-1][0], tfilter[-1][3])
        avgs = {}
        for variant in variants:
            vfilter = list(filter(lambda x: x[0] == variant, pcfilter))
            avg = sum(x[3] for x in vfilter)/len(vfilter)
            avgs[variant] = avg
        top = sorted(avgs, key=lambda x: avgs[x])[-5:][::-1]
        for v in top:
            print(pc, v, avgs[v])


    #print(COLORS.HEADER + "Best per benchmark:" + COLORS.ENDC)
    #print_best_per_benchmark(results)
    #print(COLORS.HEADER + "Top 10 Varaints:" + COLORS.ENDC)
    #print_top_n_variants(results)
