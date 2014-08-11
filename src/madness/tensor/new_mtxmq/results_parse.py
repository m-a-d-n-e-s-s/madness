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

# { (shape, i, j, k) : [(variant, gflops), ] , }

def get_data():
    results = {}
    for line in fileinput.input():
        line = line.split()
        test = tuple([line[1]] + [int(x) for x in line[2:-1]])
        gflops = float(line[-1])
        version = line[0][:-1] # Remove semicolon
        if test not in results:
            results[test] = []
        results[test].append((version, gflops))
    return results

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
    results = get_data()
    print(COLORS.HEADER + "Best per benchmark:" + COLORS.ENDC)
    print_best_per_benchmark(results)
    print(COLORS.HEADER + "Top 10 Varaints:" + COLORS.ENDC)
    print_top_n_variants(results)
