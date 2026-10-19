#!/usr/bin/env python3
import argparse
import gzip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("depth")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--min-depth", type=int, required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.fai) as fh:
        contigs = [(f[0], int(f[1])) for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) >= 2]
    lengths = dict(contigs)
    intervals = {name: [] for name, _ in contigs}
    last = {}

    def add(contig, start, end):
        if end <= start:
            return
        ivs = intervals[contig]
        if ivs and ivs[-1][1] >= start:
            ivs[-1][1] = max(ivs[-1][1], end)
        else:
            ivs.append([start, end])

    with gzip.open(args.depth, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3 or f[0] not in lengths:
                continue
            contig, pos, depth = f[0], int(f[1]), int(f[2])
            prev = last.get(contig, 0)
            if pos > prev + 1:
                add(contig, prev, pos - 1)
            if depth < args.min_depth:
                add(contig, pos - 1, pos)
            last[contig] = max(prev, pos)

    with open(args.output, "w") as out:
        for contig, length in contigs:
            add(contig, last.get(contig, 0), length)
            for start, end in intervals[contig]:
                out.write(f"{contig}\t{start}\t{end}\n")


if __name__ == "__main__":
    main()
