#!/usr/bin/env python3
import argparse
import math
import random

import numpy as np
import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam")
    ap.add_argument("output")
    ap.add_argument("--target", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    per_strand = max(1, math.ceil(args.target / 2))

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        lengths = list(bam.lengths)
        spans = []
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.reference_end is None:
                spans.append(None)
            else:
                spans.append((int(aln.is_reverse), aln.reference_id, aln.reference_start, aln.reference_end))

    depth = [[np.zeros(length, dtype=np.int32) for length in lengths] for _ in range(2)]
    order = [i for i, span in enumerate(spans) if span is not None]
    random.Random(args.seed).shuffle(order)
    keep = np.zeros(len(spans), dtype=bool)
    for i in order:
        strand, tid, start, end = spans[i]
        window = depth[strand][tid][start:end]
        if window.size and window.min() < per_strand:
            window += 1
            keep[i] = True

    with pysam.AlignmentFile(args.bam, "rb") as bam, \
            pysam.AlignmentFile(args.output, "wb", template=bam) as out:
        for i, aln in enumerate(bam.fetch(until_eof=True)):
            if keep[i]:
                out.write(aln)
    pysam.index(args.output)


if __name__ == "__main__":
    main()
