#!/usr/bin/env python3
import argparse
import csv
import re
import sys

SITE = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("barcodes")
    ap.add_argument("reference")
    ap.add_argument("--min-match", type=float, default=0.98)
    args = ap.parse_args()

    with open(args.barcodes, newline="") as fh:
        header = next(csv.reader(fh))
    sites = [(m.group(1), int(m.group(2))) for m in (SITE.match(h) for h in header[1:]) if m]
    if not sites:
        sys.stderr.write("no SNP barcode columns found\n")
        return 1
    best_name, best = None, 0.0
    for name, seq in read_fasta(args.reference).items():
        frac = sum(1 for base, pos in sites if pos <= len(seq) and seq[pos - 1] == base) / len(sites)
        if frac > best:
            best_name, best = name, frac
    if best < args.min_match:
        sys.stderr.write(f"reference does not match the barcode coordinates "
                         f"({100 * best:.1f}% of {len(sites)} sites match {best_name})\n")
        return 1
    print(best_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
