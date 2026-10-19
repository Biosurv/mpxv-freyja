#!/usr/bin/env python3
import argparse
import bisect
from collections import defaultdict


def load_mask(path):
    starts, ends = defaultdict(list), defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3:
                starts[f[0]].append(int(f[1]))
                ends[f[0]].append(int(f[2]))
    return starts, ends


def masked(mask, contig, start, end):
    starts, ends = mask
    s = starts.get(contig)
    if not s:
        return False
    i = bisect.bisect_right(s, end - 1) - 1
    return i >= 0 and ends[contig][i] > start


def variant_type(ref, alt):
    if len(ref) == len(alt):
        return "SNP" if len(ref) == 1 else "MNP"
    return "INS" if len(alt) > len(ref) else "DEL"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("--mask", required=True)
    ap.add_argument("--min-depth", type=int, required=True)
    ap.add_argument("--min-qual", type=float, required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-mask", required=True)
    args = ap.parse_args()

    mask = load_mask(args.mask)
    seen = set()
    with open(args.vcf) as vin, open(args.out_vcf, "w") as vout, open(args.out_tsv, "w") as tout, \
            open(args.out_mask, "w") as mout:
        tout.write("contig\tpos\tref\talt\tqual\tdepth\ttype\tstatus\n")
        for line in vin:
            if line.startswith("#"):
                vout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or (f[0], f[1], f[3], f[4]) in seen:
                continue
            seen.add((f[0], f[1], f[3], f[4]))
            contig, pos, ref, alt = f[0], int(f[1]), f[3], f[4].split(",")[0]
            info = dict(item.partition("=")[::2] for item in f[7].split(";") if item)
            depth = int(float(info["DP"])) if info.get("DP") not in (None, "", ".") else 0
            qual = float(f[5]) if f[5] not in ("", ".") else 0.0
            reasons = []
            if f[6] not in ("PASS", "."):
                reasons.append(f[6])
            if depth < args.min_depth:
                reasons.append("low_depth")
            if qual < args.min_qual:
                reasons.append("low_qual")
            if masked(mask, contig, pos - 1, pos - 1 + len(ref)):
                reasons.append("masked")
            if not reasons:
                vout.write(line if line.endswith("\n") else line + "\n")
            elif reasons == ["low_qual"]:
                mout.write(f"{contig}\t{pos - 1}\t{pos - 1 + len(ref)}\n")
            status = "PASS" if not reasons else ";".join(reasons)
            tout.write(f"{contig}\t{pos}\t{ref}\t{alt}\t{qual:.2f}\t{depth}\t{variant_type(ref, alt)}\t{status}\n")


if __name__ == "__main__":
    main()
