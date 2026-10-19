#!/usr/bin/env python3
import argparse
import re
import sys
from collections import defaultdict

NAME_RE = re.compile(
    r"^(?P<base>.+?)[_-](?P<num>\d+)[_-](?P<side>LEFT|RIGHT|L|R|F|FWD|FW|REV|RV)(?:[_-](?:alt\d*|\d+))?$",
    re.IGNORECASE,
)
LEFT_SIDES = {"LEFT", "L", "F", "FWD", "FW"}


def die(msg):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def warn(msg):
    sys.stderr.write(f"WARNING: {msg}\n")


def read_fai(path):
    with open(path) as fh:
        return {f[0]: int(f[1]) for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) >= 2}


def data_lines(path):
    with open(path, encoding="utf-8-sig") as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) == 1:
                fields = line.split()
            yield n, [x.strip() for x in fields]


def is_int(value):
    return re.fullmatch(r"-?\d+", value) is not None


def is_masterfile_row(f):
    return len(f) >= 10 and all(is_int(f[i]) for i in (1, 2, 4, 5, 7, 8))


def detect_format(path):
    for _, f in data_lines(path):
        if is_masterfile_row(f):
            return "masterfile"
        if len(f) >= 3 and is_int(f[1]) and is_int(f[2]):
            return "bed"
    die(f"{path}: not a primer BED or an xGen masterfile")


def interval(path, n, start, end, label):
    if start > end:
        warn(f"{path} line {n}: {label} has start > end; swapping")
        start, end = end, start
    if start == end:
        warn(f"{path} line {n}: {label} has zero length; skipping")
        return None
    return start, end


def parse_bed(path):
    primers = []
    for n, f in data_lines(path):
        if len(f) < 3 or not (is_int(f[1]) and is_int(f[2])):
            die(f"{path} line {n}: expected 'chrom start end name ...'")
        name = f[3] if len(f) > 3 and f[3] else f"primer_{n}"
        iv = interval(path, n, int(f[1]), int(f[2]), name)
        if iv is None:
            continue
        m = NAME_RE.match(name)
        strand = f[5] if len(f) > 5 and f[5] in ("+", "-") else ""
        if not strand:
            strand = ("+" if m.group("side").upper() in LEFT_SIDES else "-") if m else "."
        primers.append({"chrom": f[0], "start": iv[0], "end": iv[1], "name": name, "strand": strand})
    return primers


def parse_masterfile(path):
    primers, amplicons = [], []
    for n, f in data_lines(path):
        if not is_masterfile_row(f):
            if n == 1:
                continue
            die(f"{path} line {n}: expected an xGen masterfile row")
        chrom, amp = f[0], f[3] or f"amplicon_{len(amplicons) + 1}"
        fwd = interval(path, n, int(f[4]), int(f[5]), f"{amp} forward primer")
        rev = interval(path, n, int(f[7]), int(f[8]), f"{amp} reverse primer")
        ins = interval(path, n, int(f[1]), int(f[2]), f"{amp} insert")
        if fwd:
            primers.append({"chrom": chrom, "start": fwd[0], "end": fwd[1], "name": f[6] or f"{amp}_F", "strand": "+"})
        if rev:
            primers.append({"chrom": chrom, "start": rev[0], "end": rev[1], "name": f[9] or f"{amp}_R", "strand": "-"})
        if ins:
            amplicons.append({"chrom": chrom, "start": ins[0], "end": ins[1], "name": amp})
    return primers, amplicons


def pair_amplicons(primers):
    groups = defaultdict(lambda: {"left": [], "right": [], "base": set()})
    for p in primers:
        m = NAME_RE.match(p["name"])
        if not m:
            continue
        g = groups[(p["chrom"], int(m.group("num")))]
        g["left" if m.group("side").upper() in LEFT_SIDES else "right"].append(p)
        g["base"].add(m.group("base"))
    amplicons = []
    for (chrom, num), g in sorted(groups.items()):
        if not g["left"] or not g["right"]:
            continue
        start = max(p["end"] for p in g["left"])
        end = min(p["start"] for p in g["right"])
        if end > start:
            base = next(iter(g["base"])) if len(g["base"]) == 1 else "amplicon"
            amplicons.append({"chrom": chrom, "start": start, "end": end, "name": f"{base}_{num}"})
    return amplicons


def unique(records):
    seen, out = set(), []
    for r in records:
        key = (r["chrom"], r["start"], r["end"], r["name"])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


def validate(records, contigs, path):
    unknown = sorted({r["chrom"] for r in records} - set(contigs))
    if unknown:
        die(f"{path}: contig(s) {', '.join(unknown)} not in the reference ({', '.join(sorted(contigs)[:5])}); "
            "use the reference the primer scheme was designed on")
    for r in records:
        if r["start"] < 0 or r["end"] > contigs[r["chrom"]]:
            die(f"{path}: {r['name']} lies outside {r['chrom']} (length {contigs[r['chrom']]})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scheme")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--primers-out", required=True)
    ap.add_argument("--amplicons-out", required=True)
    args = ap.parse_args()

    contigs = read_fai(args.fai)
    if detect_format(args.scheme) == "masterfile":
        primers, amplicons = parse_masterfile(args.scheme)
    else:
        primers = parse_bed(args.scheme)
        amplicons = pair_amplicons(primers)
    primers, amplicons = unique(primers), unique(amplicons)
    if not primers:
        die(f"{args.scheme}: no primers found")
    if not amplicons:
        die(f"{args.scheme}: no amplicons could be paired from the primer names")
    validate(primers + amplicons, contigs, args.scheme)

    with open(args.primers_out, "w") as out:
        for p in sorted(primers, key=lambda r: (r["chrom"], r["start"], r["end"])):
            out.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t{p['name']}\t.\t{p['strand']}\n")
    with open(args.amplicons_out, "w") as out:
        for a in sorted(amplicons, key=lambda r: (r["chrom"], r["start"], r["end"])):
            out.write(f"{a['chrom']}\t{a['start']}\t{a['end']}\t{a['name']}\n")


if __name__ == "__main__":
    main()
