#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from array import array

CLADE_ORDER = ["IIb", "Ib", "Ia", "IIa", "I", "II"]


def seqkit_stats(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    row = rows[0] if rows else {}

    def num(key, cast=int):
        try:
            return cast(float(row.get(key, "0")))
        except ValueError:
            return cast(0)

    return {"reads": num("num_seqs"), "bases": num("sum_len"), "mean_length": num("avg_len", float), "n50": num("N50")}


def read_counts(path):
    with open(path) as fh:
        return {k: int(v) for k, v in (line.strip().split("\t") for line in fh if "\t" in line)}


def read_fai(path):
    with open(path) as fh:
        return [(f[0], int(f[1])) for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) >= 2]


def read_depth(path, contigs):
    depth = {name: array("l", [0]) * length for name, length in contigs}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3 and f[0] in depth:
                pos = int(f[1]) - 1
                if 0 <= pos < len(depth[f[0]]):
                    depth[f[0]][pos] = int(f[2])
    return depth


def read_amplicons(path):
    with open(path) as fh:
        return [{"contig": f[0], "start": int(f[1]), "end": int(f[2]), "name": f[3]}
                for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) >= 4]


def window_stats(values, min_depth):
    n = len(values)
    if not n:
        return 0, 0.0, 0.0
    return n, round(sum(values) / n, 2), round(100.0 * sum(1 for d in values if d >= min_depth) / n, 2)


def merged_targets(amplicons):
    merged = {}
    for a in sorted(amplicons, key=lambda r: (r["contig"], r["start"])):
        ivs = merged.setdefault(a["contig"], [])
        if ivs and a["start"] <= ivs[-1][1]:
            ivs[-1][1] = max(ivs[-1][1], a["end"])
        else:
            ivs.append([a["start"], a["end"]])
    return merged


def consensus_stats(path):
    length = n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq = line.strip().upper()
                length += len(seq)
                n += seq.count("N")
    return length, n


def variant_counts(path):
    snps = indels = 0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["status"] == "PASS":
                if row["type"] in ("SNP", "MNP"):
                    snps += 1
                else:
                    indels += 1
    return snps, indels


def clade_of(lineage):
    name = re.sub(r"^[A-Za-z0-9]+-", "", lineage, count=1)
    name = re.sub(r"-like$", "", name)
    return name.split(".")[0] or lineage


def read_freyja(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            key, _, value = line.rstrip("\n").partition("\t")
            rows[key] = value
    names = rows.get("lineages", "").split()
    try:
        abundances = [float(x) for x in rows.get("abundances", "").split()]
    except ValueError:
        abundances = []
    lineages = sorted(zip(names, abundances), key=lambda x: -x[1])
    return rows.get("status", "ok"), lineages


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--raw-stats", required=True)
    ap.add_argument("--filtered-stats", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--depth", required=True)
    ap.add_argument("--amplicons", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--freyja")
    ap.add_argument("--min-depth", type=int, required=True)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    raw, filt = seqkit_stats(args.raw_stats), seqkit_stats(args.filtered_stats)
    counts = read_counts(args.counts)
    contigs = read_fai(args.fai)
    depth = read_depth(args.depth, contigs)
    amplicons = read_amplicons(args.amplicons)

    genome = [d for name, _ in contigs for d in depth[name]]
    genome_len, genome_mean, genome_cov = window_stats(genome, args.min_depth)
    target = [d for contig, ivs in merged_targets(amplicons).items() for s, e in ivs for d in depth[contig][s:e]]
    target_len, target_mean, target_cov = window_stats(target, args.min_depth)

    amplicon_rows = []
    for a in amplicons:
        _, mean, cov = window_stats(depth[a["contig"]][a["start"]:a["end"]], args.min_depth)
        amplicon_rows.append({"sample": args.sample, "amplicon": a["name"], "contig": a["contig"],
                              "start": a["start"] + 1, "end": a["end"], "mean_depth": mean,
                              "coverage_pct": cov})

    cons_len, cons_n = consensus_stats(args.consensus)
    snps, indels = variant_counts(args.variants)
    with open(args.model) as fh:
        model = fh.read().strip()
    mapped = counts.get("mapped_reads", 0)

    has_amplicons = bool(amplicons)
    summary = {
        "sample": args.sample,
        "total_reads": raw["reads"],
        "total_bases": raw["bases"],
        "mean_read_length": round(raw["mean_length"], 1),
        "read_n50": raw["n50"],
        "reads_after_length_filter": filt["reads"],
        "mapped_reads": mapped,
        "mapped_pct": round(100.0 * mapped / filt["reads"], 2) if filt["reads"] else 0.0,
        "reads_after_primer_trim": counts.get("trimmed_reads", ""),
        "min_depth": args.min_depth,
        "genome_length": genome_len,
        "genome_mean_depth": genome_mean,
        "genome_coverage_pct": genome_cov,
        "target_length": target_len if has_amplicons else "",
        "target_mean_depth": target_mean if has_amplicons else "",
        "target_coverage_pct": target_cov if has_amplicons else "",
        "amplicons_total": len(amplicons) if has_amplicons else "",
        "amplicons_below_min_depth": sum(1 for r in amplicon_rows if r["mean_depth"] < args.min_depth) if has_amplicons else "",
        "consensus_length": cons_len,
        "consensus_n_count": cons_n,
        "snps": snps,
        "indels": indels,
        "medaka_model": model,
    }

    if args.freyja:
        status, lineages = read_freyja(args.freyja)
        clades = {}
        for name, abundance in lineages:
            clades[clade_of(name)] = clades.get(clade_of(name), 0.0) + abundance
        ranked = sorted(clades.items(), key=lambda x: -x[1])
        summary["freyja_status"] = status
        summary["freyja_clades"] = "; ".join(f"{c} {100 * a:.1f}%" for c, a in ranked if a >= 0.001)
        summary["freyja_top_lineage"] = lineages[0][0] if lineages else ""
        known = [c for c in CLADE_ORDER if c in clades or any(n.startswith("MPX-") for n, _ in lineages)]
        for clade in known + sorted(set(clades) - set(known)):
            summary[f"freyja_{clade}_pct"] = round(100 * clades.get(clade, 0.0), 2)
        with open(f"{args.prefix}.lineages.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample", "lineage", "clade", "abundance_pct"])
            for name, abundance in lineages:
                w.writerow([args.sample, name, clade_of(name), round(100 * abundance, 3)])

    with open(f"{args.prefix}.summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary))
        w.writeheader()
        w.writerow(summary)
    if has_amplicons:
        with open(f"{args.prefix}.amplicons.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample", "amplicon", "contig", "start", "end", "mean_depth", "coverage_pct"])
            w.writeheader()
            w.writerows(amplicon_rows)


if __name__ == "__main__":
    main()
