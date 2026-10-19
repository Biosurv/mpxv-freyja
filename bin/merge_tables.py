#!/usr/bin/env python3
import csv
import sys


def main():
    output, inputs = sys.argv[1], sys.argv[2:]
    columns, rows = [], []
    for path in inputs:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for name in reader.fieldnames or []:
                if name not in columns:
                    columns.append(name)
            rows.extend(reader)
    clade_columns = [c for c in columns if c.startswith("freyja_") and c.endswith("_pct")]
    for row in rows:
        if row.get("freyja_status") == "ok":
            for c in clade_columns:
                row.setdefault(c, "0.0")
    with open(output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, restval="")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
