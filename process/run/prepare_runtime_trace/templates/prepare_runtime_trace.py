#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict

RUNTIME_PROCESSES = {
    "DESEQ2",
    "LIMMA",
    "EDGER",
    "WILCOX",
    "TTEST",
    "MEMENTO",
    "NEBULA",
    "BOOTSTRAP2",
    "BOOTSTRAP3",
}

PROCESS_COLUMNS = [
    "deseq2",
    "limma",
    "edger",
    "wilcox",
    "ttest",
    "memento",
    "nebula",
    "bootstrap2",
    "bootstrap3",
]

PROCESS_NAME_TO_COLUMN = {
    "DESEQ2": "deseq2",
    "LIMMA": "limma",
    "EDGER": "edger",
    "WILCOX": "wilcox",
    "TTEST": "ttest",
    "MEMENTO": "memento",
    "NEBULA": "nebula",
    "BOOTSTRAP2": "bootstrap2",
    "BOOTSTRAP3": "bootstrap3",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Build runtime benchmark tables from a Nextflow trace.")
    parser.add_argument("--trace-file", required=True)
    parser.add_argument("--fixed-cells", default="fixed_cells.tsv")
    parser.add_argument("--fixed-genes", default="fixed_genes.tsv")
    return parser.parse_args()


def parse_realtime_seconds(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == "-":
        return None

    total = 0.0
    for number, unit in re.findall(r"([0-9]*\.?[0-9]+)\s*([a-zA-Z]+)", text):
        num = float(number)
        unit = unit.lower()
        if unit == "ms":
            total += num / 1e3
        elif unit in {"s", "sec", "secs", "second", "seconds"}:
            total += num
        elif unit in {"m", "min", "mins", "minute", "minutes"}:
            total += num * 60.0
        elif unit in {"h", "hr", "hrs", "hour", "hours"}:
            total += num * 3600.0
        elif unit in {"d", "day", "days"}:
            total += num * 86400.0
        else:
            return None

    if total > 0:
        return total

    try:
        return float(text)
    except ValueError:
        return None


def parse_tag(process_name, tag):
    parts = tag.split("_")
    if process_name not in RUNTIME_PROCESSES or len(parts) < 5:
        return None

    scenario = "_".join(parts[:-4])
    n_cells, n_genes, run, method = parts[-4], parts[-3], parts[-2], parts[-1]
    expected_method = PROCESS_NAME_TO_COLUMN[process_name]
    if method != expected_method:
        return None
    return scenario, n_cells, n_genes, run, method


def read_records(trace_file):
    records = defaultdict(list)
    with open(trace_file, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            process_full_name = (row.get("process") or row.get("name") or "").strip()
            if not process_full_name.startswith("RUNTIME:"):
                continue

            process_name = process_full_name.split(":")[-1]
            if process_name not in RUNTIME_PROCESSES:
                continue

            status = (row.get("status") or "").upper()
            if status not in {"COMPLETED", "CACHED"}:
                continue

            parsed = parse_tag(process_name, row.get("tag", ""))
            if parsed is None:
                continue

            realtime = row.get("realtime") or row.get("duration") or row.get("time")
            realtime_value = parse_realtime_seconds(realtime)
            if realtime_value is None:
                continue

            scenario, n_cells, n_genes, run, method = parsed
            records[(scenario, n_cells, n_genes, run, method)].append(realtime_value)

    return records


def build_rows(records, target_scenario, x_field):
    grouped = defaultdict(lambda: {column: [] for column in PROCESS_COLUMNS})

    for (scenario, n_cells, n_genes, _run, method), values in records.items():
        if scenario != target_scenario:
            continue
        grouped[(n_cells, n_genes)][method].extend(values)

    rows = []
    for (n_cells, n_genes), method_values in grouped.items():
        row = {"n_cells": int(n_cells), "n_genes": int(n_genes)}
        for column in PROCESS_COLUMNS:
            values = method_values[column]
            row[column] = sum(values) / len(values) if values else None
        rows.append(row)

    rows.sort(key=lambda item: item[x_field])
    return rows


def write_table(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["n_cells", "n_genes"] + PROCESS_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    records = read_records(args.trace_file)
    write_table(args.fixed_cells, build_rows(records, "fixed_cells", "n_genes"))
    write_table(args.fixed_genes, build_rows(records, "fixed_genes", "n_cells"))


if __name__ == "__main__":
    main()
