#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
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
    parser.add_argument("--memory-runs", default="memory_runs.tsv")
    parser.add_argument("--memory-summary", default="memory_summary.tsv")
    return parser.parse_args()


def parse_realtime_seconds(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == "-":
        return None

    if re.fullmatch(r"[0-9]*\.?[0-9]+", text):
        return float(text) / 1000.0

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

    return None


def parse_memory_gib(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == "-":
        return None
    if re.fullmatch(r"[0-9]*\.?[0-9]+", text):
        return float(text) / 1024**3

    match = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*([KMGT]?B)", text, flags=re.IGNORECASE)
    if match is None:
        return None
    number = float(match.group(1))
    unit = match.group(2).upper()
    factors = {"B": 1 / 1024**3, "KB": 1 / 1024**2, "MB": 1 / 1024, "GB": 1, "TB": 1024}
    return number * factors[unit]


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
    memory_records = defaultdict(list)
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
            scenario, n_cells, n_genes, run, method = parsed
            key = (scenario, n_cells, n_genes, run, method)
            if realtime_value is not None:
                records[key].append(realtime_value)
            peak_rss_value = parse_memory_gib(row.get("peak_rss"))
            if peak_rss_value is not None:
                memory_records[key].append(peak_rss_value)

    return records, memory_records


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


def build_memory_rows(records):
    parsed = []
    for (scenario, n_cells, n_genes, run, method), values in records.items():
        if len(values) != 1:
            raise ValueError(
                "Expected one completed method task for "
                f"{scenario}/{n_cells}/{n_genes}/run={run}/{method}, found {len(values)}"
            )
        parsed.append({
            "scenario": scenario,
            "n_cells": int(n_cells),
            "n_genes": int(n_genes),
            "run": int(run),
            "method": method,
            "peak_rss_gib": values[0],
        })

    parsed.sort(key=lambda row: (
        row["scenario"], row["n_genes"], row["n_cells"],
        PROCESS_COLUMNS.index(row["method"]), row["run"],
    ))
    return parsed


def write_memory_tables(runs_path, summary_path, rows):
    run_fields = ["scenario", "n_cells", "n_genes", "run", "method", "peak_rss_gib"]
    with open(runs_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=run_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grouped = defaultdict(list)
    for row in rows:
        key = (row["scenario"], row["n_cells"], row["n_genes"], row["method"])
        grouped[key].append(row["peak_rss_gib"])

    summary_rows = []
    for (scenario, n_cells, n_genes, method), values in grouped.items():
        summary_rows.append({
            "scenario": scenario,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "method": method,
            "mean_peak_rss_gib": sum(values) / len(values),
            "sd_peak_rss_gib": statistics.stdev(values) if len(values) > 1 else 0.0,
            "n_runs": len(values),
        })
    summary_rows.sort(key=lambda row: (
        row["scenario"], row["n_genes"], row["n_cells"],
        PROCESS_COLUMNS.index(row["method"]),
    ))
    summary_fields = [
        "scenario", "n_cells", "n_genes", "method",
        "mean_peak_rss_gib", "sd_peak_rss_gib", "n_runs",
    ]
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


def main():
    args = parse_args()
    records, memory_records = read_records(args.trace_file)
    write_table(args.fixed_cells, build_rows(records, "fixed_cells", "n_genes"))
    write_table(args.fixed_genes, build_rows(records, "fixed_genes", "n_cells"))
    write_memory_tables(
        args.memory_runs,
        args.memory_summary,
        build_memory_rows(memory_records),
    )


if __name__ == "__main__":
    main()
