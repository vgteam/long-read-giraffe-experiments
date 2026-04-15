#!/usr/bin/env python3
"""Combine `vg stats -a` reports (and optionally a total-errors TSV) into a
headered TSV on stdout.

Each positional argument is a file containing the output of `vg stats -a` for
one alignment condition. The condition name for a row is the input file's
basename with any trailing ".stats" stripped (the benchmark harness writes
these reports with a .stats suffix by convention).

If --errors is given, it is read as a two-column condition<TAB>total_errors
TSV; conditions without a matching row get an empty total_errors field.
"""

import argparse
import re
import sys
from pathlib import Path


COLUMNS = [
    "condition",
    "total_errors",
    "total_aligned",
    "total_perfect",
    "total_gapless",
    "align_score_mean",
    "align_score_median",
    "align_score_stdev",
    "mapq_mean",
    "mapq_median",
    "mapq_stdev",
    "mapq_max_count",
    "insertions_bp",
    "insertions_events",
    "deletions_bp",
    "deletions_events",
    "substitutions_bp",
    "substitutions_events",
    "matches_bp",
    "matches_bp_per_aln",
    "softclips_bp",
    "softclips_pct",
    "softclips_bp_per_aln",
    "softclips_events",
    "total_time_sec",
    "speed_reads_per_sec",
]


# Each entry maps a regex to the output columns its capture groups fill, in
# order. Keeping the regexes together lets a new `vg stats -a` line be picked
# up by appending one row rather than extending an if/elif chain.
LINE_PARSERS = [
    (re.compile(r"^Total aligned:\s*(\S+)"), ["total_aligned"]),
    (re.compile(r"^Total perfect:\s*(\S+)"), ["total_perfect"]),
    (re.compile(r"^Total gapless[^:]*:\s*(\S+)"), ["total_gapless"]),
    (
        re.compile(
            r"^Alignment score: mean (\S+?), median (\S+?), stdev (\S+?), max"
        ),
        ["align_score_mean", "align_score_median", "align_score_stdev"],
    ),
    (
        re.compile(
            r"^Mapping quality: mean (\S+?), median (\S+?), stdev (\S+?),"
            r" max \S+ \((\S+) reads\)"
        ),
        ["mapq_mean", "mapq_median", "mapq_stdev", "mapq_max_count"],
    ),
    (
        re.compile(r"^Insertions:\s*(\S+)\s*bp in\s*(\S+)\s*read events"),
        ["insertions_bp", "insertions_events"],
    ),
    (
        re.compile(r"^Deletions:\s*(\S+)\s*bp in\s*(\S+)\s*read events"),
        ["deletions_bp", "deletions_events"],
    ),
    (
        re.compile(r"^Substitutions:\s*(\S+)\s*bp in\s*(\S+)\s*read events"),
        ["substitutions_bp", "substitutions_events"],
    ),
    (
        re.compile(r"^Matches:\s*(\S+)\s*bp \((\S+)\s*bp/alignment\)"),
        ["matches_bp", "matches_bp_per_aln"],
    ),
    (
        re.compile(
            r"^Softclips:\s*(\S+)\s*bp \((\S+)% of bases,"
            r"\s*(\S+)\s*bp/alignment\) in\s*(\S+)\s*read events"
        ),
        [
            "softclips_bp",
            "softclips_pct",
            "softclips_bp_per_aln",
            "softclips_events",
        ],
    ),
    (re.compile(r"^Total time:\s*(\S+)"), ["total_time_sec"]),
    (re.compile(r"^Speed:\s*(\S+)"), ["speed_reads_per_sec"]),
]


def parse_stats(path):
    fields = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            for pattern, cols in LINE_PARSERS:
                m = pattern.match(line)
                if m:
                    for col, val in zip(cols, m.groups()):
                        fields[col] = val
                    break
    return fields


def load_errors(path):
    errors = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                errors[parts[0]] = parts[1]
    return errors


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "stats",
        nargs="+",
        help="files containing `vg stats -a` output (one per condition)",
    )
    ap.add_argument(
        "--errors",
        help="TSV of condition<TAB>total_errors to join against the stats rows",
    )
    args = ap.parse_args()

    errors = load_errors(args.errors) if args.errors else {}

    out = sys.stdout
    out.write("\t".join(COLUMNS) + "\n")
    for stats_path in args.stats:
        condition = Path(stats_path).name
        if condition.endswith(".stats"):
            condition = condition[: -len(".stats")]
        fields = parse_stats(stats_path)
        fields["condition"] = condition
        fields["total_errors"] = errors.get(condition, "")
        out.write("\t".join(fields.get(col, "") for col in COLUMNS) + "\n")


if __name__ == "__main__":
    main()
