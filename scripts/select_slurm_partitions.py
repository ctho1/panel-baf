#!/usr/bin/env python3
"""Rank allowed SLURM CPU partitions by their currently idle CPUs."""

import argparse
import subprocess
import sys


def parse_sinfo(text, candidates):
    """Return available candidates ranked by aggregate idle CPU count."""
    candidate_order = {name: index for index, name in enumerate(candidates)}
    idle_cpus = {name: 0 for name in candidates}
    available = set()

    for line in text.splitlines():
        fields = [field.strip() for field in line.split("|")]
        if len(fields) != 3:
            continue
        partition = fields[0].rstrip("*")
        if partition not in candidate_order or fields[1].lower() != "up":
            continue

        cpu_states = fields[2].split("/")
        if len(cpu_states) != 4:
            continue
        try:
            idle = int(cpu_states[1])
        except ValueError:
            continue

        available.add(partition)
        idle_cpus[partition] += idle

    ranked = sorted(
        (name for name in available if idle_cpus[name] > 0),
        key=lambda name: (-idle_cpus[name], candidate_order[name]),
    )
    return ranked, idle_cpus


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--partitions", required=True,
                        help="Comma-separated allowlist of public CPU partitions")
    parser.add_argument("--fallback", default="normal")
    parser.add_argument("--sinfo", default="sinfo")
    args = parser.parse_args(argv)

    candidates = [name.strip() for name in args.partitions.split(",") if name.strip()]
    if not candidates:
        print(args.fallback)
        print("Keine SLURM-Kandidaten konfiguriert; nutze Fallback.", file=sys.stderr)
        return

    try:
        result = subprocess.run(
            [
                args.sinfo,
                "--noheader",
                f"--partition={','.join(candidates)}",
                "--format=%P|%a|%C",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        ranked, idle_cpus = parse_sinfo(result.stdout, candidates)
    except (OSError, subprocess.CalledProcessError) as exc:
        print(args.fallback)
        print(f"SLURM-Auslastung nicht abfragbar ({exc}); nutze {args.fallback}.",
              file=sys.stderr)
        return

    if not ranked:
        print(args.fallback)
        print(f"Keine aktive Kandidatenpartition gefunden; nutze {args.fallback}.",
              file=sys.stderr)
        return

    summary = ", ".join(f"{name}: {idle_cpus[name]} freie CPUs" for name in ranked)
    print(f"SLURM-Auslastung: {summary}", file=sys.stderr)
    print(",".join(ranked))


if __name__ == "__main__":
    main()
