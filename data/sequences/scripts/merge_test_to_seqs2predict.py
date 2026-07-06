#!/usr/bin/env python3
import sys
import random
from pathlib import Path

random.seed(42)


def read_fasta(path):
    records = []
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records


def write_fasta(records, path):
    with open(path, "w") as out:
        for h, s in records:
            out.write(f"{h}\n")
            out.write(f"{s}\n")


def sample_records(records, n):
    if n > len(records):
        n = len(records)
    return random.sample(records, n)


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: merge_tests_to_seqs2predict.py <nc_test.fa> <pc_test.fa>\n"
            "Example: merge_tests_to_seqs2predict.py ncRNA_Sus_scrofa_test.fasta pcRNA_Sus_scrofa_test.fasta\n"
        )
        sys.exit(1)

    nc_test_path = Path(sys.argv[1])
    pc_test_path = Path(sys.argv[2])

    if not nc_test_path.is_file() or not pc_test_path.is_file():
        sys.stderr.write(
            f"Missing {nc_test_path} or {pc_test_path} in current directory\n"
        )
        sys.exit(1)

    # Wczytaj pełne zbiory test
    nc_records = read_fasta(nc_test_path)
    pc_records = read_fasta(pc_test_path)

    count_nc = len(nc_records)
    count_pc = len(pc_records)
    sys.stderr.write(
        f"[INFO] nc_test: {count_nc} sequences, pc_test: {count_pc} sequences\n"
    )

    # === 1) seqs2predict_full: wszystkie testy, nagłówki seqX_nc / seqY_pc ===
    all_full = []
    idx = 1

    # najpierw nc, potem pc – ID rosną po kolei od seq1
    for _, seq in nc_records:
        seq_id = f"seq{idx}_nc"
        all_full.append((f">{seq_id}", seq))
        idx += 1
    for _, seq in pc_records:
        seq_id = f"seq{idx}_pc"
        all_full.append((f">{seq_id}", seq))
        idx += 1

    write_fasta(all_full, "seqs2predict_full.fasta")
    sys.stderr.write(
        "[INFO] wrote seqs2predict_full.fasta (all test sequences, renamed)\n"
    )

    # === 2) seqs2predict_small: 500 nc + 500 pc, nagłówki seqX_nc / seqY_pc ===
    target_nc_small = 500
    target_pc_small = 500

    if count_nc < target_nc_small:
        sys.stderr.write(
            f"[WARN] nc_test has only {count_nc} sequences, cannot sample 500; using {count_nc}\n"
        )
        target_nc_small = count_nc
    if count_pc < target_pc_small:
        sys.stderr.write(
            f"[WARN] pc_test has only {count_pc} sequences, cannot sample 500; using {count_pc}\n"
        )
        target_pc_small = count_pc

    sys.stderr.write(
        f"[INFO] seqs2predict_small targets: nc = {target_nc_small}, pc = {target_pc_small}\n"
    )

    nc_small = sample_records(nc_records, target_nc_small)
    pc_small = sample_records(pc_records, target_pc_small)

    all_small = []
    idx = 1

    # najpierw nc, potem pc – ID od seq1 w górę
    for _, seq in nc_small:
        seq_id = f"seq{idx}_nc"
        all_small.append((f">{seq_id}", seq))
        idx += 1
    for _, seq in pc_small:
        seq_id = f"seq{idx}_pc"
        all_small.append((f">{seq_id}", seq))
        idx += 1

    write_fasta(all_small, "seqs2predict_small.fasta")
    sys.stderr.write(
        "[INFO] wrote seqs2predict_small.fasta (500 nc + 500 pc, renamed)\n"
    )


if __name__ == "__main__":
    main()