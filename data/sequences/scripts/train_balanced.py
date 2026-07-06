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
            "Usage: train_balanced.py <nc_train.fa> <pc_train.fa>\n"
            "Example: train_balanced.py ncRNA_Sus_scrofa_train.fasta pcRNA_Sus_scrofa_train.fasta\n"
        )
        sys.exit(1)

    nc_train_path = Path(sys.argv[1])
    pc_train_path = Path(sys.argv[2])

    if not nc_train_path.is_file() or not pc_train_path.is_file():
        sys.stderr.write(
            f"Missing {nc_train_path} or {pc_train_path} in current directory\n"
        )
        sys.exit(1)

    nc_records = read_fasta(nc_train_path)
    pc_records = read_fasta(pc_train_path)

    count_nc = len(nc_records)
    count_pc = len(pc_records)
    sys.stderr.write(
        f"[INFO] nc_train: {count_nc} sequences, pc_train: {count_pc} sequences\n"
    )

    min_class = min(count_nc, count_pc)
    sys.stderr.write(f"[INFO] minority class size (balanced_full) = {min_class}\n")

    bal_nc = sample_records(nc_records, min_class)
    bal_pc = sample_records(pc_records, min_class)

    bal_nc_path = Path("nc_train_balanced.fasta")
    bal_pc_path = Path("pc_train_balanced.fasta")

    write_fasta(bal_nc, bal_nc_path)
    write_fasta(bal_pc, bal_pc_path)
    sys.stderr.write(
        f"[INFO] wrote balanced trains: {bal_nc_path}, {bal_pc_path}\n"
    )


if __name__ == "__main__":
    main()