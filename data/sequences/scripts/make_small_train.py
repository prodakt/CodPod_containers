#!/usr/bin/env python3

import argparse
import random
from pathlib import Path


def read_fasta(path):
    records = []
    header = None
    seq_lines = []

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
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


def write_fasta(records, path, wrap=80):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i+wrap] + "\n")


def sample_n_records(records, n, seed):
    if n > len(records):
        raise ValueError(
            f"Requested {n} sequences, but file contains only {len(records)} sequences."
        )
    rng = random.Random(seed)
    return rng.sample(records, n)


def make_output_name(input_path):
    p = Path(input_path)
    if p.stem.endswith("_train"):
        out_name = p.stem + "_small" + p.suffix
    else:
        out_name = p.stem + "_train_small" + p.suffix
    return p.with_name(out_name)


def main():
    parser = argparse.ArgumentParser(
        description="Randomly sample N sequences from pc_train and nc_train FASTA files."
    )
    parser.add_argument("pc_fasta", help="Input PC train FASTA")
    parser.add_argument("nc_fasta", help="Input NC train FASTA")
    parser.add_argument("seed", type=int, help="Random seed")
    parser.add_argument("n_seqs", type=int, help="Number of sequences to sample")

    args = parser.parse_args()

    pc_path = Path(args.pc_fasta)
    nc_path = Path(args.nc_fasta)

    if not pc_path.exists():
        raise FileNotFoundError(f"PC FASTA not found: {pc_path}")
    if not nc_path.exists():
        raise FileNotFoundError(f"NC FASTA not found: {nc_path}")

    pc_records = read_fasta(pc_path)
    nc_records = read_fasta(nc_path)

    pc_selected = sample_n_records(pc_records, args.n_seqs, args.seed)
    nc_selected = sample_n_records(nc_records, args.n_seqs, args.seed + 1)

    pc_out = make_output_name(pc_path)
    nc_out = make_output_name(nc_path)

    write_fasta(pc_selected, pc_out)
    write_fasta(nc_selected, nc_out)

    print(f"PC input : {pc_path}")
    print(f"PC total : {len(pc_records)}")
    print(f"PC output: {pc_out}")
    print(f"PC wrote : {len(pc_selected)} sequences")
    print()
    print(f"NC input : {nc_path}")
    print(f"NC total : {len(nc_records)}")
    print(f"NC output: {nc_out}")
    print(f"NC wrote : {len(nc_selected)} sequences")


if __name__ == "__main__":
    main()