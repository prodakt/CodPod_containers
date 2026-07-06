#!/usr/bin/env python3
import sys
import gzip
import random

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

def read_fasta(path):
    """Czyta FASTA -> lista (header, seq)."""
    records = []
    with open_maybe_gzip(path, "rt") as fh:
        header = None
        seq_chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records

def write_fasta(path, records):
    """Zapisuje listę (header, seq) do FASTA."""
    with open(path, "wt") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")

def split_train_test(records, train_frac=0.6, seed=42):
    """Losowy podział records na train/test (train_frac / (1 - train_frac))."""
    random.seed(seed)
    indices = list(range(len(records)))
    random.shuffle(indices)
    n_train = int(round(train_frac * len(records)))
    train_idx = indices[:n_train]
    test_idx = indices[n_train:]
    train = [records[i] for i in train_idx]
    test = [records[i] for i in test_idx]
    return train, test

def main():
    if len(sys.argv) < 4:
        print(
            "Użycie: python codpot_split_train_test.py <species_tag> <pc_labeled.fa> <nc_labeled.fa> [train_frac] [seed]\n"
            "Domyślnie: train_frac=0.6, seed=42\n"
            "Przykład:\n"
            "  python codpot_split_train_test.py pig "
            "Sus_scrofa.Sscrofa11.1.cds.all_pc_labeled.fa "
            "Sus_scrofa.Sscrofa11.1.ncrna_nc_labeled.fa 0.6 42"
        )
        sys.exit(1)

    species   = sys.argv[1]   # np. human, mouse, pig
    pc_path   = sys.argv[2]
    nc_path   = sys.argv[3]
    train_frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
    seed       = int(sys.argv[5]) if len(sys.argv) > 5 else 42

    print(f"[INFO] Gatunek: {species}")
    print(f"[INFO] Czytam pcRNA labeled: {pc_path}")
    pc_all = read_fasta(pc_path)
    print(f"[INFO] pc: {len(pc_all)} sekwencji")

    print(f"[INFO] Czytam ncRNA labeled: {nc_path}")
    nc_all = read_fasta(nc_path)
    print(f"[INFO] nc: {len(nc_all)} sekwencji")

    print(f"[INFO] Dzielę pc na train/test (train_frac={train_frac})")
    pc_train, pc_test = split_train_test(pc_all, train_frac=train_frac, seed=seed)
    print(f"[INFO] pc train: {len(pc_train)}, pc test: {len(pc_test)}")

    print(f"[INFO] Dzielę nc na train/test (train_frac={train_frac})")
    nc_train, nc_test = split_train_test(nc_all, train_frac=train_frac, seed=seed + 1)
    print(f"[INFO] nc train: {len(nc_train)}, nc test: {len(nc_test)}")

    pc_train_out = f"pcRNA_{species}_train.fasta"
    pc_test_out  = f"pcRNA_{species}_test.fasta"
    nc_train_out = f"ncRNA_{species}_train.fasta"
    nc_test_out  = f"ncRNA_{species}_test.fasta"

    print("[INFO] Zapisuję pliki train/test")
    write_fasta(pc_train_out, pc_train)
    write_fasta(pc_test_out,  pc_test)
    write_fasta(nc_train_out, nc_train)
    write_fasta(nc_test_out,  nc_test)

    print("[INFO] Gotowe.")
    print("Outputy:")
    print(f"  {pc_train_out}")
    print(f"  {pc_test_out}")
    print(f"  {nc_train_out}")
    print(f"  {nc_test_out}")

if __name__ == "__main__":
    main()