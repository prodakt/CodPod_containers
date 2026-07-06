#!/usr/bin/env python3
import sys
import os
import gzip

def open_maybe_gzip(path, mode="rt"):
    """Otwiera zwykły plik lub .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

def read_fasta(path):
    """
    Czyta FASTA i zwraca listę (header, seq).
    Nagłówek bez znaku '>'.
    """
    records = []
    with open_maybe_gzip(path, "rt") as fh:
        header = None
        seq_chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # zapisz poprzedni rekord
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:]  # bez '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # ostatni rekord
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records

def write_fasta(path, records):
    """Zapisuje listę (header, seq) do FASTA."""
    with open(path, "wt") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            # można łamać na linie 60 nt, ale prościej całość w jednej
            fh.write(f"{seq}\n")

def label_records(records, label_suffix):
    """
    Dokleja sufiks do ID w nagłówku.
    Jeśli nagłówek ma formę 'ID cośtam', zachowujemy tylko ID.
    """
    labeled = []
    for header, seq in records:
        # weź pierwszy token jako ID
        parts = header.split()
        base_id = parts[0]
        new_header = f"{base_id}{label_suffix}"
        labeled.append((new_header, seq))
    return labeled

def main():
    if len(sys.argv) != 3:
        print(
            "Użycie: python label_fasta_pc_nc.py <cds_fasta> <nc_fasta>\n"
            "Przykład:\n"
            "  python label_fasta_pc_nc.py Sus_scrofa.Sscrofa11.1.cds.all.fa Sus_scrofa.Sscrofa11.1.ncrna.fa"
        )
        sys.exit(1)

    cds_path = sys.argv[1]
    nc_path = sys.argv[2]

    # nazwy outputów (możesz dopasować do własnego schematu)
    cds_out = os.path.splitext(cds_path)[0] + "_pc_labeled.fa"
    nc_out  = os.path.splitext(nc_path)[0] + "_nc_labeled.fa"

    print(f"[INFO] Czytam CDS z: {cds_path}")
    cds_records = read_fasta(cds_path)
    print(f"[INFO] Znaleziono {len(cds_records)} sekwencji CDS")

    print(f"[INFO] Czytam ncRNA z: {nc_path}")
    nc_records = read_fasta(nc_path)
    print(f"[INFO] Znaleziono {len(nc_records)} sekwencji ncRNA")

    # Doklejanie sufiksów
    cds_labeled = label_records(cds_records, "_pc")
    nc_labeled  = label_records(nc_records, "_nc")

    print(f"[INFO] Zapisuję CDS z sufiksami '_pc' do: {cds_out}")
    write_fasta(cds_out, cds_labeled)

    print(f"[INFO] Zapisuję ncRNA z sufiksami '_nc' do: {nc_out}")
    write_fasta(nc_out, nc_labeled)

    print("[INFO] Gotowe.")

if __name__ == "__main__":
    main()