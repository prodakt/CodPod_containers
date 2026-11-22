# KROK 1: Wczytanie i deduplikacja (sprawdzamy Ścieżkę A: extract_lnc = TRUE)
dane_surowe <- load_and_dedup(
  cds_file    = "Mus_musculus.GRCm39.cds.all.fa.gz",
  nc_file     = "Mus_musculus.GRCm39.ncrna.fa.gz",
  ex_file     = "merged_transcriptome.fasta", # lub NULL
  extract_lnc = TRUE,
  lnc_pattern = "lncRNA|lnc_rna",
  dedup_mode  = "exact_seq"
)

# KROK 2: Podgląd liczebności klas w konsoli
inspect_counts(dane_surowe)

# Jeśli widzisz, że 'lnc' jest za mało (np. < 600), uruchamiasz ponownie load_and_dedup z extract_lnc = FALSE.

# KROK 3: Podział z osobnymi proporcjami i limitami dla CP i ML
wynik <- split_datasets(
  df            = dane_surowe,
  prefix        = "Mus",
  seed          = 123,
  cp_stage_prop = 0.4,        # np. 40% danych idzie na CP, a 60% na ML
  cp_tt_prop    = 0.8,        # 80/20 w etapie CP
  ml_tt_prop    = 0.75,       # 75/25 w etapie ML
  max_n_cp      = 2000,       # max 2000 na klasę w CP
  max_n_ml      = 3000,       # max 3000 na klasę w ML
  balance_train = TRUE,
  ml_neg_mode   = "cds_only"
)

# KROK 4: Wyciąganie gotowych zbiorów sekwencji
cp_train <- extract_subset(wynik, stage_name = "cp", tt_name = "train")
cp_test  <- extract_subset(wynik, stage_name = "cp", tt_name = "test")
ml_train <- extract_subset(wynik, stage_name = "ml", tt_name = "train")
ml_test  <- extract_subset(wynik, stage_name = "ml", tt_name = "test")
ex_pred  <- extract_subset(wynik, stage_name = "inference")

# Szybka tabela kontrolna zachowanych próbek (keep == TRUE)
df_ok <- wynik$legend[wynik$legend$keep, ]
table(Etap_Zbior = paste(df_ok$stage, df_ok$tt_set, sep = "_"), Rola = df_ok$role)
