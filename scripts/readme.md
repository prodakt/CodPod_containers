# Dwuetapowy Pipeline Identyfikacji lncRNA (`CP` -> `ML` -> `EX`): Metodologia Podziału Danych i Skrypt R

Dokumentacja modułu **`dataset preparation`** przeznaczonego do przygotowania zbiorów transkryptomicznych (RNA-seq) do dwuetapowej klasyfikacji długich niekodujących RNA (`lncRNA`).

---

## 1. Architektura i Założenia Metodologiczne

Pipeline składa się z dwóch zależnych od siebie etapów modelowania oraz etapu predykcji na danych eksperymentalnych:
1. **Etap 1 (`cp` – Coding Potential):** Ocena potencjału kodowania sekwencji za pomocą zestawu narzędzi bioinformatycznych.
2. **Etap 2 (`ml` – Machine Learning):** Meta-klasyfikator łączący wyniki predykcji z etapu `cp` z cechami strukturalnymi i sekwencyjnymi transkryptów (długość, liczba eksonów, obecność ogona polyA itp.).
3. **Etap 3 (`inference` – `ex`):** Identyfikacja nowych, niezaadnotowanych transkryptów `lncRNA` ze złożenia transkryptomu (np. StringTie `merged_transcriptome.fasta`).

### Kluczowe zasady czystości danych (Zapobieganie *Data Leakage*)
* **Rozłączność bloków `cp` i `ml`:** Ponieważ wyniki z narzędzi `cp` stanowią cechy wejściowe dla modelu `ml`, oba etapy muszą być trenowane i testowane na całkowicie rozłącznych pulach transkryptów. Trenowanie `ml` na danych widzianych wcześniej przez `cp` prowadziłoby do sztucznego zawyżenia skuteczności cech z pierwszego etapu (*data leakage*).
* **Czyszczenie zbioru eksperymentalnego (`ex`):** Zbiór `ex` zawiera mieszaninę znanych transkryptów kodujących, znanych niekodujących oraz transkryptów nowych. Przed predykcją musi zostać oczyszczony z wszelkich rekordów identycznych z referencją (`cds`, `nc`, `lnc`) i nigdy nie bierze udziału w treningu ani walidacji.
* **Analiza na poziomie transkryptów:** Ze względu na zjawisko alternatywnego splicingu (wiele genów kodujących białka posiada niekodujące izoformy, np. *retained intron* lub *processed transcript*), podział wykonywany jest na poziomie pojedynczych transkryptów z zachowaniem rygorystycznej deduplikacji.

---

## 2. Podwójna Ścieżka Przetwarzania (`Ścieżka A` vs `Ścieżka B`)

W zależności od jakości adnotacji badanego gatunku (organizmy modelowe vs niemodelowe) pipeline oferuje dwie ścieżki konfiguracji sterowane parametrem `extract_lnc`:

* **Ścieżka A (`extract_lnc = TRUE` – organizmy dobrze zaadnotowane):**
  * Z ogólnego pliku niekodującego (`nc_file`) wyodrębniana jest podklasa `lnc` (na podstawie wzorca regex `lnc_pattern` lub osobnego pliku `lnc_file`), a pozostałe transkrypty oznaczane są jako `nc_other`.
  * W etapie `cp` model uczy się odróżniać pełną pulę niekodującą (`lnc` + `nc_other`) od kodującej (`cds`), co zachowuje naturalny rozkład biologiczny.
  * W etapie `ml` klasą pozytywną jest wyłącznie `lnc` (przeciwko `cds` lub `cds + nc_other`).
* **Ścieżka B (`extract_lnc = FALSE` – ubogie organizmy niemodelowe):**
  * Gdy liczba znanych `lncRNA` w bazie jest zbyt mała na podział na dwa etapy, nie wyodrębniamy podklasy `lnc`.
  * Zarówno w `cp`, jak i w `ml` klasą pozytywną jest cały zbiór `nc` przeciwko `cds`.

---

## 3. Teoria Statystyki a Dobór Liczebności Próbek

Dobór progów liczbowych i limitów (`max_n_cp`, `max_n_ml`) opiera się na trzech zasadach statystyki i uczenia maszynowego:

1. **Reguła EPV (*Events Per Variable* – reguła Peduzziego / Harrella):**
   * Na każdy predyktor (cechę wejściową w modelu `ml`) powinno przypadać minimum **15–20 obserwacji z klasy mniejszościowej** w zbiorze treningowym (`Train`).
   * *Przykład:* Przy 15 cechach wejściowych w etapie `ml` absolutne minimum wynosi $15 \times 20 = 300$ próbek `lncRNA` w zbiorze `ml_train` (czyli ok. 375 próbek w całym bloku `ml` przed podziałem 80/20, a ~750 próbek łącznie w bazie przed podziałem 50/50 na `cp` i `ml`). Poniżej tej wartości gwałtownie rośnie wariancja modelu i ryzyko przeuczenia (*overfitting*).
2. **Wypłaszczenie krzywej uczenia (*Power Law Plateau*):**
   * Dla klasyfikatorów tabelarycznych (Random Forest, XGBoost, LightGBM) operujących na kilkunastu/kilkudziesięciu cechach błąd generalizacji spada wykładniczo tylko do poziomu ok. **1 500 – 2 500 próbek na klasę** w zbiorze treningowym. Powyżej tego progu krzywa uczenia ulega wypłaszczeniu (przyrost AUC/F1 < 0.2%), a rośnie jedynie koszt obliczeniowy.
3. **Szerokość 95% przedziałów ufności (CI) w zbiorze Testowym:**
   * Przy 100 próbkach testowych błąd standardowy dla skuteczności 90% wynosi $\pm 6.0\%$ (95% CI: 84.0%–96.0%).
   * Przy **400–500 próbkach testowych** błąd spada do $\pm 2.6\%$ (95% CI: 87.4%–92.6%), zapewniając wysoką precyzję raportowanych metryk.

### Tabela Decyzyjna (Wybór Ścieżki i Limitów)

| Liczba dostępnych `lncRNA` po deduplikacji | Rekomendowana ścieżka | Sugerowane limity (`max_n_cp`, `max_n_ml`) | Przykładowy bilans: Mysz (~5k `lnc`, ~22k `cds`, ~40k `ex`) | Przykładowy bilans: Gatunek niemodelowy (~300 `lnc`, ~15k `cds`, ~50k `ex`) |
| :--- | :--- | :--- | :--- | :--- |
| **< 600** | **Ścieżka B** (`extract_lnc = FALSE`) | `Inf` (lub limit 2500 jeśli całe `nc` jest bardzo duże) | *Nie dotyczy* | **Całe `nc` (np. 1800) + `cds`:**<br>• **CP:** 900 poz / 900 neg<br>• **ML:** 900 poz / 900 neg<br>• **EX:** ~50k (predykcja) |
| **600 – 1 500** | **Ścieżka A** (`extract_lnc = TRUE`) | `max_n_cp = Inf`<br>`max_n_ml = Inf` | *Nie dotyczy* | **Wszystkie `lnc` (np. 1000) + `cds`:**<br>• **CP:** 500 `nc` / 500 `cds`<br>• **ML:** 500 `lnc` / 500 `cds`<br>• **EX:** ~50k (predykcja) |
| **> 1 500 – 5 000+** | **Ścieżka A** (`extract_lnc = TRUE`) | `max_n_cp = 2500`<br>`max_n_ml = 2500` | **Przycięte do limitów:**<br>• **CP:** 2500 `nc` / 2500 `cds`<br>• **ML:** 2500 `lnc` / 2500 `cds`<br>• **EX:** ~40k (predykcja) | *Nie dotyczy* |

---

## 4. Mapa Przepływu Danych przez Skrypt R

Poniższa tabela odwzorowuje dokładne działanie funkcji `load_and_dedup()` oraz `split_datasets()` dla **Ścieżki A (`extract_lnc = TRUE`)** przy przykładowych parametrach: `cp_stage_prop = 0.4`, `cp_tt_prop = 0.8`, `ml_tt_prop = 0.75`, `max_n_cp = 2000`, `max_n_ml = 3000`, `ml_neg_mode = "cds_only"`.

| Plik wejściowy (`load_and_dedup`) | `raw_group` | Priorytet deduplikacji (`dedup_mode`) | Przypisana `class` (`extract_lnc = TRUE`) | Podział etapów (`cp_stage_prop = 0.4`) | Etap 1: `stage == "cp"` (`max_n_cp = 2000`, `cp_tt_prop = 0.8`) | Etap 2: `stage == "ml"` (`max_n_ml = 3000`, `ml_tt_prop = 0.75`) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`nc_file`** *(lub `lnc_file`)*<br>Pasuje do `lnc_pattern` | `"nc"` *(lub `"lnc_ext"`)* | Priorytet 2 / 3 | **`"lnc"`**<br>*(np. 5 000)* | • **40% -> `"cp"`** (2 000)<br>• **60% -> `"ml"`** (3 000) | **`role = "pos"`** *(razem z `nc_other`)*:<br>Pula przycięta do `max_n_cp = 2000`:<br>• **`cp_train` (`keep=TRUE`):** 1 600<br>• **`cp_test` (`keep=TRUE`):** 400 | **`role = "pos"`**:<br>Pula 3 000 mieści się w `max_n_ml = 3000`:<br>• **`ml_train` (`keep=TRUE`):** 2 250 (75%)<br>• **`ml_test` (`keep=TRUE`):** 750 (25%) |
| **`nc_file`**<br>Brak wzorca `lnc_pattern` | `"nc"` | Priorytet 3 | **`"nc_other"`**<br>*(np. 25 000)* | • **40% -> `"cp"`** (10 000)<br>• **60% -> `"ml"`** (15 000) | **`role = "pos"`** *(wchodzi do wspólnej puli `"pos"` w `cp`)*:<br>Część powyższych 2 000 próbek (`keep=TRUE`) | **Gdy `ml_neg_mode = "cds_only"`:**<br>• `role = NA`, **`keep = FALSE`** *(odrzucone)*<br>**Gdy `"cds_and_nc_other"`:**<br>• `role = "neg"` *(razem z `cds`)* |
| **`cds_file`** | `"cds"` | **Priorytet 1**<br>*(wzorzec nadrzędny)* | **`"cds"`**<br>*(np. 22 000)* | • **40% -> `"cp"`** (8 800)<br>• **60% -> `"ml"`** (13 200) | **`role = "neg"`**:<br>Przycięte do `max_n_cp = 2000`:<br>• **`cp_train` (`keep=TRUE`):** 1 600<br>• **`cp_test` (`keep=TRUE`):** 400 | **`role = "neg"`**:<br>Przycięte do `max_n_ml = 3000`:<br>• **`ml_train` (`keep=TRUE`):** 2 250 (75%)<br>• **`ml_test` (`keep=TRUE`):** 750 (25%) |
| **`ex_file`** *(StringTie)* | `"ex"` | Priorytet 4<br>*(usuwa duplikaty obecne w `cds`, `lnc`, `nc`)* | **`"ex"`**<br>*(np. 40 000)* | **`stage = "inference"`**<br>*(100% do predykcji)* | *Nieużywany w etapie `cp`* | **`role = "unknown"`**, `tt_set = "predict"`:<br>• **`ex_pred` (`keep=TRUE`):** 40 000 |

### Porównanie ról (`role`) w obu ścieżkach

| Etap (`stage`) | Rola w modelu (`role`) | Ścieżka A (`extract_lnc = TRUE`) | Ścieżka B (`extract_lnc = FALSE`) |
| :--- | :--- | :--- | :--- |
| **`stage == "cp"`** | **`"pos"`** | `class %in% c("lnc", "nc_other")` | `class == "nc"` |
| **`stage == "cp"`** | **`"neg"`** | `class == "cds"` | `class == "cds"` |
| **`stage == "ml"`** | **`"pos"`** | `class == "lnc"` | `class == "nc"` |
| **`stage == "ml"`** | **`"neg"`** | `class == "cds"` *(lub `c("cds", "nc_other")`)* | `class == "cds"` |
| **`stage == "inference"`** | **`"unknown"`** | `class == "ex"` | `class == "ex"` |

### Wyciąganie gotowych obiektów (`extract_subset()`)

| Obiekt w R | Wywołanie funkcji | Warunki filtrowania w `wynik$legend` | Przeznaczenie |
| :--- | :--- | :--- | :--- |
| **`cp_train`** | `extract_subset(wynik, "cp", "train")` | `stage == "cp" & tt_set == "train" & keep == TRUE` | Trening narzędzi `cp` |
| **`cp_test`** | `extract_subset(wynik, "cp", "test")` | `stage == "cp" & tt_set == "test" & keep == TRUE` | Walidacja i dobór progów `cp` |
| **`ml_train`** | `extract_subset(wynik, "ml", "train")` | `stage == "ml" & tt_set == "train" & keep == TRUE` | Trening meta-klasyfikatora `ml` |
| **`ml_test`** | `extract_subset(wynik, "ml", "test")` | `stage == "ml" & tt_set == "test" & keep == TRUE` | Końcowa ewaluacja pipeline'u `cp` + `ml` |
| **`ex_pred`** | `extract_subset(wynik, "inference")` | `stage == "inference" & keep == TRUE` | Predykcja nowych `lncRNA` w danych RNA-seq |

---

## 5. Implementacja w R (`dataset_preparation.R`)

Wymagany pakiet: `Biostrings` (Bioconductor).

```r
library(Biostrings)

# FUNKCJA 1: Wczytanie plików FASTA, tagowanie lncRNA i opcjonalna deduplikacja
load_and_dedup <- function(cds_file,
                           nc_file,
                           lnc_file    = NULL,
                           ex_file     = NULL,
                           extract_lnc = TRUE,
                           lnc_pattern = "lncRNA|lnc_rna|long_non_coding",
                           dedup_mode  = "exact_seq") {

  read_fasta_df <- function(filepath, grp) {
    if (is.null(filepath)) return(NULL)
    seqs <- readDNAStringSet(filepath)
    data.frame(
      old_name  = names(seqs),
      raw_group = grp,
      seq_char  = as.character(seqs),
      stringsAsFactors = FALSE
    )
  }

  # Priorytet deduplikacji: cds > lnc_ext > nc > ex
  df <- rbind(
    read_fasta_df(cds_file, "cds"),
    read_fasta_df(lnc_file, "lnc_ext"),
    read_fasta_df(nc_file,  "nc"),
    read_fasta_df(ex_file,  "ex")
  )

  # Opcjonalna deduplikacja
  if (dedup_mode == "exact_seq") {
    df <- df[!duplicated(df$seq_char), ]
  } else if (dedup_mode == "header_id") {
    clean_ids <- sub("\\s+.*$", "", df$old_name)
    df <- df[!duplicated(clean_ids), ]
  } else if (dedup_mode != "none") {
    stop("Nieznany dedup_mode. Wybierz: 'none', 'exact_seq' lub 'header_id'.")
  }

  # Przypisanie klas: Ścieżka A (extract_lnc = TRUE) vs Ścieżka B (FALSE)
  is_lnc <- (df$raw_group == "lnc_ext") |
            (df$raw_group == "nc" & grepl(lnc_pattern, df$old_name, ignore.case = TRUE))

  df$class <- df$raw_group
  if (extract_lnc) {
    df$class[is_lnc] <- "lnc"
    df$class[df$raw_group == "nc" & !is_lnc] <- "nc_other"
  } else {
    df$class[df$raw_group %in% c("nc", "lnc_ext")] <- "nc"
  }

  return(df)
}

# FUNKCJA 2: Szybka inspekcja liczebności przed podziałem
inspect_counts <- function(df) {
  cat("\n--- LICZEBNOŚĆ KLAS PO WCZYTANIU I DEDUPLIKACJI ---\n")
  print(table(Klasa = df$class))
  cat("---------------------------------------------------\n\n")
}

# FUNKCJA 3: Podział na etapy (cp / ml), zbiory (train / test) oraz limity
split_datasets <- function(df,
                           prefix        = "Species",
                           seed          = 123,
                           cp_stage_prop = 0.5,
                           cp_tt_prop    = 0.8,
                           ml_tt_prop    = 0.8,
                           max_n_cp      = Inf,
                           max_n_ml      = Inf,
                           balance_train = TRUE,
                           ml_neg_mode   = "cds_only") {
  set.seed(seed)

  # Nadanie unikalnych ID
  df$ID <- ave(df$class, df$class, FUN = function(x) {
    paste(prefix, x, sprintf("%07d", seq_along(x)), sep = "_")
  })

  df_ref <- df[df$class != "ex", ]
  df_exp <- df[df$class == "ex", ]

  split_exact <- function(n, labels, prop) {
    if (n == 0) return(character(0))
    n1 <- round(n * prop)
    sample(c(rep(labels[1], n1), rep(labels[2], n - n1)))
  }

  # 1. Podział na etapy: cp vs ml (stratyfikowany po 'class')
  df_ref$stage <- ave(df_ref$class, df_ref$class, FUN = function(x) {
    split_exact(length(x), c("cp", "ml"), cp_stage_prop)
  })

  # 2. Podział train vs test z osobnymi proporcjami dla cp i ml
  df_ref$tt_set <- NA_character_
  for (st in c("cp", "ml")) {
    prop <- if (st == "cp") cp_tt_prop else ml_tt_prop
    idx  <- which(df_ref$stage == st)
    df_ref$tt_set[idx] <- ave(df_ref$class[idx], df_ref$class[idx], FUN = function(x) {
      split_exact(length(x), c("train", "test"), prop)
    })
  }

  df_ref$cohort <- paste(df_ref$class, df_ref$stage, df_ref$tt_set, sep = "_")
  df_ref$keep   <- TRUE

  # 3. Przypisanie ról (pos / neg)
  df_ref$role <- NA_character_
  has_lnc <- "lnc" %in% df_ref$class

  if (has_lnc) {
    # Ścieżka A
    df_ref$role[df_ref$stage == "cp" & df_ref$class %in% c("lnc", "nc_other")] <- "pos"
    df_ref$role[df_ref$stage == "cp" & df_ref$class == "cds"]                    <- "neg"
    df_ref$role[df_ref$stage == "ml" & df_ref$class == "lnc"]                    <- "pos"
    if (ml_neg_mode == "cds_only") {
      df_ref$role[df_ref$stage == "ml" & df_ref$class == "cds"]                  <- "neg"
      df_ref$keep[df_ref$stage == "ml" & df_ref$class == "nc_other"]             <- FALSE
    } else {
      df_ref$role[df_ref$stage == "ml" & df_ref$class %in% c("cds", "nc_other")] <- "neg"
    }
  } else {
    # Ścieżka B
    df_ref$role[df_ref$class == "nc"]  <- "pos"
    df_ref$role[df_ref$class == "cds"] <- "neg"
  }

  # 4. Osobne limity (max_n_cp, max_n_ml) i balansowanie 1:1
  for (st in c("cp", "ml")) {
    lim_total <- if (st == "cp") max_n_cp else max_n_ml
    tt_p      <- if (st == "cp") cp_tt_prop else ml_tt_prop

    for (tt in c("train", "test")) {
      lim <- if (tt == "train") floor(lim_total * tt_p) else floor(lim_total * (1 - tt_p))
      idx_pos <- which(df_ref$stage == st & df_ref$tt_set == tt & df_ref$role == "pos" & df_ref$keep)
      idx_neg <- which(df_ref$stage == st & df_ref$tt_set == tt & df_ref$role == "neg" & df_ref$keep)

      if (tt == "train" && balance_train) {
        target_n <- min(length(idx_pos), length(idx_neg), lim)
      } else {
        target_n <- min(max(length(idx_pos), length(idx_neg)), lim)
      }

      if (length(idx_pos) > target_n) {
        df_ref$keep[idx_pos[sample.int(length(idx_pos), length(idx_pos) - target_n)]] <- FALSE
      }
      if (length(idx_neg) > target_n) {
        df_ref$keep[idx_neg[sample.int(length(idx_neg), length(idx_neg) - target_n)]] <- FALSE
      }
    }
  }

  # 5. Dołączenie zbioru ex (jeśli istnieje)
  if (nrow(df_exp) > 0) {
    df_exp$stage  <- "inference"
    df_exp$tt_set <- "predict"
    df_exp$cohort <- "ex_inference"
    df_exp$keep   <- TRUE
    df_exp$role   <- "unknown"
    df_final <- rbind(df_ref, df_exp)
  } else {
    df_final <- df_ref
  }

  all_seqs <- DNAStringSet(df_final$seq_char)
  names(all_seqs) <- df_final$ID
  df_final$seq_char <- NULL

  return(list(seqs = all_seqs, legend = df_final))
}

# FUNKCJA 4: Pomocnicze wyciąganie gotowych podzbiorów sekwencji (FASTA)
extract_subset <- function(wynik, stage_name, tt_name = NULL, role_name = NULL, only_kept = TRUE) {
  df <- wynik$legend
  mask <- df$stage == stage_name
  if (!is.null(tt_name))   mask <- mask & (df$tt_set == tt_name)
  if (!is.null(role_name)) mask <- mask & (df$role == role_name)
  if (only_kept)           mask <- mask & df$keep
  
  return(wynik$seqs[df$ID[mask]])
}
```

---

## 6. Przykład Użycia i Testy Spójności

```r
# 1. Wczytanie i deduplikacja (sprawdzamy Ścieżkę A: extract_lnc = TRUE)
dane_surowe <- load_and_dedup(
  cds_file    = "Mus_musculus.GRCm39.cds.all.fa.gz",
  nc_file     = "Mus_musculus.GRCm39.ncrna.fa.gz",
  ex_file     = "merged_transcriptome.fasta",
  extract_lnc = TRUE,
  lnc_pattern = "lncRNA|lnc_rna",
  dedup_mode  = "exact_seq"
)

# 2. Podgląd liczebności klas w konsoli
inspect_counts(dane_surowe)

# Decyzja: jeśli 'lnc' < 600 -> przełączamy na Ścieżkę B (extract_lnc = FALSE)

# 3. Podział na zbiory z własnymi proporcjami i limitami
wynik <- split_datasets(
  df            = dane_surowe,
  prefix        = "Mus",
  seed          = 123,
  cp_stage_prop = 0.4,        # 40% referencji do CP, 60% do ML
  cp_tt_prop    = 0.8,        # 80/20 w etapie CP
  ml_tt_prop    = 0.75,       # 75/25 w etapie ML
  max_n_cp      = 2000,       # max 2000 na klasę w CP
  max_n_ml      = 3000,       # max 3000 na klasę w ML
  balance_train = TRUE,
  ml_neg_mode   = "cds_only"
)

# 4. Ekstrakcja gotowych obiektów DNAStringSet
cp_train <- extract_subset(wynik, "cp", "train")
cp_test  <- extract_subset(wynik, "cp", "test")
ml_train <- extract_subset(wynik, "ml", "train")
ml_test  <- extract_subset(wynik, "ml", "test")
ex_pred  <- extract_subset(wynik, "inference")

# 5. Kontrola spójności i braku wycieku danych (wszystkie przecięcia muszą zwrócić 0)
df_keep <- wynik$legend[wynik$legend$keep, ]
table(Etap_Zbior = paste(df_keep$stage, df_keep$tt_set, sep = "_"), Rola = df_keep$role)

stopifnot(
  length(intersect(names(cp_train), names(ml_train))) == 0,
  length(intersect(as.character(cp_train), as.character(ml_test))) == 0,
  length(intersect(as.character(ml_train), as.character(ex_pred))) == 0
)
```
