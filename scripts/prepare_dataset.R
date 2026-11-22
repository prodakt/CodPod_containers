library(Biostrings)

# ==============================================================================
# FUNKCJA 1: Wczytanie plików FASTA, tagowanie lncRNA i opcjonalna deduplikacja
# ==============================================================================
load_and_dedup <- function(cds_file,
                           nc_file,
                           lnc_file    = NULL,
                           ex_file     = NULL,
                           extract_lnc = TRUE,
                           lnc_pattern = "lncRNA|lnc_rna|long_non_coding",
                           dedup_mode  = "exact_seq") { # "none", "exact_seq", "header_id"

  read_to_df <- function(file, group_name) {
    if (is.null(file)) return(NULL)
    s <- readDNAStringSet(file)
    data.frame(old_name = names(s), raw_group = group_name, seq_char = as.character(s), stringsAsFactors = FALSE)
  }

  # Kolejność priorytetów przy deduplikacji: cds > lnc_ext > nc > ex
  df <- rbind(
    read_to_df(cds_file, "cds"),
    read_to_df(lnc_file, "lnc_ext"),
    read_to_df(nc_file,  "nc"),
    read_to_df(ex_file,  "ex")
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


# ==============================================================================
# FUNKCJA 2: Szybka inspekcja liczebności przed podziałem
# ==============================================================================
inspect_counts <- function(df) {
  cat("\n--- LICZEBNOŚĆ KLAS PO WCZYTANIU I DEDUPLIKACJI ---\n")
  print(table(Klasa = df$class))
  cat("---------------------------------------------------\n\n")
}


# ==============================================================================
# FUNKCJA 3: Podział na etapy (cp / ml), zbiory (train / test) oraz limity
# ==============================================================================
split_datasets <- function(df,
                           prefix        = "Species",
                           seed          = 123,
                           cp_stage_prop = 0.5,        # Jaka część referencji trafia do etapu 'cp' (reszta do 'ml')
                           cp_tt_prop    = 0.8,        # Proporcja train/test w etapie 'cp'
                           ml_tt_prop    = 0.8,        # Proporcja train/test w etapie 'ml'
                           max_n_cp      = Inf,        # Limit próbek na rolę (pos/neg) w etapie 'cp' (train + test)
                           max_n_ml      = Inf,        # Limit próbek na rolę (pos/neg) w etapie 'ml' (train + test)
                           balance_train = TRUE,       # Wyrównanie 1:1 w zbiorach train
                           ml_neg_mode   = "cds_only") # W Ścieżce A: "cds_only" lub "cds_and_nc_other"
{
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


# ==============================================================================
# FUNKCJA 4: Pomocnicze wyciąganie gotowych podzbiorów sekwencji (FASTA)
# ==============================================================================
extract_subset <- function(wynik, stage_name, tt_name = NULL, role_name = NULL, only_kept = TRUE) {
  df <- wynik$legend
  mask <- df$stage == stage_name
  if (!is.null(tt_name))   mask <- mask & (df$tt_set == tt_name)
  if (!is.null(role_name)) mask <- mask & (df$role == role_name)
  if (only_kept)           mask <- mask & df$keep
  
  return(wynik$seqs[df$ID[mask]])
}
