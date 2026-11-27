# ///////////////////////////////////////////////////////////////////////////////
# prepare_dataset_v4.1.R | CodPod_containers / lncRNA | wersja robocza
# Charakterystyka transkryptow i przygotowanie zbiorow CP -> ML -> EX
# Dokumentacja metodologiczna: prepare_dataset_v2.md (sekcje 00-16).
# V4 (robocza): minimalna modyfikacja V3 uzytkownika, zachowane sekcje 00-16,
# definicje cech, klasy biologiczne, funkcje CP/ML/EX i dotychczasowe eksporty.
# V4.1: poprawka importu REF; jeden lub wiele FASTA + wspolny GTF.
# Test Mm115: cDNA 116358, ncRNA 158776, brak wspolnych ID, razem 275134.
# Wszystkie 152133 lncRNA z GTF sa w ncRNA FASTA, zadne w cDNA FASTA.
# Nie zmieniamy biotypow, cech, algorytmow splicingu ani splitu CP/ML/EX.
# Cel 1: Ensembl REF-only (pelne transkrypty FASTA + GTF, bez EX).
# Cel 2: Ensembl REF + EX (pelne cDNA/GTF + eksperymentalne FASTA/GTF).
# W obu przypadkach genom DNA jest OPCJONALNY i niezalezny od wyboru REF/EX.
# V3 REF-only przetestowano na rzeczywistych danych; V4 wymaga ponownego testu.
# Windows: scanFa() w Rsamtools 2.28.0 zglaszal bledy i przerwal sesje R
# przy duzym genomie GRCm39. Dokladnej przyczyny nie ustalono. Biostrings
# odczytal ten sam region zgodnie z samtools na Linuksie (test uzytkownika).
# V4: biostrings / rsamtools / none. Domyslnie biostrings, bez indeksowania.
# Import genomu myszy: 55.33 s, obiekt ok. 2.5 GiB (pomiar uzytkownika).
# Nie jest to pomiar szczytowej pamieci calego przebiegu REF + EX.
# Testy syntetyczne 00A: wykonaj przed uruchomieniem duzych danych.
# Nie wymaga Rcpp ani kompilatora na typowej instalacji binarnej pakietow.
# Uruchamiaj kolejne sekcje w RStudio; ustaw sciezki w sekcji 00.
# Nie instaluje pakietow, nie pobiera danych i nie zapisuje nic na GitHubie.
# Przyklad korzysta z pelnego cDNA, nie z wycietych sekwencji CDS.
# ///////////////////////////////////////////////////////////////////////////////

# 00. PAKIETY I KONFIGURACJA ----------------------------------------------------
# Instalacja jednorazowa (odkomentuj w razie potrzeby):
# install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "rtracklayer", "GenomicRanges",
#                        "GenomicFeatures", "txdbmaker", "Rsamtools"))
# install.packages("ggplot2")

# Ustaw katalog roboczy samodzielnie (RStudio: Session > Set Working Directory).
# Przyklad uzytkownika; odkomentuj tylko na odpowiednim komputerze:
# setwd("D:/lncRna/PrepDataset/")

# Ustawienia wejscia: "ensembl" (jeden lub wiele FASTA pelnych transkryptow i GTF), "separate" (oddzielne pliki PC, NC, opcjonalnie lncRNA i GTF) lub "experimental_only" (transkryptom eksperymentalny i opcjonalna adnotacja referencyjna).
input_mode <- "ensembl"

# Dwa testowane scenariusze korzystaja z tego samego input_mode="ensembl".
# Roznia sie TYLKO dolaczeniem EX, nie sposobem klasyfikacji referencji.
#   1. REF-only: include_ex=FALSE. Cztery kohorty CP/ML, ex_pred pusta.
#   2. REF+EX:   include_ex=TRUE. Te same role REF, dodatkowo ex_pred.
# W obu scenariuszach run_split=TRUE, jezeli przygotowujemy dane do CP i ML.
# Uwaga: split na poziomie transkryptu, nie genu; gene_overlap.tsv to audyt.
# Nie oznacza to niezaleznej walidacji na innym organizmie lub wydaniu bazy.
#
# Przyklad nowych plikow Ensembl geneset (REF-only, sprawdzony w V3):
# ZASTAP aktywne wartosci ponizej tym przykladem (nie tylko odkomentuj).
# ref_fasta <- "cdna.fa.bgz"
# ref_gtf <- "genes.gtf.bgz"
# ref_id_mode <- "ensembl"
# reference_release <- "Ensembl geneset 2026_04"
# assembly <- "GRCm39"                 # sprawdz zgodnosc z wybranymi plikami
# Przyklad aktywny: Ensembl release 115, REF + EX (test V4.1).
# W klasycznym Ensembl cDNA i ncRNA sa osobnymi plikami. Sam cdna.all
# nie wystarcza do treningu lncRNA. Nie uzywamy cds.all (wyciete CDS).
# ref_fasta: pojedyncza sciezka LUB wektor c(...) w dowolnej kolejnosci.
# Wszystkie pliki musza pochodzic z tego samego wydania i assembly co ref_gtf.
# Etykiety PC/NC nadal wyznacza transcript_biotype w GTF, NIE nazwa pliku.

setwd("D:/lncRna/PrepDataset/")

ref_fasta <- c("Mm115/Mus_musculus.GRCm39.cdna.all.fa.gz",
               "Mm115/Mus_musculus.GRCm39.ncrna.fa.gz")
ref_gtf <- "Mm115/Mus_musculus.GRCm39.115.gtf.gz"
ref_id_mode <- "ensembl"
reference_release <- "Ensembl release 115"

assembly <- "GRCm39"
prefix <- "Mus"

# Alternatywa: osobne FASTA. Najlepiej pelne transkrypty kodujace, NIE CDS.
pc_file <- "protein_coding.fa.gz"
nc_file <- "non_coding.fa.gz"
lnc_file <- NULL                     # opcjonalny, moze byc podzbiorem nc_file
ambiguous_file <- NULL               # opcjonalny
separate_gtf <- NULL                 # np. "genes.gtf.gz"; NULL = cechy FASTA
separate_representation <- "full_transcript" # lub "CDS_only" (nie do pelnego ML)

# JEDEN przelacznik dla scenariuszy 1 i 2. Nie nadpisuj go drugi raz.
# FALSE: tylko dane z bazy, ex_file i ex_gtf nie sa otwierane ani sprawdzane.
# TRUE:  dane z bazy + dane eksperymentalne (aktualny test uzytkownika).
# EX nigdy nie dostaje etykiet PC/NC z GTF; sluzy do pozniejszej predykcji.
# Dokladnie identyczne sekwencje EX i REF sa wykluczane z ex_pred (sekcja 12).
include_ex <- TRUE
ex_file <- "Mm115/super_transcriptome.fa"
ex_gtf <- "Mm115/super_merged.gtf"

ex_id_mode <- "exact"                # NIE usuwamy kropek z MSTRG.123.1
# Przy experimental_only mozna podac ref_gtf jako niezalezna adnotacje PC.
# Dla braku takiej referencji ustaw ref_gtf <- NULL.

# GTF bezposrednio = zachowuje metadane i obsluguje GTF bez wierszy transcript.
# "txdb" = alternatywa zgodna z poprzednim kodem, kontroluje zgodnosc eksonow.
structure_engine <- "gtf"            # "gtf" lub "txdb"
# V4: wybor silnika analizy DNA genomowego, NIE wybor zrodla transkryptow.
# "biostrings": oryginalny FASTA/FASTA.gz; bez FAI/GZI, genom raz w RAM.
#   Domyslnie dla Windows. Wczytanie 2.5 GiB genomu myszy + pamiec tymczasowa
#   i pozostale obiekty REF/EX wymagaja odpowiedniego zapasu RAM.
# "rsamtools": lokalny FASTA z .fai lub BGZF z .fai i .gzi; odczyt z dysku.
#   Mniejsze zapotrzebowanie na RAM, ale zaleznosc od indeksow i bibliotek.
#   Na Windows przetestuj scanFa() osobno; tryCatch nie chroni przed awaria R.
# "none": nie otwieramy genomu DNA; cechy splicingu = NA, nie 0.
#   Nadal liczymy cechy eksonowo-intronowe i kontekst PC z GTF (sekcje 05-06).
# NULL w genome_fasta jest zgodnym z V3 skrotem wylaczajacym analize DNA.
# Zmiana silnika nie zmienia definicji cech ani klasyfikacji biologicznej.
# Nie przelaczamy automatycznie silnika po bledzie, aby nie ukryc problemu.
genome_engine <- "biostrings"         # "biostrings", "rsamtools", "none"
genome_fasta <- "Mm115/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
# Dla rsamtools podaj np. ...fa (nieskompresowany) lub ...fa.bgz (BGZF).
# Dla REF-only bez genomu wystarczy genome_engine <- "none".
# W obydwu scenariuszach mozna uzywac dowolnego z trzech silnikow.
# FASTA genomu i GTF musza miec zgodny assembly i nazwy chromosomow.
# Nie dopasowujemy automatycznie 1/chr1, MT/chrM ani roznych wersji assembly.
splice_chunk_size <- 100000L          # maks. liczba intronow w porcji, NIE limit analizy
# Genom jest usuwany z pamieci po sekcji 07; pozostaja cechy transkryptow.
# Zmniejsz splice_chunk_size, jezeli porcje zuzywaja zbyt duzo RAM.
# Ten parametr nie zmniejsza rozmiaru samego genomu w RAM (biostrings).
calc_kmers <- FALSE                   # 16 dinukleotydow zawsze opisowo; TRUE dodaje 64 trimery
pas_window <- 100L                    # ostatnie 100 nt; motywy w sekwencji DNA
min_orf_nt <- 6L                     # ATG + STOP, wliczamy STOP

# V3: kontrola zasobow. Limity dotycza porcji, NIE liczby analizowanych sekwencji.
# Zmniejsz oba limity przy niedoborze RAM; zwieksz po pomiarze czasu.
sequence_chunk_size <- 3000L         # maksymalna liczba transkryptow w porcji
sequence_chunk_nt <- 8000000         # limit sumy nt (wyjatek: pojedynczy dlugi transkrypt)
progress_enabled <- TRUE             # komunikaty postepu w konsoli RStudio
save_sequence_checkpoint <- TRUE     # RDS po zakonczeniu kazdego zrodla (gdy write_outputs)
# Nie wczytujemy automatycznie checkpointow: zmienione FASTA/parametry moglyby
# spowodowac uzycie nieaktualnych cech. Patrz komentarz w sekcji 04.

# Ustawienia do pozniejszego podzialu (sekcje 12-14).
run_split <- TRUE                    # FALSE: tylko charakterystyka i eksport, brak podzialu na train/test
extract_lnc <- TRUE                  # TRUE = A, FALSE = B; decyzja reczna!
dedup_mode <- "exact_seq"            # "exact_seq", "header_id", "none"
exclude_length_mismatch <- TRUE       # wyklucz niezgodne FASTA/GTF z treningu ref
# Proporcje stosujemy osobno do kazdej klasy biologicznej, przed limitowaniem.
# Seed umozliwia powtorzenie losowania na tych samych uporzadkowanych danych.
seed <- 123L
cp_stage_prop <- 0.4                  # 40% REF -> CP, pozostale 60% -> ML
cp_tt_prop <- 0.8                     # CP: 80% train, 20% test
ml_tt_prop <- 0.75                    # ML: 75% train, 25% test
# Limity NA ROLE (pos/neg osobno), lacznie train+test, a nie liczba wszystkich
# transkryptow. Inf = bez tego limitu. Balansowanie moze zmniejszyc train.
# Wartosci 10000/20000 zachowane z V3 uzytkownika: robocze, nie rekomendacja
# statystyczna do publikacji. Dobor liczebnosci wymaga pozniejszej walidacji.
max_n_cp <- 10000
max_n_ml <- 20000
balance_train <- TRUE                # rowna liczba pos/neg tylko w train
# A (extract_lnc=TRUE): CP: cale NC vs PC; ML: lnc vs PC (cds_only),
# albo lnc vs PC + pozostale NC (cds_and_nc_other). B: NC vs PC w obu etapach.
ml_neg_mode <- "cds_only"
# EPV jest jedynie heurystyka, nie gwarancja mocy. Uwzglednij przyszle wyniki CP.
planned_ml_parameters <- 25L
epv_values <- c(10, 15, 20)
test_target_sensitivity <- 0.90
test_ci_half_width <- 0.03

# Rozne katalogi dla REF-only i REF+EX oraz roznych silnikow. Przy zmianie
# wydania bazy lub powtorzeniu eksperymentu ustaw nowy lub pusty katalog.
# Przyklad dla nowego geneset: out_dir <- "Mus_2026_04_REF_none_v4.1"
out_dir <- paste0("Mm115_analysis_",if(include_ex) "REF_EX" else "REF", "_",
                  if(is.null(genome_fasta)) "none" else genome_engine,"_v4.1")
write_outputs <- TRUE
make_plots <- TRUE                   # ggplot2; wykresy opisowe
run_self_tests <- FALSE              # TRUE uruchamia male testy przed importem
self_tests_only <- FALSE             # TRUE: tylko testy; source() zatrzyma sie po nich
# Opcjonalny test zgodnosci dwoch silnikow na MALYM, syntetycznym FASTA.
# Dziala tylko przy run_self_tests=TRUE lub self_tests_only=TRUE. Nie dotyka
# genomu uzytkownika. Domyslnie FALSE, aby nie uruchamiac scanFa na Windows.
self_test_rsamtools <- FALSE
# Zachowanie calej sesji R jest opcjonalne, moze zajac wiele GB na dysku.
# TRUE zapisuje RData w out_dir po sekcji 16, tylko gdy write_outputs=TRUE.
save_workspace <- FALSE

# Podstawowe pakiety; stosujemy jawne namespace, by uniknac konfliktu select().
required <- c("Biostrings", "rtracklayer", "GenomicRanges", "IRanges", "S4Vectors")
if (!self_tests_only && structure_engine == "txdb") required <- c(required, "GenomicFeatures", "txdbmaker")
stopifnot(is.logical(include_ex),length(include_ex)==1L,!is.na(include_ex),
          length(genome_engine)==1L,genome_engine %in% c("biostrings","rsamtools","none"),
          length(splice_chunk_size)==1L,is.finite(splice_chunk_size),
          splice_chunk_size>=1,splice_chunk_size==floor(splice_chunk_size),
          splice_chunk_size<=.Machine$integer.max)
# Efektywny silnik zapisujemy w run_config.R. NULL ma pierwszenstwo jak w V3.
if (is.null(genome_fasta)) genome_engine <- "none"
if ((!self_tests_only && genome_engine=="rsamtools") ||
    ((run_self_tests || self_tests_only) && self_test_rsamtools))
  required <- c(required, "Rsamtools")
if (!self_tests_only && make_plots) required <- c(required, "ggplot2")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Brak pakietow: ", paste(missing_packages, collapse = ", "))
# S4 metody width/start/end/strand itd. dostepne po zaladowaniu pakietow.
library(Biostrings)
library(GenomicRanges)
# Nie dolaczamy Rsamtools do search path. Uzywamy Rsamtools:: w jego silniku.
# rtracklayer moze nadal importowac Rsamtools jako zaleznosc pakietowa.
stopifnot(input_mode %in% c("ensembl", "separate", "experimental_only"),
          structure_engine %in% c("gtf", "txdb"),
          separate_representation %in% c("full_transcript", "CDS_only"))
if (input_mode == "separate" && separate_representation == "CDS_only" && run_split)
  stop("CDS_only: ustaw run_split=FALSE lub dostarcz pelne transkrypty kodujace.")
if (input_mode == "experimental_only") run_split <- FALSE
if (self_tests_only) { run_self_tests <- TRUE; write_outputs <- FALSE }
# Prosta kontrola sciezki przed kosztownym importem. Odczyt DNA dopiero w 07.
if (!self_tests_only && genome_engine!="none") {
  if (!is.character(genome_fasta) || length(genome_fasta)!=1L ||
      is.na(genome_fasta) || !file.exists(genome_fasta) || dir.exists(genome_fasta))
    stop("Brak lokalnego FASTA genomu: ",genome_fasta,". Ustaw poprawna sciezke lub genome_engine='none'.")
  if (genome_engine=="rsamtools" && !file.exists(paste0(genome_fasta,".fai")))
    stop("Brak FAI genomu. Uzyj indeksowanego FASTA/BGZF albo silnika biostrings.")
}
stopifnot(sequence_chunk_size >= 1, sequence_chunk_size == floor(sequence_chunk_size),
          is.finite(sequence_chunk_nt), sequence_chunk_nt >= 1,
          pas_window >= 1, pas_window == floor(pas_window), min_orf_nt >= 6,
          min_orf_nt == floor(min_orf_nt))
if (write_outputs) {
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    stop("Katalog wynikowy nie jest pusty. Wybierz nowy out_dir (ochrona przed nadpisaniem).")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Słowniki biotypow: operacyjne klasy projektu, nie nowa ontologia.
# Nie uzywamy negacji protein_coding jako definicji non_coding.
pc_biotypes <- c("protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
                 "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene")
lncrna_biotypes <- c("lncRNA", "lincRNA", "antisense", "sense_intronic", "sense_overlapping",
                     "bidirectional_promoter_lncRNA", "3prime_overlapping_ncRNA", "macro_lncRNA",
                     "long_non_coding_RNA")
short_ncrna_biotypes <- c("miRNA", "snRNA", "snoRNA", "rRNA", "sRNA", "scaRNA", "scRNA",
                         "misc_RNA", "ribozyme", "vault_RNA", "Y_RNA", "Mt_rRNA", "Mt_tRNA", "tRNA")
all_ncrna_biotypes <- unique(c(lncrna_biotypes, short_ncrna_biotypes))

# ///////////////////////////////////////////////////////////////////////////////
# MALE FUNKCJE POMOCNICZE ---- (definicje przed kodem glownym 01-16)
# ///////////////////////////////////////////////////////////////////////////////

# 02: Pierwsze pole FASTA, a wersje usuwamy TYLKO z rozpoznanych ID Ensembl.
extract_transcript_ids <- function(x, mode = "exact") {
  mode <- match.arg(mode, c("exact", "ensembl"))
  ids <- sub("[[:space:]].*$", "", trimws(as.character(x)))
  if (mode == "ensembl") ids <- sub("^(ENS[A-Z]*T[0-9]+)\\.[0-9]+$", "\\1", ids)
  ids
}

# Lokalny plik lub URL. Zalecamy lokalny cache i przypieta wersje Ensembl.
read_fasta <- function(file) {
  if (is.null(file)) return(DNAStringSet())
  if (!grepl("^(https?|ftp)://", file) && !file.exists(file)) stop("Brak FASTA: ", file)
  s <- readDNAStringSet(file)
  if (!length(s) || is.null(names(s)) || any(!nzchar(names(s)))) stop("Pusty FASTA lub brak nazw: ", file)
  if (any(width(s) == 0L)) stop("Puste sekwencje: ", file)
  s
}

# V4.1 / 01-02: adapter wielu FASTA tylko dla input_mode="ensembl".
# read_fasta() pozostaje czytnikiem pojedynczego pliku, takze dla genomu.
# Zachowujemy kolejnosc plikow i rekordow; bez konwersji calych FASTA na tekst.
# Deduplikacja ID tutaj NIE jest deduplikacja biologiczna sekwencji z sekcji 12.
# Ten sam znormalizowany ID w dwoch plikach mozna scalic TYLKO gdy oryginalny
# ID (wlacznie z wersja) i cala sekwencja sa identyczne. Rozne naglowki opisowe
# sa dopuszczalne; zostaje pierwszy naglowek, a fasta_files zachowuje wszystkie
# sciezki pochodzenia (separator " | "). Konflikt lub duplikat w jednym pliku
# zatrzymuje analize. Rozne ID o tej samej sekwencji zachowujemy do sekcji 12.
read_reference_fastas <- function(files, mode="ensembl") {
  if (!is.character(files) || !length(files) || anyNA(files) ||
      any(!nzchar(files)) || anyDuplicated(files))
    stop("ref_fasta: podaj jedna lub kilka roznych, niepustych sciezek FASTA.")
  ss <- lapply(files,read_fasta)
  s <- if(length(ss)==1L) ss[[1L]] else do.call(c,unname(ss))
  origins <- rep(files,lengths(ss)); ids <- extract_transcript_ids(names(s),mode)
  if(anyNA(ids) || any(!nzchar(ids))) stop("Brak ID w referencyjnym FASTA.")
  for(j in seq_along(ss)) {
    if(anyDuplicated(extract_transcript_ids(names(ss[[j]]),mode)))
      stop("Nieunikalne ID po normalizacji wewnatrz FASTA: ",files[j])
  }
  repeated <- duplicated(ids) | duplicated(ids,fromLast=TRUE)
  if(any(repeated)) {
    groups <- split(which(repeated),ids[repeated])
    for(i in groups) {
      if(length(unique(extract_transcript_ids(names(s)[i],"exact")))!=1L)
        stop("Konflikt wersji ID pomiedzy FASTA: ",ids[i[1L]])
      if(length(unique(as.character(s[i])))!=1L)
        stop("Ten sam ID, rozne sekwencje pomiedzy FASTA: ",ids[i[1L]])
      origins[i[1L]] <- paste(unique(origins[i]),collapse=" | ")
    }
  }
  keep <- !duplicated(ids)
  message("REF FASTA: ",length(s)," rekordow, ",sum(keep)," unikalnych ID; scalono ",sum(!keep)," kopii.")
  list(seqs=s[keep],fasta_files=origins[keep])
}

# 03: Pobranie jednej kolumny GTF, bez bledu gdy atrybut nie wystepuje.
gtf_field <- function(gr, field) {
  if (!field %in% names(S4Vectors::mcols(gr))) return(rep(NA_character_, length(gr)))
  x <- as.character(S4Vectors::mcols(gr)[[field]])
  x[!nzchar(x)] <- NA_character_
  x
}

# 03: Metadane transkryptu. Wiersze transcript maja pierwszenstwo, exony sa
# alternatywa dla StringTie i GTF bez jawnych rekordow transcript.
# Konflikty biotypow lub gene_id dla tego samego ID zatrzymuja analize.
get_tx_meta <- function(gr, mode) {
  types <- gtf_field(gr, "type")
  q <- gr[which(types %in% c("transcript", "mRNA", "exon"))]
  raw <- gtf_field(q, "transcript_id")
  ok <- !is.na(raw) & nzchar(raw)
  q <- q[ok]; raw <- raw[ok]
  if (!length(q)) stop("GTF nie ma transcript_id dla transcript/exon. GFF3 wymaga mapowania ID/Parent.")
  ids <- extract_transcript_ids(raw, mode)
  # V3: zamiast vapply() dla kazdego ID i kazdego pola, redukujemy cala
  # kolumne do unikalnych par ID-wartosc. Brak wartosci nie jest konfliktem.
  # Kolejnosc alfabetyczna ID jak w v2 (split); zachowujemy wszystkie metadane.
  keys <- sort(unique(ids)); owner <- match(ids, keys)
  pick <- function(x, field) {
    i <- which(!is.na(x) & nzchar(x))
    pairs <- unique(data.frame(owner=owner[i], value=x[i], stringsAsFactors=FALSE))
    dup <- anyDuplicated(pairs$owner)
    if (dup) stop("Konflikt ", field, " dla ", keys[pairs$owner[dup]])
    ans <- rep(NA_character_, length(keys))
    ans[pairs$owner] <- pairs$value
    ans
  }
  out <- data.frame(transcript_id=keys, gtf_transcript_id=pick(raw,"transcript_id"), stringsAsFactors=FALSE)
  for (f in c("gene_id","gene_name","gene_biotype","transcript_biotype","transcript_version")) {
    v <- gtf_field(q,f)
    if (f %in% c("gene_biotype","transcript_biotype")) {
      missing <- is.na(v)
      if (any(missing)) v[missing] <- gtf_field(q,sub("biotype$","type",f))[missing]
    }
    out[[f]] <- pick(v,f)
  }
  out
}

# 03: Eksony bezposrednio z GTF. Usuwamy tylko powtorzone zakresy tego samego
# transkryptu, nie eksony wspolne dla roznych izoform.
get_exons_by_tx <- function(gr, mode, engine = "gtf") {
  if (engine == "txdb") {
    txdb <- txdbmaker::makeTxDbFromGRanges(gr)
    ex <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
    names(ex) <- extract_transcript_ids(names(ex),mode)
    if (anyNA(names(ex)) || anyDuplicated(names(ex))) stop("TxDb: nieunikalne ID")
    # Ta sama walidacja i kolejnosc co w adapterze GTF.
    keys <- names(ex); e <- unlist(ex,use.names=FALSE)
    owner <- rep(seq_along(ex),lengths(ex))
  } else {
    e <- gr[which(gtf_field(gr,"type")=="exon")]
    ids <- extract_transcript_ids(gtf_field(e,"transcript_id"),mode)
    ok <- !is.na(ids) & nzchar(ids)
    e <- e[ok]; ids <- ids[ok]
    if (!length(e)) stop("Brak eksonow z transcript_id w GTF.")
    keys <- sort(unique(ids)); owner <- match(ids,keys)
  }
  if (!length(e)) stop("Brak eksonow.")
  # V3: jedna operacja sortowania na plaskich wektorach, zadnego endoapply.
  # Ekson jest unikalny W TRANSKRYPCIE, nie w calym genomie.
  S4Vectors::mcols(e) <- NULL
  chr <- as.character(seqnames(e)); st <- as.character(strand(e))
  pairs <- unique(data.frame(owner=owner,chr=chr,st=st))
  if (anyDuplicated(pairs$owner)) stop("Transkrypt obejmuje rozne chromosomy lub nici.")
  keep <- !duplicated(data.frame(owner=owner,start=start(e),end=end(e)))
  e <- e[keep]; owner <- owner[keep]; st <- st[keep]
  # Uklad 5' -> 3' dla obu nici. Gwiazdka oznacza nieznana orientacje.
  o <- order(owner,ifelse(st=="-",-start(e),start(e)),end(e),method="radix")
  counts <- tabulate(owner,nbins=length(keys))
  if (any(counts==0L)) stop("Transkrypt bez eksonow w adapterze.")
  relist(e[o],IRanges::PartitioningByEnd(cumsum(counts),names=keys))
}

# V4.1 / 03A: pokrycie GTF przez FASTA wedlug biotypu transkryptu.
# Liczymy unikalne ID po normalizacji, nie geny i nie wiersze eksonow GTF.
# Unia ID: GTF-only, FASTA-only oraz matched; brak biotypu = osobna grupa.
# Pokrycie = matched / gtf_n. NA gdy brak rekordow GTF w danej grupie.
# Nie jest to jeszcze kontrola wersji, dlugosci, etykiet ani deduplikacji
# sekwencji: te pozostaja w 08, 09 i 12. Dla EX biotyp pozostaje opisowy.
calc_fasta_coverage <- function(f, m, source) {
  gtf_ids <- if(is.null(m)) character() else m$transcript_id
  ids <- union(f$transcript_id,gtf_ids)
  has_fasta <- ids %in% f$transcript_id; has_gtf <- ids %in% gtf_ids
  b <- if(is.null(m)) rep(NA_character_,length(ids)) else m$transcript_biotype[match(ids,gtf_ids)]
  b[is.na(b) | !nzchar(b)] <- "(missing_biotype)"
  status <- ifelse(!has_gtf,"fasta_only",ifelse(!has_fasta,"gtf_only","matched"))
  t <- table(b,factor(status,levels=c("matched","gtf_only","fasta_only")))
  out <- data.frame(source=rep(source,nrow(t)),transcript_biotype=rownames(t),
    matched_n=as.integer(t[,"matched"]),gtf_only_n=as.integer(t[,"gtf_only"]),
    fasta_only_n=as.integer(t[,"fasta_only"]),stringsAsFactors=FALSE)
  out$gtf_n <- out$matched_n+out$gtf_only_n
  out$fasta_n <- out$matched_n+out$fasta_only_n
  out$coverage_percent <- ifelse(out$gtf_n>0,100*out$matched_n/out$gtf_n,NA_real_)
  out
}

# V4.1: wczesny warunek konieczny dla splitu w trybie ensembl. Nie zmieniamy
# extract_lnc ani slownikow, nie uzupelniamy klas na podstawie nazwy FASTA.
# Dodatnie liczebnosci nie gwarantuja poprawnego splitu po QC/deduplikacji.
check_reference_classes <- function(coverage,extract_lnc) {
  pc <- sum(coverage$matched_n[coverage$transcript_biotype %in% pc_biotypes])
  nc <- sum(coverage$matched_n[coverage$transcript_biotype %in% all_ncrna_biotypes])
  lnc <- sum(coverage$matched_n[coverage$transcript_biotype %in% lncrna_biotypes])
  cat("03A: REF FASTA+GTF | PC=",pc," | NC=",nc," | lncRNA=",lnc,"\n",sep="")
  if(pc==0L || (if(extract_lnc) lnc==0L else nc==0L))
    stop("Brak wymaganej klasy REF przed sekcja 04. Sprawdz ref_fasta i ref_gtf; klasyczny Ensembl wymaga zwykle cdna.all + ncrna. Nie zmieniaj extract_lnc tylko w celu ominiecia bledu.",call.=FALSE)
  invisible(c(pc=pc,nc=nc,lnc=lnc))
}

# 04: Najdluzszy kompletny ORF: ATG do PIERWSZEGO STOP w tej samej ramce.
# Wszystkie starty liczone osobno (ORF zagniezdzone tez); 3 ramki na nici cDNA.
# Standardowy kod jadrowy. Nie jest to dowod translacji; nie liczymy ORF bez STOP.
calc_longest_orf_bp <- function(dna_set, min_nt = 6L) {
  n <- length(dna_set)
  out <- data.frame(max_orf_len=integer(n),orf_count=integer(n),
                    orf_start=rep(NA_integer_,n),orf_end=rep(NA_integer_,n))
  if (!n) return(out)
  # Dopasowania wzorcow wykonuje Biostrings. Wektor pozycji + wlasciciel
  # pozwalaja uniknac petli po kazdym transkrypcie i kazdej ramce odczytu.
  hits <- function(pattern) {
    x <- start(vmatchPattern(pattern,dna_set))
    list(pos=unlist(x,use.names=FALSE),owner=rep(seq_len(n),lengths(x)))
  }
  a <- hits("ATG"); if (!length(a$pos)) return(out)
  b <- lapply(c("TAA","TAG","TGA"),hits)
  sp <- unlist(lapply(b,`[[`,"pos"),use.names=FALSE)
  so <- unlist(lapply(b,`[[`,"owner"),use.names=FALSE)
  if (!length(sp)) return(out)
  # Kodujemy (transkrypt, ramka, pozycja) w rozlacznych przedzialach liczbowych.
  # Double zapobiega przepelnieniu 32-bit; w typowej porcji klucze sa dokladne.
  stride <- as.double(width(dna_set))+1
  offset <- c(0,cumsum(3*stride))[seq_len(n)]
  if (sum(3*stride)>=2^53) stop("Za duza porcja dla kluczy ORF.")
  ak <- offset[a$owner]+(a$pos%%3L)*stride[a$owner]+a$pos
  bk <- offset[so]+(sp%%3L)*stride[so]+sp
  o <- order(bk,method="radix"); bk <- bk[o]; sp <- sp[o]; so <- so[o]
  # findInterval + 1 = pierwszy STOP scisle za START; sprawdzamy wlasciciela
  # i ramke, zeby nie przejsc przypadkiem do nastepnej ramki/transkryptu.
  j <- findInterval(ak,bk)+1L
  i <- which(j<=length(bk)); j <- j[i]
  good <- so[j]==a$owner[i] & sp[j]%%3L==a$pos[i]%%3L
  i <- i[good]; j <- j[good]
  len <- sp[j]-a$pos[i]+3L
  good <- len>=min_nt; i <- i[good]; len <- len[good]
  if (!length(i)) return(out)
  owner <- a$owner[i]; beg <- a$pos[i]
  out$orf_count <- tabulate(owner,nbins=n)
  # Remis dlugosci: najwczesniejszy START, identycznie jak w v2.
  o <- order(owner,-len,beg,method="radix"); best <- o[!duplicated(owner[o])]
  k <- owner[best]
  out$max_orf_len[k] <- len[best]; out$orf_start[k] <- beg[best]
  out$orf_end[k] <- beg[best]+len[best]-1L
  out
}

# 04: Cechy FASTA obliczane niezaleznie od etykiety. Sklad jako ulamek 0-1.
# GC = (G+C)/(A+C+G+T); niejednoznaczne symbole raportujemy osobno.
calc_sequence_features <- function(s, ids, pas_window=100L, min_orf_nt=6L, kmers=FALSE) {
  stopifnot(length(s)==length(ids),all(width(s)>0))
  n <- width(s); acgt <- letterFrequency(s,letters=c("A","C","G","T","N"))
  denom <- rowSums(acgt[,1:4,drop=FALSE]); p <- acgt[,1:4,drop=FALSE]/ifelse(denom==0,NA_real_,denom)
  ent <- -rowSums(ifelse(p>0,p*log2(p),0)); ent[denom==0] <- NA_real_
  out <- data.frame(transcript_id=ids,seq_length=n,gc_content=p[,"G"]+p[,"C"],
    a_content=p[,"A"],c_content=p[,"C"],g_content=p[,"G"],t_content=p[,"T"],
    n_fraction=acgt[,"N"]/n,ambiguous_fraction=1-denom/n,sequence_entropy=ent)
  out <- cbind(out,calc_longest_orf_bp(s,min_orf_nt)); out$orf_coverage <- out$max_orf_len/n
  # PAS: wyszukiwanie obu motywow w Biostrings, rowniez nakladajacych sie.
  tail_seq <- subseq(s,start=pmax(1L,n-pas_window+1L)); tail_n <- width(tail_seq)
  last_hit <- function(pattern) {
    x <- start(vmatchPattern(pattern,tail_seq)); counts <- lengths(x)
    ans <- rep(-Inf,length(s)); i <- which(counts>0L)
    ans[i] <- unlist(x,use.names=FALSE)[cumsum(counts)[i]]
    ans
  }
  last <- pmax(last_hit("AATAAA"),last_hit("ATTAAA"))
  out$pas_present <- is.finite(last)
  out$pas_distance_from_end <- ifelse(out$pas_present,tail_n-last-5L,NA_real_)
  # Odwracamy sekwencje bez komplementowania; kotwica ^A+ nie skanuje
  # calego tekstu w poszukiwaniu wewnetrznych ciagow A. Nie jest to dowod poly(A).
  tail_a <- regexpr("^A+",as.character(reverse(s)),perl=TRUE)
  out$terminal_A_run <- pmax(0L,attr(tail_a,"match.length"))
  for (k in if(kmers)c(2L,3L) else 2L) {
    freq <- oligonucleotideFrequency(s,width=k,as.prob=FALSE)
    den <- rowSums(freq); freq <- freq/ifelse(den==0,NA_real_,den)
    colnames(freq) <- paste0("k",k,"_",colnames(freq))
    out <- cbind(out,as.data.frame(freq))
  }
  out
}

# 05: Statystyki eksonow i intronow z rzeczywistych zakresow.
# Brak intronow: liczba i suma = 0, srednia/min/max/SD/ratio = NA (nie dotyczy).
# Nakladajace sie eksony w jednym transkrypcie sa bledem struktury, nie sumujemy
# ich po cichu. Nieznana nic: nie interpretujemy pierwszego/ostatniego eksonu.
# V3: sumy, mediany, min/max i SD bez osobnych obiektow GRanges dla kazdego ID.
# Sortujemy wartosci w grupach; mediana = srodkowe indeksy w plaskim wektorze.
# Dla pustych grup suma=0, reszta NA. SD probkowe, tak jak stats::sd().
grouped_width_stats <- function(w, owner, n) {
  counts <- tabulate(owner,nbins=n); active <- which(counts>0L)
  total <- numeric(n); squares <- numeric(n)
  if (length(w)) {
    sums <- rowsum(cbind(w=as.double(w),sq=as.double(w)^2),owner,reorder=TRUE)
    idx <- as.integer(rownames(sums)); total[idx] <- sums[,1]; squares[idx] <- sums[,2]
  }
  o <- order(owner,w,method="radix"); sorted <- w[o]
  hi <- cumsum(counts); lo <- hi-counts+1L
  mn <- mx <- med <- avg <- dev <- rep(NA_real_,n)
  mn[active] <- sorted[lo[active]]; mx[active] <- sorted[hi[active]]
  med[active] <- (as.double(sorted[floor((lo[active]+hi[active])/2)])+
                  as.double(sorted[ceiling((lo[active]+hi[active])/2)]))/2
  avg[active] <- total[active]/counts[active]
  j <- which(counts>1L)
  dev[j] <- sqrt(pmax(0,(squares[j]-total[j]^2/counts[j])/(counts[j]-1L)))
  list(count=counts,total=total,mean=avg,median=med,min=mn,max=mx,sd=dev)
}

# Wspolny plaski uklad rosnaco po wspolrzednych genomowych.
# Uwaga: pierwszego eksonu biologicznego nie mylimy z lewym na nici minus.
flat_exons <- function(ex) {
  n <- lengths(ex); z <- unlist(ex,use.names=FALSE)
  owner <- rep(seq_along(ex),n)
  o <- order(owner,start(z),end(z),method="radix")
  z <- z[o]; owner <- owner[o]
  hi <- cumsum(n); lo <- hi-n+1L
  i <- seq_len(max(0L,length(z)-1L)); i <- i[owner[i]==owner[i+1L]]
  if (any(start(z)[i+1L]-end(z)[i]-1L<=0L)) stop("Nakladajace lub przylegajace eksony w transkrypcie.")
  list(z=z,owner=owner,n=n,lo=lo,hi=hi,pairs=i)
}

calc_exon_features <- function(ex, introns=NULL) {
  # Jezeli introny zostaly juz wyznaczone w 03, nie wyznaczamy ich drugi raz.
  f <- flat_exons(ex); z <- f$z; w <- width(z); n <- length(ex)
  a <- grouped_width_stats(w,f$owner,n)
  if (is.null(introns)) {
    i <- f$pairs; gaps <- start(z)[i+1L]-end(z)[i]-1L
    b <- grouped_width_stats(gaps,f$owner[i],n)
  } else {
    stopifnot(identical(names(ex),names(introns)),all(lengths(introns)==pmax(f$n-1L,0L)))
    b <- grouped_width_stats(width(unlist(introns,use.names=FALSE)),rep(seq_len(n),lengths(introns)),n)
  }
  st <- as.character(strand(z))[f$lo]
  first <- ifelse(st=="-",w[f$hi],w[f$lo]); last <- ifelse(st=="-",w[f$lo],w[f$hi])
  first[st=="*"] <- NA; last[st=="*"] <- NA
  span <- end(z)[f$hi]-start(z)[f$lo]+1L
  data.frame(transcript_id=names(ex),exon_count=a$count,tx_length=a$total,
    mean_exon_len=a$mean,median_exon_len=a$median,min_exon_len=a$min,max_exon_len=a$max,sd_exon_len=a$sd,
    first_exon_len=first,last_exon_len=last,genomic_span=span,exon_density=a$total/span,
    first_exon_fraction=first/a$total,last_exon_fraction=last/a$total,
    intron_count=b$count,total_intron_len=b$total,mean_intron_len=b$mean,median_intron_len=b$median,
    min_intron_len=b$min,max_intron_len=b$max,sd_intron_len=b$sd,
    exon_intron_ratio=ifelse(b$total>0,a$total/b$total,NA_real_),row.names=NULL)
}

# 05: Intron jako przerwa miedzy kolejnymi eksonami tej samej izoformy.
get_introns_by_tx <- function(ex) {
  f <- flat_exons(ex); z <- f$z; i <- f$pairs
  a <- end(z)[i]+1L; b <- start(z)[i+1L]-1L
  # Jeden GRanges, potem relist; jednoeksonowe transkrypty = puste grupy.
  intr <- GRanges(seqnames=seqnames(z)[i],IRanges::IRanges(a,b),strand=strand(z)[i],seqinfo=seqinfo(z))
  relist(intr,IRanges::PartitioningByEnd(cumsum(pmax(f$n-1L,0L)),names=names(ex)))
}

# 05: Zakresy transkryptow od pierwszej do ostatniej pozycji genomowej.
get_tx_ranges <- function(ex) {
  f <- flat_exons(ex); z <- f$z
  if (any(f$n==0L)) stop("Transkrypt bez eksonow.")
  ans <- GRanges(seqnames(z)[f$lo],IRanges::IRanges(start(z)[f$lo],end(z)[f$hi]),
                 strand=strand(z)[f$lo],seqinfo=seqinfo(z))
  names(ans) <- names(ex); ans
}

# 06: Słownik referencyjnych loci PC. W pierwszej kolejnosci gene_biotype;
# gdy brak gene_biotype, jedynie jawny coding transcript jako zrodlo genu PC.
# Do zakresu genu wlaczamy wszystkie jego znane izoformy. TSS odnosimy do
# startow poszczegolnych transkryptow tych genow, a nie jednego startu genu.
build_pc_reference <- function(meta, ex, tx, introns) {
  m <- meta[match(names(tx), meta$transcript_id), , drop=FALSE]
  gene <- m$gene_id
  pc <- !is.na(gene) & (m$gene_biotype %in% pc_biotypes |
                       (is.na(m$gene_biotype) & m$transcript_biotype %in% pc_biotypes))
  pc_genes <- unique(gene[pc]); idx <- which(gene %in% pc_genes)
  if (!length(idx)) return(NULL)
  # Zakres genu: min start i max end wszystkich izoform, bez petli po genach.
  z <- tx[idx]; keys <- sort(unique(gene[idx])); own <- match(gene[idx],keys)
  pairs <- unique(data.frame(own=own,chr=as.character(seqnames(z)),st=as.character(strand(z))))
  if (anyDuplicated(pairs$own)) stop("Gen PC na wielu kontigach/ niciach; sprawdz adnotacje.")
  lo <- order(own,start(z),method="radix"); lo <- lo[!duplicated(own[lo])]
  hi <- order(own,-end(z),method="radix"); hi <- hi[!duplicated(own[hi])]
  gene_ranges <- GRanges(seqnames(z)[lo],IRanges::IRanges(start(z)[lo],end(z)[hi]),strand=strand(z)[lo],seqinfo=seqinfo(z))
  names(gene_ranges) <- keys
  e <- ex[names(tx)[idx]]; intr <- introns[names(tx)[idx]]
  # Unia eksonow per nic redukuje liczbe redundantnych trafien izoform.
  list(genes=gene_ranges, tx=tx[idx], exons=reduce(unlist(e,use.names=FALSE),ignore.strand=FALSE),
       introns=unique(unlist(intr,use.names=FALSE)),
       contigs=unique(as.character(seqnames(tx))))
}

# 06: Obliczanie kontekstu genomowego. UWAGA: opisowe, nie Core ML.
# Dystans do zakresu genu = liczba baz pomiedzy zakresami; 0 takze dla sasiednich.
# Nakladanie sprawdzamy findOverlaps, nigdy przez sam warunek distance==0.
# Flagi obu nici moga byc jednoczesnie TRUE. Nieznana nic = NA dla orientacji.
calc_genomic_features <- function(tx, ex, meta, pc_ref, source="reference") {
  n <- length(tx); ids <- names(tx)
  out <- data.frame(transcript_id=ids, dist_to_pc=rep(NA_real_,n), nearest_pc_gene_id=rep(NA_character_,n),
    dist_to_nearest_tss=rep(NA_real_,n), genomic_class=rep(NA_character_,n), neighbor_orientation=rep(NA_character_,n),
    overlap_same_strand=rep(NA,n), overlap_opposite_strand=rep(NA,n), overlap_pc=rep(NA,n),
    exonic_overlap_bp=rep(NA_real_,n), exonic_overlap_fraction=rep(NA_real_,n), overlapping_pc_genes=rep(NA_integer_,n))
  if (is.null(pc_ref)) return(out)
  pc <- pc_ref$genes
  covered <- as.character(seqnames(tx)) %in% pc_ref$contigs
  known <- as.character(strand(tx)) %in% c("+","-")
  idx <- which(covered); if (!length(idx)) return(out)
  # Niezalezne flagi i liczba trafionych genow, bez wnioskowania z odleglosci.
  h <- findOverlaps(tx, pc, ignore.strand=TRUE)
  q <- S4Vectors::queryHits(h); s <- S4Vectors::subjectHits(h)
  qs <- as.character(strand(tx))[q]; ss <- as.character(strand(pc))[s]
  good <- qs %in% c("+","-") & ss %in% c("+","-")
  out$overlap_pc[idx] <- FALSE; out$overlap_pc[unique(q)] <- TRUE
  out$overlapping_pc_genes[idx] <- tabulate(q, n)[idx]
  out$overlap_same_strand[covered & known] <- FALSE
  out$overlap_opposite_strand[covered & known] <- FALSE
  out$overlap_same_strand[unique(q[good & qs==ss])] <- TRUE
  out$overlap_opposite_strand[unique(q[good & qs!=ss])] <- TRUE
  # Nieznana nic genu PC nie jest dowodem zgodnosci/ przeciwnosci nici.
  uncertain <- unique(q[!good]); out$overlap_same_strand[uncertain] <- NA
  out$overlap_opposite_strand[uncertain] <- NA
  near <- distanceToNearest(tx, pc, ignore.strand=TRUE)
  nq <- S4Vectors::queryHits(near); ns <- S4Vectors::subjectHits(near)
  out$dist_to_pc[nq] <- S4Vectors::mcols(near)$distance
  out$nearest_pc_gene_id[nq] <- names(pc)[ns]
  # TSS jako pojedyncza pozycja 5'; dystans = roznica wspolrzednych (nie gap).
  qt <- which(known); rt <- which(as.character(strand(pc_ref$tx)) %in% c("+","-"))
  if(length(qt) && length(rt)) {
    a <- promoters(tx[qt], upstream=0, downstream=1)
    b <- promoters(pc_ref$tx[rt], upstream=0, downstream=1)
    # distanceToNearest ma remis 0 przy overlap i sasiedztwie; dla TSS
    # obliczamy dokladna odleglosc punktow na kazdym chromosomie.
    ac <- as.character(seqnames(a)); bc <- as.character(seqnames(b))
    for(chr in intersect(unique(ac),unique(bc))) {
      i <- which(ac==chr); av <- start(a)[i]; bv <- sort(unique(start(b)[bc==chr]))
      j <- findInterval(av,bv)
      left <- ifelse(j>0L,abs(av-bv[pmax(j,1L)]),Inf)
      right <- ifelse(j<length(bv),abs(av-bv[pmin(j+1L,length(bv))]),Inf)
      out$dist_to_nearest_tss[qt[i]] <- pmin(left,right)
    }
  }
  # Orientacja wzgledem najblizszego genu; dla nakladajacych sie pozostaw opis overlap.
  a <- as.character(strand(tx))[nq]; b <- as.character(strand(pc))[ns]
  ok <- a %in% c("+","-") & b %in% c("+","-") & !out$overlap_pc[nq]
  ori <- ifelse(a==b,"tandem", ifelse((start(tx)[nq]<start(pc)[ns] & a=="-") |
       (start(tx)[nq]>start(pc)[ns] & a=="+"),"divergent","convergent"))
  out$neighbor_orientation[nq[ok]] <- ori[ok]
  # Rzeczywiste eksony query vs eksony PC. Sumujemy UNIE przecięc, nie trafienia.
  qe <- unlist(ex[ids], use.names=FALSE); owner <- rep(seq_len(n), lengths(ex[ids]))
  eh <- findOverlaps(qe, pc_ref$exons, ignore.strand=TRUE)
  eq <- S4Vectors::queryHits(eh); es <- S4Vectors::subjectHits(eh)
  out$exonic_overlap_bp[idx] <- 0
  same_ex <- integer()
  if(length(eq)) {
    a <- as.character(strand(qe))[eq]; b <- as.character(strand(pc_ref$exons))[es]
    same_ex <- unique(owner[eq[a==b & a %in% c("+","-")]])
  }
  # Do liczby bp stosujemy UNIE eksonow PC obu nici. Dzieki temu kazda
  # baza jest liczona tylko raz, nawet przy nakladaniu genow na obu niciach.
  # Eksony query sa rozlaczne w obrebie transkryptu (kontrola w 03).
  pc_union <- reduce(pc_ref$exons,ignore.strand=TRUE)
  uh <- findOverlaps(qe,pc_union,ignore.strand=TRUE)
  uq <- S4Vectors::queryHits(uh); us <- S4Vectors::subjectHits(uh)
  if(length(uq)) {
    bp <- pmin(end(qe)[uq],end(pc_union)[us])-pmax(start(qe)[uq],start(pc_union)[us])+1L
    sums <- rowsum(matrix(as.double(bp),ncol=1L),owner[uq],reorder=FALSE)
    out$exonic_overlap_bp[as.integer(rownames(sums))] <- sums[,1]
  }
  txlen <- as.numeric(sum(width(ex[ids])))
  out$exonic_overlap_fraction <- out$exonic_overlap_bp/txlen
  # Sense intronic: CALY zakres query miesci sie w konkretnym intronie PC.
  ih <- findOverlaps(tx,pc_ref$introns,type="within",ignore.strand=FALSE)
  intronic <- unique(S4Vectors::queryHits(ih))
  out$genomic_class[covered & !out$overlap_pc] <- "intergenic"
  out$genomic_class[covered & out$overlap_pc] <- "overlap_unresolved"
  opp <- which(out$overlap_opposite_strand %in% TRUE)
  out$genomic_class[opp] <- "antisense"
  same <- which(out$overlap_same_strand %in% TRUE)
  out$genomic_class[same] <- "sense_other_overlap"
  out$genomic_class[intersect(same,intronic)] <- "sense_intronic"
  out$genomic_class[intersect(same,same_ex)] <- "sense_overlapping"
  both <- intersect(same,opp); out$genomic_class[both] <- "mixed_strand_overlap"
  # Wlasny gen PC jest jawna informacja z adnotacji; NIGDY nie podajemy tej
  # kolumny do Core ML. Dla EX identyfikatory genow StringTie nie sa rownowazne.
  if(source=="reference") {
    g <- meta$gene_id[match(ids,meta$transcript_id)]
    out$genomic_class[which(g %in% names(pc))] <- "protein_coding_locus"
  }
  out
}

# 07: OPCJONALNA ANALIZA MOTYWOW SPLICINGOWYCH -------------------------------
# Wspolne definicje biologiczne V3/V4; rozni sie wylacznie dostep do DNA.
# Wspolrzedne GTF/GRanges sa 1-based inclusive. Intron to przerwa miedzy
# eksonami; pobieramy 2 nt na lewym i prawym koncu INTRONU (nie eksonu).
# Dla nici +: donor=lewy koniec, akceptor=prawy koniec.
# Dla nici -: donor=reverseComplement(prawy), akceptor=RC(lewy).
# GT-AG, GC-AG, AT-AC raportujemy oddzielnie; ich suma to operacyjna
# definicja canonical_junctions w projekcie. Nie przypisujemy typu spliceosomu.
# Inne pary A/C/G/T = noncanonical. Brak kontigu, pozycja poza genomem,
# nieznana nic, intron <4 nt lub zasady niejednoznaczne = unknown.
# To cechy sekwencji genomowej, nie dowod eksperymentalnego splicingu.
# Liczymy wystapienia intronow w transkryptach (wspolny intron izoform
# liczy sie wielokrotnie). Unikalne junctions i dalsze cechy to kolejny etap.

# Biostrings: jeden import bez indeksow. Zachowujemy tylko pierwsze pole
# naglowka, np. "7 dna:chromosome chromosome:GRCm39:7:..." -> "7".
# Nie konwertujemy calego genomu na character: to znacznie zwiekszyloby RAM.
read_genome_biostrings <- function(genome_file) {
  g <- read_fasta(genome_file)
  names(g) <- extract_transcript_ids(names(g),"exact")
  if(anyNA(names(g)) || any(!nzchar(names(g))) || anyDuplicated(names(g)))
    stop("Genom: brak lub powtorzone nazwy kontigow po obcieciu naglowkow.")
  g
}

# Wspolny szablon zachowuje wszystkie kolumny rowniez przy engine="none".
# NA = nie wykonano analizy; 0 = wykonano i stwierdzono brak zdarzen.
empty_splice_features <- function(introns) {
  n <- length(introns)
  data.frame(transcript_id=as.character(names(introns)),
    gt_ag=rep(NA_integer_,n),gc_ag=rep(NA_integer_,n),at_ac=rep(NA_integer_,n),
    canonical_junctions=rep(NA_integer_,n),noncanonical_junctions=rep(NA_integer_,n),
    unknown_junctions=rep(NA_integer_,n),noncanonical_fraction=rep(NA_real_,n))
}

# Wspolny etap interpretacji: oba silniki dostarczaja fragmenty w orientacji
# genomowej (+). Funkcja NIE otrzymuje etykiet PC/NC ani informacji o EX.
# Konwertujemy do DNAStringSet tylko krotkie fragmenty na nici minus.
classify_splice_motifs <- function(left,right,st) {
  stopifnot(length(left)==length(right),length(left)==length(st),
            all(st %in% c("+","-")))
  d <- as.character(left); ac <- as.character(right)
  minus <- st=="-"
  if(any(minus)) {
    d[minus] <- as.character(reverseComplement(DNAStringSet(right[minus])))
    ac[minus] <- as.character(reverseComplement(DNAStringSet(left[minus])))
  }
  good <- grepl("^[ACGT]{2}$",d) & grepl("^[ACGT]{2}$",ac)
  ans <- rep("unknown",length(st))
  ans[good] <- paste(d[good],ac[good],sep="-")
  ans
}

# Dla obu silnikow zachowujemy wspolne zasady poprawnosci zakresow oraz
# tabulacje per transkrypt. Nie przelaczamy automatycznie silnika po bledzie.
# Argument genome pozwala przekazac juz wczytany DNAStringSet (REF i EX
# korzystaja z jednego obiektu). Bez niego funkcja moze wykonac osobny import.
calc_splice_features <- function(introns, genome_file=NULL, engine="biostrings",
                                 genome=NULL, chunk_size=100000L) {
  engine <- match.arg(engine,c("biostrings","rsamtools","none"))
  stopifnot(length(chunk_size)==1L,is.finite(chunk_size),chunk_size>=1,
            chunk_size==floor(chunk_size),chunk_size<=.Machine$integer.max)
  n <- length(introns); counts <- lengths(introns)
  out <- empty_splice_features(introns)
  # NULL bez obiektu genomu zachowuje sie jak w V3 (brak analizy).
  if(engine=="none" || (is.null(genome_file) && is.null(genome))) return(out)
  if(engine=="biostrings") {
    if(is.null(genome)) genome <- read_genome_biostrings(genome_file)
    stopifnot(is(genome,"DNAStringSet"),length(genome)>0L,
              !is.null(names(genome)),!anyNA(names(genome)),
              all(nzchar(names(genome))),!anyDuplicated(names(genome)))
    lens <- setNames(width(genome),names(genome))
  } else {
    if(is.null(genome_file) || !file.exists(genome_file)) stop("Brak genomu dla rsamtools: ",genome_file)
    fai <- paste0(genome_file,".fai"); gzi <- paste0(genome_file,".gzi")
    if(!file.exists(fai)) stop("Brak FAI: ",fai)
    # Rozpoznajemy gzip po sygnaturze, nie rozszerzeniu. Sam zwykly gzip
    # nie jest BGZF; plik i indeksy musza byc przygotowane poza tym skryptem.
    gz <- identical(readBin(genome_file,what="raw",n=2L),as.raw(c(31,139)))
    if(gz && !file.exists(gzi)) stop("Kompresja: wymagany BGZF z GZI i FAI; zwykly gzip -> biostrings.")
    fa <- Rsamtools::FaFile(genome_file,index=fai,gzindex=gzi)
    idx <- Rsamtools::scanFaIndex(fa)
    lens <- setNames(width(idx),as.character(seqnames(idx)))
    if(!length(lens) || anyDuplicated(names(lens))) stop("Nieprawidlowy indeks genomu.")
  }
  z <- unlist(introns,use.names=FALSE); owner <- rep(seq_len(n),counts)
  motif <- rep("unknown",length(z)); chr <- as.character(seqnames(z)); st <- as.character(strand(z))
  # Kontrola nazw i dlugosci kontigow. Nie naprawiamy nazw ani wspolrzednych.
  chr_missing <- setdiff(unique(chr),names(lens))
  if(length(chr_missing)) warning("Brak kontigow w genomie (motywy unknown): ",
                                  paste(head(chr_missing,20L),collapse=", "),
                                  if(length(chr_missing)>20L) " ..." else "")
  limit <- unname(lens[chr]); in_bounds <- !is.na(limit) & start(z)>=1L & end(z)<=limit
  outside <- sum(!is.na(limit) & !in_bounds)
  if(outside>0L) warning("Introny poza zakresem kontigu: ",outside,"; sprawdz assembly.")
  valid <- which(st %in% c("+","-") & width(z)>=4L & in_bounds)
  if(length(z) && !length(valid)) warning("Brak intronow z poprawnymi zakresami i znana nicia; wszystkie motywy unknown.")
  if(length(valid)) {
    if(engine=="biostrings") {
      # Views dziala na pojedynczym chromosomie. Grupujemy indeksy intronow,
      # aby uniknac petli po milionach intronow. Views nie kopiuje chromosomu;
      # materializujemy jedynie krotkie fragmenty porcji. Wynik wpisujemy
      # do pozycji vi, wiec zachowana jest oryginalna kolejnosc transkryptow.
      groups <- split(valid,chr[valid]); processed <- 0
      for(chrom in names(groups)) {
        ii <- groups[[chrom]]; subject <- genome[[chrom]]
        for(j in seq.int(1L,length(ii),by=chunk_size)) {
          vi <- ii[j:min(length(ii),as.double(j)+chunk_size-1)]
          left <- as.character(Views(subject,start=start(z)[vi],width=2L))
          right <- as.character(Views(subject,start=end(z)[vi]-1L,width=2L))
          motif[vi] <- classify_splice_motifs(left,right,st[vi])
        }
        processed <- processed+length(ii)
        if(progress_enabled) { cat(sprintf("  %s | %d/%d intronow (%.1f%%)\n",chrom,processed,length(valid),100*processed/length(valid))); flush.console() }
      }
    } else {
      # Zachowana sciezka V3: indeksowany odczyt fragmentow z dysku.
      # GRanges bez podanej nici => odczyt genomowy; RC tylko w funkcji wspolnej.
      open(fa); on.exit(close(fa),add=TRUE)
      for(j in seq.int(1L,length(valid),by=chunk_size)) {
        vi <- valid[j:min(length(valid),as.double(j)+chunk_size-1)]
        v <- z[vi]
        left <- GRanges(as.character(seqnames(v)),IRanges::IRanges(start(v),start(v)+1L))
        right <- GRanges(as.character(seqnames(v)),IRanges::IRanges(end(v)-1L,end(v)))
        a <- Rsamtools::scanFa(fa,param=left); b <- Rsamtools::scanFa(fa,param=right)
        motif[vi] <- classify_splice_motifs(a,b,st[vi])
        if(progress_enabled) { cat(sprintf("  %d/%d intronow\n",min(length(valid),as.double(j)+chunk_size-1),length(valid))); flush.console() }
      }
    }
  }
  count <- function(m) tabulate(owner[m],nbins=n)
  out$gt_ag <- count(motif=="GT-AG"); out$gc_ag <- count(motif=="GC-AG"); out$at_ac <- count(motif=="AT-AC")
  out$canonical_junctions <- out$gt_ag+out$gc_ag+out$at_ac
  out$unknown_junctions <- count(motif=="unknown")
  out$noncanonical_junctions <- counts-out$canonical_junctions-out$unknown_junctions
  tested <- counts-out$unknown_junctions
  out$noncanonical_fraction <- ifelse(tested>0,out$noncanonical_junctions/tested,NA_real_)
  # Jednoeksonowy transkrypt: liczniki 0, fraction=NA (mianownik 0).
  # Dla transkryptu z samymi unknown fraction=NA, nie 0.
  stopifnot(all(out$canonical_junctions>=0L),all(out$noncanonical_junctions>=0L),
            all(out$unknown_junctions>=0L),
            all(out$canonical_junctions+out$noncanonical_junctions+out$unknown_junctions==counts))
  out
}

# Jeden import DNA dla wszystkich zrodel. Genom jest lokalny w funkcji;
# po zwroceniu tabel nie pozostaje w globalnym srodowisku. Nie usuwamy
# obiektu genome, ktory uzytkownik ewentualnie wczytal samodzielnie.
run_splice_analysis <- function(introns_by_tx,genome_file,engine,chunk_size) {
  engine <- match.arg(engine,c("biostrings","rsamtools","none"))
  if(is.null(genome_file)) engine <- "none"
  g <- NULL; ans <- list()
  if(!length(introns_by_tx)) return(ans)
  if(engine=="biostrings") {
    g <- timed_step("07: import genomu Biostrings",read_genome_biostrings(genome_file))
    cat("  kontigi: ",length(g)," | nt: ",sum(as.double(width(g))),
        " | rozmiar obiektu: ",format(object.size(g),units="auto"),"\n",sep="")
  }
  for(src in names(introns_by_tx))
    ans[[src]] <- timed_step(paste(src,"07: motywy splicingowe [",engine,"]",sep=" "),
      calc_splice_features(introns_by_tx[[src]],genome_file,engine,g,chunk_size))
  ans
}

# Kontrola wynikow per zrodlo, nie per unikalny junction. NA przy wylaczonej
# analizie, zeby sum(NA,na.rm=TRUE) nie zasugerowalo braku noncanonical.
# Mianownik procentu: tylko introny z rozpoznanym motywem (ACGT na obu koncach).
summarize_splice_features <- function(features,introns,source,engine) {
  intron_n <- sum(as.double(lengths(introns)))
  sm <- function(v) if(engine=="none") NA_real_ else sum(as.double(v))
  ans <- data.frame(source=source,genome_engine=engine,transcripts=length(introns),
    introns=intron_n,gt_ag=sm(features$gt_ag),gc_ag=sm(features$gc_ag),at_ac=sm(features$at_ac),
    canonical=sm(features$canonical_junctions),noncanonical=sm(features$noncanonical_junctions),
    unknown=sm(features$unknown_junctions))
  ans$resolved <- ans$canonical+ans$noncanonical
  ans$noncanonical_percent <- ifelse(ans$resolved>0,100*ans$noncanonical/ans$resolved,NA_real_)
  ans$resolved_percent <- ifelse(intron_n>0,100*ans$resolved/intron_n,NA_real_)
  if(engine!="none") stopifnot(ans$resolved+ans$unknown==intron_n)
  ans
}

# 08: Bezpieczne dolaczanie 1:1, bez mnozenia rekordow i zmiany kolejnosci.
join_by_tx <- function(x,y) {
  if(anyNA(y$transcript_id) || anyDuplicated(y$transcript_id)) stop("Nieunikalny klucz tabeli cech.")
  cols <- setdiff(names(y),"transcript_id")
  if(length(intersect(names(x),cols))) stop("Powtorzone nazwy cech: ",paste(intersect(names(x),cols),collapse=", "))
  cbind(x, y[match(x$transcript_id,y$transcript_id),cols,drop=FALSE], row.names=NULL)
}

# 10: Statystyki grup. Dane GTF i FASTA raportowane oddzielnie; statystyki
# sekwencji/struktury dla rekordow obecnych w FASTA (porownywalny mianownik).
calc_group_stats <- function(df, group_label, total_fasta) {
  x <- df[df$has_fasta %in% TRUE,,drop=FALSE]
  med <- function(z) if(any(!is.na(z))) median(z,na.rm=TRUE) else NA_real_
  avg <- function(z) if(any(!is.na(z))) mean(z,na.rm=TRUE) else NA_real_
  data.frame(group=group_label,n_gtf=sum(df$has_gtf),n_fasta=nrow(x),
    percent_fasta=if(total_fasta>0) 100*nrow(x)/total_fasta else NA_real_,
    n_genes=length(unique(na.omit(x$gene_id))),n_unique_sequences=length(unique(na.omit(x$sequence_group))),
    median_seq_length=med(x$seq_length),mean_seq_length=avg(x$seq_length),median_exons=med(x$exon_count),
    mean_exons=avg(x$exon_count),single_exon_percent=100*avg(x$exon_count==1),median_gc=med(x$gc_content),
    median_orf=med(x$max_orf_len),median_orf_coverage=med(x$orf_coverage),median_dist_pc=med(x$dist_to_pc),
    overlap_pc_percent=100*avg(x$overlap_pc),missing_structure_percent=if(nrow(x))100*mean(is.na(x$exon_count)) else NA_real_)
}

# 12: Dawna funkcja 1, rozszerzona o dane juz scharakteryzowane.
# Wariant z df: record_id, old_name, raw_group, seq_char (oraz opcjonalne metadane).
# raw_group: cds (historyczna nazwa roli PC), lnc_ext, nc, ex, ambiguous, annotation_missing.
# Deduplikacja NIE usuwa informacji z final_features. Zwraca audyt jako atrybut.
load_and_dedup <- function(cds_file=NULL,nc_file=NULL,lnc_file=NULL,ex_file=NULL,
                           extract_lnc=TRUE,lnc_pattern="lncRNA|lnc_rna|long_non_coding",
                           dedup_mode="exact_seq",df=NULL) {
  dedup_mode <- match.arg(dedup_mode,c("exact_seq","header_id","none"))
  if(is.null(df)) {
    read_to_df <- function(file,group) {
      if(is.null(file)) return(NULL)
      s <- read_fasta(file)
      data.frame(old_name=names(s),raw_group=group,seq_char=as.character(s))
    }
    df <- do.call(rbind,list(read_to_df(cds_file,"cds"),read_to_df(lnc_file,"lnc_ext"),
                            read_to_df(nc_file,"nc"),read_to_df(ex_file,"ex")))
    if(is.null(df)) stop("Nie podano sekwencji.")
    df$record_id <- paste0("legacy_",seq_len(nrow(df)))
  }
  stopifnot(all(c("record_id","old_name","raw_group","seq_char") %in% names(df)), !anyDuplicated(df$record_id))
  if(!all(df$raw_group %in% c("cds","lnc_ext","nc","ex","ambiguous","annotation_missing"))) stop("Nieznana raw_group")
  is_lnc <- if("is_lncrna" %in% names(df)) df$is_lncrna %in% TRUE else
    df$raw_group=="lnc_ext" | (df$raw_group=="nc" & grepl(lnc_pattern,df$old_name,ignore.case=TRUE))
  df$class <- df$raw_group
  df$class[df$raw_group %in% c("lnc_ext","nc")] <- if(extract_lnc) ifelse(is_lnc[df$raw_group %in% c("lnc_ext","nc")],"lnc","nc_other") else "nc"
  # Wszystkie referencje (tez ambiguous/missing) sluza do wykluczenia
  # identycznych EX, niezaleznie od tego czy zostana uzyte do uczenia.
  is_ex <- df$raw_group=="ex"
  reason <- rep("eligible",nrow(df))
  reason[is_ex & df$seq_char %in% df$seq_char[!is_ex]] <- "known_reference_exact"
  reason[df$raw_group %in% c("ambiguous","annotation_missing")] <- "not_labeled_for_training"
  # Identyczna sekwencja o sprzecznych etykietach PC vs NC: nie rozstrzygamy
  # biologii arbitralnie. Wykluczamy obie z treningu, zachowujac w audycie.
  pos <- df$class %in% c("lnc","nc_other","nc"); neg <- df$class=="cds"
  conflict <- intersect(df$seq_char[pos],df$seq_char[neg])
  reason[(pos|neg) & df$seq_char %in% conflict] <- "conflicting_labels"
  # W A identyczne lnc vs nc_other ma niejednoznaczna etykiete w ML;
  # zachowujemy dotychczasowy priorytet lnc, ale zapisujemy odrzucenie kopii.
  eligible <- which(reason=="eligible")
  ord <- eligible[order(match(df$raw_group[eligible],c("cds","lnc_ext","nc","ex")))]
  key <- if(dedup_mode=="header_id") extract_transcript_ids(df$old_name,"exact") else df$seq_char
  # Dla header_id uwzgledniamy source jesli dostepne: MSTRG i Ensembl to osobne przestrzenie.
  if(dedup_mode=="header_id" && "source" %in% names(df)) key <- paste(df$source,key,sep="::")
  representative <- rep(NA_character_,nrow(df)); representative[ord] <- df$record_id[ord]
  if(dedup_mode!="none") {
    first <- match(key[ord],key[ord]); dup <- duplicated(key[ord]); representative[ord] <- df$record_id[ord[first]]
    reason[ord[dup]] <- "duplicate"
  }
  audit <- data.frame(record_id=df$record_id,raw_group=df$raw_group,class=df$class,
                      representative_record_id=representative,reason=reason,stringsAsFactors=FALSE)
  out <- df[ord[reason[ord]=="eligible"],,drop=FALSE]
  attr(out,"dedup_audit") <- audit; attr(out,"extract_lnc") <- extract_lnc
  out
}

# 12: Dawna funkcja 2; szybka inspekcja przed reczna decyzja A/B.
inspect_counts <- function(df) {
  print(table(Klasa=df$class,useNA="ifany"))
  audit <- attr(df,"dedup_audit"); if(!is.null(audit)) print(table(audit$reason,audit$class))
  invisible(table(df$class))
}

# 13: Orientacyjny kalkulator liczebnosci. EPP dotyczy parametrow modelu,
# nie prostego liczenia kolumn. Zwraca spodziewana pule ograniczajacej klasy.
# Nie symuluje limitow ani rzeczywistych zaokraglen w split_datasets().
estimate_sample_size <- function(p,epv=c(10,15,20),cp_stage_prop=.4,ml_tt_prop=.75) {
  stopifnot(p>0,all(epv>0),cp_stage_prop>0,cp_stage_prop<1,ml_tt_prop>0,ml_tt_prop<1)
  target <- ceiling(p*epv)
  data.frame(parameters=p,epv=epv,target_ml_train_per_class=target,
             approximate_input_per_limiting_class=ceiling(target/((1-cp_stage_prop)*ml_tt_prop)))
}

# 14: Dawna funkcja 3. Ten sam schemat stratyfikowanego losowania, limitow
# i balansowania. Jawna sciezka A/B, kontrola parametrow i pustych kohort.
split_datasets <- function(df,prefix="Species",seed=123,cp_stage_prop=.5,cp_tt_prop=.8,ml_tt_prop=.8,
                           max_n_cp=Inf,max_n_ml=Inf,balance_train=TRUE,ml_neg_mode="cds_only",
                           extract_lnc=attr(df,"extract_lnc")) {
  ml_neg_mode <- match.arg(ml_neg_mode,c("cds_only","cds_and_nc_other"))
  props <- c(cp_stage_prop,cp_tt_prop,ml_tt_prop)
  stopifnot(all(is.finite(props)&props>0&props<1),all(c(max_n_cp,max_n_ml)>0),
            all(c(max_n_cp,max_n_ml)==floor(c(max_n_cp,max_n_ml))),
            is.logical(extract_lnc),length(extract_lnc)==1L,!is.na(extract_lnc),
            is.logical(balance_train),length(balance_train)==1L,!is.na(balance_train))
  expected <- if(extract_lnc) c("cds","lnc","nc_other","ex") else c("cds","nc","ex")
  if(!all(df$class %in% expected)) stop("Niezgodna sciezka A/B.")
  if(!all(c("cds",if(extract_lnc) "lnc" else "nc") %in% df$class)) stop("Brak wymaganej klasy. Wybierz recznie inna sciezke lub dane.")
  set.seed(seed)
  df$ID <- ave(df$class,df$class,FUN=function(x) paste(prefix,x,sprintf("%07d",seq_along(x)),sep="_"))
  df$stage <- ifelse(df$class=="ex","inference",NA_character_)
  df$tt_set <- ifelse(df$class=="ex","predict",NA_character_)
  df$keep <- TRUE; df$role <- ifelse(df$class=="ex","unknown",NA_character_)
  df$drop_reason <- NA_character_
  split_exact <- function(n,labels,prop) {
    if(n==0L) return(character())
    n1 <- round(n*prop); sample(c(rep(labels[1],n1),rep(labels[2],n-n1)))
  }
  for(cl in setdiff(unique(df$class),"ex")) {
    i <- which(df$class==cl); df$stage[i] <- split_exact(length(i),c("cp","ml"),cp_stage_prop)
    for(st in c("cp","ml")) {
      j <- i[df$stage[i]==st]; p <- if(st=="cp") cp_tt_prop else ml_tt_prop
      df$tt_set[j] <- split_exact(length(j),c("train","test"),p)
    }
  }
  df$cohort <- ifelse(df$class=="ex","ex_inference",paste(df$class,df$stage,df$tt_set,sep="_"))
  df$role[df$class=="cds"] <- "neg"
  if(extract_lnc) {
    df$role[df$stage=="cp" & df$class %in% c("lnc","nc_other")] <- "pos"
    df$role[df$stage=="ml" & df$class=="lnc"] <- "pos"
    i <- which(df$stage=="ml" & df$class=="nc_other")
    if(ml_neg_mode=="cds_only") {df$keep[i]<-FALSE; df$drop_reason[i]<-"ml_nc_other_excluded"} else df$role[i]<-"neg"
  } else df$role[df$class=="nc"] <- "pos"
  for(st in c("cp","ml")) {
    limit <- if(st=="cp") max_n_cp else max_n_ml; p <- if(st=="cp") cp_tt_prop else ml_tt_prop
    for(tt in c("train","test")) {
      lim <- floor(limit*(if(tt=="train") p else (1-p)))
      pos <- which(df$stage==st & df$tt_set==tt & df$role=="pos" & df$keep)
      neg <- which(df$stage==st & df$tt_set==tt & df$role=="neg" & df$keep)
      target <- if(tt=="train" && balance_train) min(length(pos),length(neg),lim) else min(max(length(pos),length(neg)),lim)
      for(i in list(pos,neg)) if(length(i)>target) {
        drop <- i[sample.int(length(i),length(i)-target)]
        df$keep[drop] <- FALSE; df$drop_reason[drop] <- "cap_or_balance"
      }
    }
  }
  for(st in c("cp","ml")) for(tt in c("train","test")) {
    roles <- df$role[df$stage==st & df$tt_set==tt & df$keep]
    if(!all(c("pos","neg") %in% roles)) stop("Pusta klasa w ",st,"_",tt,". Zmien proporcje/limity/dane.")
  }
  all_seqs <- DNAStringSet(df$seq_char); names(all_seqs) <- df$ID; df$seq_char <- NULL
  list(seqs=all_seqs,legend=df)
}

# 14: Dawna funkcja 4. which() zapobiega zwracaniu NA z wykluczonych rol.
extract_subset <- function(wynik,stage_name,tt_name=NULL,role_name=NULL,only_kept=TRUE) {
  df <- wynik$legend; mask <- df$stage==stage_name
  if(!is.null(tt_name)) mask <- mask & df$tt_set==tt_name
  if(!is.null(role_name)) mask <- mask & df$role==role_name
  if(only_kept) mask <- mask & df$keep
  wynik$seqs[df$ID[which(mask)]]
}

# 15: Wszystkie 10 par 5 zbiorow, ID oraz IDENTYCZNE SEKWENCJE.
validate_subsets <- function(sets) {
  pairs <- combn(names(sets),2,simplify=FALSE)
  do.call(rbind,lapply(pairs,function(p) data.frame(set1=p[1],set2=p[2],
    shared_ids=length(intersect(names(sets[[p[1]]]),names(sets[[p[2]]]))),
    shared_sequences=length(intersect(as.character(sets[[p[1]]]),as.character(sets[[p[2]]]))))))
}

# 16: Eksport TSV i FASTA. Puste FASTA zapisujemy jako pusty plik gzip.
write_tsv <- function(x,file) write.table(x,file,sep="\t",row.names=FALSE,quote=TRUE,na="NA",fileEncoding="UTF-8")
write_fasta <- function(x,file) {
  if(length(x)) writeXStringSet(x,file,compress=TRUE) else {con<-gzfile(file,"wb");close(con)}
}

# V3: wspolny pomiar czasu etapow; argument expr jest obliczany tylko raz.
timed_step <- function(label, expr) {
  if (progress_enabled) { cat("\n",label," ...\n",sep=""); flush.console() }
  t <- proc.time()[["elapsed"]]; result <- force(expr)
  if (progress_enabled) { cat(sprintf("  OK | %.2f s\n",proc.time()[["elapsed"]]-t)); flush.console() }
  result
}

# V3: porcje ograniczone liczba rekordow i suma nt. Jedna bardzo dluga
# sekwencja moze przekroczyc limit nt, ale nigdy nie tniemy jej na fragmenty.
# Nie zmienia to biologii ORF, PAS ani k-merow na granicach porcji.
sequence_chunks <- function(s, max_n=3000L, max_nt=8000000) {
  w <- as.double(width(s)); n <- length(w)
  if(!n) return(list())
  # Sumy prefiksowe; przeszukujemy co najwyzej max_n kolejnych rekordow.
  # Petla biegnie po porcjach, nie po transkryptach. Limit nt jest scisly,
  # z wyjatkiem pojedynczego transkryptu dluzszego niz max_nt.
  cs <- c(0,cumsum(w)); i <- 1L; k <- 0L; ans <- list()
  while (i<=n) {
    j <- min(n,i+max_n-1L)
    candidates <- i:j
    take <- max(1L,findInterval(cs[i]+max_nt,cs[candidates+1L]))
    k <- k+1L; ans[[k]] <- i:(i+take-1L); i <- i+take
  }
  ans
}

# 00A. Male testy regresyjne (opcjonalne; nie wymagaja plikow referencyjnych).
if(run_self_tests) {
  stopifnot(identical(extract_transcript_ids(c("ENSMUST0000001.2 cdna","MSTRG.123.1 x"),"ensembl"),c("ENSMUST0000001","MSTRG.123.1")))
  toy <- DNAStringSet(c("ATGAAATAA","ATGAAATAGATGCCCCCTTGA","AAAAAA","ATG"))
  o <- calc_longest_orf_bp(toy)
  stopifnot(o$max_orf_len[1]==9,o$max_orf_len[3]==0,o$max_orf_len[4]==0)
  stopifnot(o$orf_end[1]==9)
  # Przykladowe eksony minus: pierwszy w kierunku transkrypcji jest na prawo.
  ex0 <- GRangesList(t=GRanges("1",IRanges::IRanges(c(10,30),c(14,39)),strand="-"))
  z <- calc_exon_features(ex0)
  stopifnot(z$tx_length==15,z$first_exon_len==10,z$last_exon_len==5,z$total_intron_len==15)
  # Genomic context: same-strand exon overlap, intronic, antisense, adjacent.
  qex <- GRangesList(a=GRanges("1",IRanges::IRanges(10,12),strand="+"),
                     b=GRanges("1",IRanges::IRanges(21,25),strand="+"),
                     c=GRanges("1",IRanges::IRanges(11,13),strand="-"),
                     d=GRanges("1",IRanges::IRanges(40,42),strand="+"))
  pce <- GRanges("1",IRanges::IRanges(c(10,30),c(15,39)),strand="+")
  pcr <- GRanges("1",IRanges::IRanges(10,39),strand="+"); names(pcr)<-"pc_gene"
  pc0 <- list(genes=pcr,tx=pcr,exons=pce,introns=GRanges("1",IRanges::IRanges(16,29),strand="+"),contigs="1")
  g0 <- calc_genomic_features(get_tx_ranges(qex),qex,data.frame(transcript_id=names(qex),gene_id=NA_character_),pc0,"experimental")
  stopifnot(identical(g0$genomic_class,c("sense_overlapping","sense_intronic","antisense","intergenic")),
            g0$dist_to_pc[4]==0,!g0$overlap_pc[4],g0$exonic_overlap_bp[1]==3)
  # Rozne sekwencje bez konfliktu, 40 PC + 40 lnc + 40 nc, i EX identyczny z PC.
  toyseq <- vapply(0:120,function(i) paste0(c("A","C","G","T")[(i%/%4^(0:4))%%4+1],collapse=""),character(1))
  toy_df <- data.frame(record_id=paste0("r",1:121),old_name=paste0("r",1:121),
    raw_group=c(rep("cds",40),rep("lnc_ext",40),rep("nc",40),"ex"),seq_char=c(toyseq[1:120],toyseq[1]))
  for(path_a in c(TRUE,FALSE)) {
    d <- load_and_dedup(extract_lnc=path_a,df=toy_df)
    w <- split_datasets(d,extract_lnc=path_a)
    ss <- list(cp_train=extract_subset(w,"cp","train"),cp_test=extract_subset(w,"cp","test"),
               ml_train=extract_subset(w,"ml","train"),ml_test=extract_subset(w,"ml","test"),ex_pred=extract_subset(w,"inference"))
    v <- validate_subsets(ss); stopifnot(all(v$shared_ids==0),all(v$shared_sequences==0),length(ss$ex_pred)==0)
  }
  # V3: test regresji ORF przeciw niezaleznemu algorytmowi z wersji 2.
  # Wolniejszy algorytm istnieje TYLKO w bloku testowym.
  orf_reference <- function(dna_set, min_nt = 6L) {
    starts <- start(vmatchPattern("ATG", dna_set))
    stops1 <- start(vmatchPattern("TAA", dna_set))
    stops2 <- start(vmatchPattern("TAG", dna_set))
    stops3 <- start(vmatchPattern("TGA", dna_set))
    rows <- lapply(seq_along(dna_set), function(i) {
      st <- starts[[i]]; sp <- sort(unique(c(stops1[[i]], stops2[[i]], stops3[[i]])))
      len <- integer(); beg <- integer()
      for (frame in 0:2) {
        a <- st[st %% 3L == frame]; b <- sp[sp %% 3L == frame]
        if (!length(a) || !length(b)) next
        j <- findInterval(a, b) + 1L
        good <- j <= length(b)
        l <- b[j[good]] - a[good] + 3L
        keep <- l >= min_nt
        len <- c(len, l[keep]); beg <- c(beg, a[good][keep])
      }
      if (!length(len)) return(c(max_orf_len=0L, orf_count=0L, orf_start=NA_integer_, orf_end=NA_integer_))
      best <- order(-len, beg)[1L]
      c(max_orf_len=len[best], orf_count=length(len), orf_start=beg[best], orf_end=beg[best]+len[best]-1L)
    })
    as.data.frame(do.call(rbind, rows))
  }
  set.seed(42)
  random_dna <- DNAStringSet(vapply(seq_len(200),function(i) paste(sample(c("A","C","G","T","N"),sample(1:500,1),replace=TRUE),collapse=""),character(1)))
  test_dna <- c(toy,DNAStringSet(c("ATGATGTAA","ATGTAATAG","CATGAAATAA","ATGAAATGAATGAAATAA","AATAAAAATTAAA","AAAA","N")),random_dna)
  for(m in c(6L,9L,30L)) stopifnot(isTRUE(all.equal(calc_longest_orf_bp(test_dna,m),orf_reference(test_dna,m),check.attributes=FALSE)))
  test_ids <- paste0("test",seq_along(test_dna))
  full <- calc_sequence_features(test_dna,test_ids,kmers=TRUE)
  chunks <- sequence_chunks(test_dna,17L,1000)
  part <- do.call(rbind,lapply(chunks,function(i) calc_sequence_features(test_dna[i],test_ids[i],kmers=TRUE)))
  rownames(part) <- NULL
  stopifnot(isTRUE(all.equal(full,part,check.attributes=FALSE)),identical(unlist(chunks,use.names=FALSE),seq_along(test_dna)),
            all(lengths(chunks)<=17L),all(vapply(chunks,function(i) sum(width(test_dna[i]))<=1000,logical(1))))
  pas_test <- calc_sequence_features(DNAStringSet(c("AATAAACC","ATTAAA","CCCC","NNNN","AAA")),paste0("p",1:5))
  stopifnot(identical(pas_test$pas_present,c(TRUE,TRUE,FALSE,FALSE,FALSE)),
    isTRUE(all.equal(pas_test$pas_distance_from_end,c(2,0,NA,NA,NA))),
    identical(as.integer(pas_test$terminal_A_run),c(0L,3L,0L,0L,3L)))
  a <- get_introns_by_tx(ex0); stopifnot(lengths(a)==1L,start(unlist(a))==15,end(unlist(a))==29)
  b <- get_introns_by_tx(qex); stopifnot(identical(names(b),names(qex)),all(lengths(b)==0L))
  # Sprawdzenie statystyk liczbowych, w tym pustych grup i median parzystych.
  v <- grouped_width_stats(c(2,4,9,5),c(1L,1L,1L,3L),4L)
  stopifnot(v$median[1]==4,v$total[2]==0,is.na(v$mean[2]),isTRUE(all.equal(v$sd[1],sd(c(2,4,9)))))
  stopifnot(isTRUE(all.equal(calc_exon_features(ex0),calc_exon_features(ex0,a))))
  # V4: syntetyczne introny. Rozne kontigi, plus/minus, wszystkie 3 motywy
  # canonical, noncanonical, N, za krotki intron, brak kontigu, poza zakresem,
  # nieznana nic, pusty transkrypt. Celowo mala porcja sprawdza jej granice.
  # Chr1: GTAAAG | CTAAAC | GCAAAG | ATAAAC | CCAATT | NNAAAG.
  # Na minusie CT...AC odpowiada GT...AG; identyczny zakres na plus nie.
  toy_genome <- DNAStringSet(c(chr1="GTAAAGCTAAACGCAAAGATAAACCCAATTNNAAAG",chr2="GTCCAG"))
  toy_intr <- GRangesList(
    plus=GRanges("chr1",IRanges::IRanges(c(1L,13L,19L),width=6L),strand="+"),
    minus=GRanges("chr1",IRanges::IRanges(7L,width=6L),strand="-"),
    mixed=GRanges(c("chr1","chr2","chr1"),IRanges::IRanges(c(25L,1L,31L),width=6L),strand="+"),
    short=GRanges("chr1",IRanges::IRanges(1L,width=3L),strand="+"),
    absent=GRanges("absent",IRanges::IRanges(1L,width=6L),strand="+"),
    outside=GRanges("chr1",IRanges::IRanges(35L,width=6L),strand="+"),
    unstranded=GRanges("chr1",IRanges::IRanges(1L,width=6L),strand="*"),
    mono=GRanges(),
    plus_noncanonical=GRanges("chr1",IRanges::IRanges(7L,width=6L),strand="+"))
  sp1 <- suppressWarnings(calc_splice_features(toy_intr,engine="biostrings",genome=toy_genome,chunk_size=2L))
  sp2 <- suppressWarnings(calc_splice_features(toy_intr,engine="biostrings",genome=toy_genome,chunk_size=100L))
  stopifnot(identical(sp1,sp2),identical(sp1$transcript_id,names(toy_intr)),
    identical(sp1$gt_ag,c(1L,1L,1L,0L,0L,0L,0L,0L,0L)),
    identical(sp1$gc_ag,c(1L,0L,0L,0L,0L,0L,0L,0L,0L)),
    identical(sp1$at_ac,c(1L,0L,0L,0L,0L,0L,0L,0L,0L)),
    identical(sp1$noncanonical_junctions,c(0L,0L,1L,0L,0L,0L,0L,0L,1L)),
    identical(sp1$unknown_junctions,c(0L,0L,1L,1L,1L,1L,1L,0L,0L)),
    isTRUE(all.equal(sp1$noncanonical_fraction,c(0,0,.5,NA,NA,NA,NA,NA,1))))
  sp_none <- calc_splice_features(toy_intr,engine="none",genome_file="NIE_ISTNIEJE.fa")
  stopifnot(identical(names(sp1),names(sp_none)),all(is.na(sp_none[,-1L])))
  sp_both_none <- run_splice_analysis(list(reference=toy_intr,experimental=toy_intr),
                                      "NIE_ISTNIEJE.fa","none",2L)
  stopifnot(identical(sp_both_none$reference,sp_none),identical(sp_both_none$experimental,sp_none))
  # Import FASTA.gz, normalizacja naglowkow i wspolny przebieg REF/EX.
  # Plik tymczasowy jest usuwany nawet po nieudanym tescie R (nie po awarii R).
  (function() {
    toy_file <- tempfile(fileext=".fa.gz")
    on.exit(unlink(toy_file),add=TRUE)
    g <- toy_genome; names(g) <- paste(names(g),"test description")
    writeXStringSet(g,toy_file,compress=TRUE)
    both <- suppressWarnings(run_splice_analysis(list(reference=toy_intr,experimental=toy_intr),toy_file,"biostrings",2L))
    stopifnot(identical(both$reference,sp1),identical(both$experimental,sp1))
    sm <- summarize_splice_features(sp1,toy_intr,"reference","biostrings")
    sn <- summarize_splice_features(sp_none,toy_intr,"reference","none")
    stopifnot(sm$resolved+sm$unknown==sm$introns,is.na(sn$canonical),is.na(sn$noncanonical_percent))
  })()
  if(self_test_rsamtools) (function() {
    toy_file <- tempfile(fileext=".fa")
    on.exit(unlink(c(toy_file,paste0(toy_file,".fai"),paste0(toy_file,".gzi"))),add=TRUE)
    writeXStringSet(toy_genome,toy_file)
    Rsamtools::indexFa(toy_file)
    sp_r <- suppressWarnings(calc_splice_features(toy_intr,toy_file,"rsamtools",chunk_size=2L))
    stopifnot(identical(sp1,sp_r))
  })()
  # V4: REF-only oraz REF+EX z rzeczywistym niepustym zbiorem EX.
  # Identyczny REF ma takie same ID i przydzialy CP/ML niezaleznie od dolaczenia
  # EX, przy tych samych ustawieniach i seed. EX identyczny z REF odpada.
  toy_ref <- toy_df[1:120,]
  toy_ex <- data.frame(record_id="r122",old_name="r122",raw_group="ex",seq_char=toyseq[121])
  for(path_a in c(TRUE,FALSE)) {
    only <- split_datasets(load_and_dedup(df=toy_ref,extract_lnc=path_a),extract_lnc=path_a)
    both <- split_datasets(load_and_dedup(df=rbind(toy_df,toy_ex),extract_lnc=path_a),extract_lnc=path_a)
    stopifnot(nrow(only$legend)==120L,length(extract_subset(only,"inference"))==0L,
      length(extract_subset(both,"inference"))==1L,
      identical(unname(as.character(extract_subset(both,"inference"))),toyseq[121]))
    for(stage in c("cp","ml")) for(tt in c("train","test")) {
      aa <- extract_subset(only,stage,tt); bb <- extract_subset(both,stage,tt)
      stopifnot(identical(names(aa),names(bb)),identical(as.character(aa),as.character(bb)))
    }
  }
  # V4.1: syntetyczne FASTA na dysku; test importu jednego i wielu plikow,
  # scalenia zgodnych kopii, odrzucenia konfliktu wersji/sekwencji i kontroli klas.
  local({
    a <- tempfile(fileext=".fa"); b <- tempfile(fileext=".fa")
    writeLines(c(">ENSMUST0000001.1 cdna","ATGAAATAA"),a)
    writeLines(c(">ENSMUST0000002.1 ncrna","ACCCACCC",">ENSMUST0000001.1 copy","ATGAAATAA"),b)
    one <- read_reference_fastas(a); both <- read_reference_fastas(c(a,b))
    stopifnot(length(one$seqs)==1L,length(both$seqs)==2L,
      identical(unname(as.character(both$seqs)),c("ATGAAATAA","ACCCACCC")),
      both$fasta_files[1L]==paste(a,b,sep=" | "),both$fasta_files[2L]==b)
    m <- data.frame(transcript_id=paste0("ENSMUST000000",1:3),
      transcript_biotype=c("protein_coding","lncRNA","lncRNA"))
    f <- data.frame(transcript_id=extract_transcript_ids(names(both$seqs),"ensembl"))
    cov <- calc_fasta_coverage(f,m,"reference")
    stopifnot(sum(cov$matched_n)==2L,sum(cov$gtf_only_n)==1L,sum(cov$fasta_only_n)==0L,
      cov$coverage_percent[cov$transcript_biotype=="lncRNA"]==50,
      check_reference_classes(cov,TRUE)[["lnc"]]==1L)
    bad <- calc_fasta_coverage(f[1L,,drop=FALSE],m,"reference")
    stopifnot(inherits(try(check_reference_classes(bad,TRUE),silent=TRUE),"try-error"))
    writeLines(c(">ENSMUST0000001.2","ATGAAATAA"),b)
    stopifnot(inherits(try(read_reference_fastas(c(a,b)),silent=TRUE),"try-error"))
    writeLines(c(">ENSMUST0000001.1","ATGCCCTAA"),b)
    stopifnot(inherits(try(read_reference_fastas(c(a,b)),silent=TRUE),"try-error"))
    writeLines(c(">ENSMUST0000001.1","ATGAAATAA",">ENSMUST0000001.1","ATGAAATAA"),b)
    stopifnot(inherits(try(read_reference_fastas(b),silent=TRUE),"try-error"))
    unlink(c(a,b))
  })
  message("Self-tests OK (ID, ORF, porcje, struktura, A/B, REF/EX, splicing, multi-FASTA, pokrycie GTF).")
  if (self_tests_only) stop("Testy zakonczone pomyslnie; self_tests_only=TRUE, bez importu danych.",call.=FALSE)
}

# ///////////////////////////////////////////////////////////////////////////////
# KOD GLOWNY — wykonuj kolejne sekcje w RStudio
# ///////////////////////////////////////////////////////////////////////////////

# 01. DANE WEJSCIOWE -----------------------------------------------------------
# Nowy Ensembl: przypiety przyklad z plikow uzytkownika. Pobieranie opcjonalne:
# base_url <- "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/635/9/ensembl/2026_04/geneset/"
# options(timeout=1200)
# download.file(paste0(base_url,"cdna.fa.bgz"), "cdna.fa.bgz", mode="wb")
# download.file(paste0(base_url,"genes.gtf.bgz"), "genes.gtf.bgz", mode="wb")
# Nie uzywaj niepelnego adresu https://ensembl.org jako adresu FASTA/GTF.

seqs <- list(); group_hint <- list(); anno <- list(); id_mode <- list(); fasta_files <- list()
if(input_mode=="ensembl") {
  ref_input <- read_reference_fastas(ref_fasta,ref_id_mode)
  seqs$reference <- ref_input$seqs; fasta_files$reference <- ref_input$fasta_files
  group_hint$reference <- rep(NA_character_,length(seqs$reference)); rm(ref_input)
  anno$reference <- rtracklayer::import(ref_gtf,format="gtf"); id_mode$reference <- ref_id_mode
} else if(input_mode=="separate") {
  # LncRNA moze byc celowo obecne w obu plikach NC i lnc. Po identyfikatorach
  # scalimy te rekordy w 02, a nie potraktujemy jako dwa rozne transkrypty.
  files <- list(protein_coding=pc_file,non_coding=nc_file,lncRNA=lnc_file,ambiguous=ambiguous_file)
  ss <- lapply(files,read_fasta)
  seqs$reference <- do.call(c,unname(ss)); group_hint$reference <- rep(names(ss),lengths(ss))
  fasta_files$reference <- unlist(Map(function(file,s) rep(file,length(s)),files,ss),use.names=FALSE)
  if(!is.null(separate_gtf)) anno$reference <- rtracklayer::import(separate_gtf,format="gtf")
  id_mode$reference <- ref_id_mode
} else if(!is.null(ref_gtf)) {
  # Tylko EX: referencyjny GTF sluzy jako tlo genomowe, nie do nadawania klas EX.
  anno$reference <- rtracklayer::import(ref_gtf,format="gtf"); id_mode$reference <- ref_id_mode
}
if(include_ex || input_mode=="experimental_only") {
  seqs$experimental <- read_fasta(ex_file); group_hint$experimental <- rep(NA_character_,length(seqs$experimental))
  fasta_files$experimental <- rep(ex_file,length(seqs$experimental))
  anno$experimental <- rtracklayer::import(ex_gtf,format="gtf"); id_mode$experimental <- ex_id_mode
}
print(vapply(seqs,length,integer(1)))


# Inspekcja dla wszystkich AKTYWNYCH zrodel (rowniez przy REF-only).
# sum(as.double(...)) unika przepelnienia sumy szerokosci >2^31-1 nt.
print(vapply(anno,length,integer(1)))
print(vapply(seqs,function(s) sum(as.double(width(s))),numeric(1)))
print(lapply(seqs,function(s) summary(width(s))))
cat("FASTA: ",format(object.size(seqs),units="auto"),
    " | GTF: ",format(object.size(anno),units="auto"),"\n",sep="")
print(gc())


# 02. IDENTYFIKATORY ----------------------------------------------------------
# Nie zmieniamy oryginalnych naglowkow FASTA. record_id jest kluczem globalnym
# (source + transcript_id), a transcript_id jest kluczem wewnatrz danego zrodla.
fasta_meta <- list()
for(src in names(seqs)) {
  s <- seqs[[src]]; h <- names(s); ids <- extract_transcript_ids(h,id_mode[[src]])
  if(anyNA(ids)||any(!nzchar(ids))) stop("Brak ID FASTA: ",src)
  hints <- group_hint[[src]]; origins <- fasta_files[[src]]
  stopifnot(length(origins)==length(s))
  # Zachowanie wersji/oryginalnego ID. Nie gubimy naglowka w eksportach biologicznych.
  if(anyDuplicated(ids)) {
    if(input_mode!="separate" || src!="reference") stop("Nieunikalne ID FASTA po normalizacji: ",src)
    groups <- split(which(duplicated(ids)|duplicated(ids,fromLast=TRUE)), ids[duplicated(ids)|duplicated(ids,fromLast=TRUE)])
    for(i in groups) {
      if(length(unique(as.character(s[i])))!=1L) stop("Ten sam ID, rozne sekwencje: ",ids[i[1]])
      # Tolerujemy tylko identyczne NC/lnc i duplikat tej samej klasy.
      if(length(unique(hints[i]))>1L && !all(hints[i] %in% c("non_coding","lncRNA"))) stop("Konflikt plikow dla ID: ",ids[i[1]])
      if(any(hints[i]=="lncRNA")) hints[i] <- "lncRNA"
      origins[i[1L]] <- paste(unique(origins[i]),collapse=" | ")
    }
    keep <- !duplicated(ids); s <- s[keep]; h <- h[keep]; ids <- ids[keep]; hints <- hints[keep]; origins <- origins[keep]
  }
  seqs[[src]] <- s; fasta_files[[src]] <- origins
  fasta_meta[[src]] <- data.frame(transcript_id=ids,fasta_id=extract_transcript_ids(h,"exact"),original_header=h,fasta_files=origins,
                                  group_hint=hints,has_fasta=TRUE,stringsAsFactors=FALSE)
  print(head(fasta_meta[[src]],2))
}

# 03. ADNOTACJA I METADANE -----------------------------------------------------
# Ten sam kod dla Ensembl i StringTie. GFF3 ma ID/Parent; najpierw przekonwertuj
# na GTF (np. gffread -T -o genes.gtf genes.gff3) i sprawdz zachowanie biotypow.
tx_meta <- list(); exons_by_tx <- list(); tx_ranges <- list(); introns_by_tx <- list()
for(src in names(anno)) {
  tx_meta[[src]] <- timed_step(paste(src,"03: metadane"),get_tx_meta(anno[[src]],id_mode[[src]]))
  exons_by_tx[[src]] <- timed_step(paste(src,"03: eksony"),get_exons_by_tx(anno[[src]],id_mode[[src]],structure_engine))
  tx_ranges[[src]] <- timed_step(paste(src,"03: zakresy"),get_tx_ranges(exons_by_tx[[src]]))
  introns_by_tx[[src]] <- timed_step(paste(src,"03: introny"),get_introns_by_tx(exons_by_tx[[src]]))
  # Kontrola TxDb / GTF: nie pozwalamy na ciche zgubienie eksonow przez adapter.
  if(!all(names(exons_by_tx[[src]]) %in% tx_meta[[src]]$transcript_id)) stop("Eksony bez metadanych: ",src)
  print(head(tx_meta[[src]]))
}

# 03A. KONTROLA KOMPLETNOSCI FASTA/GTF (V4.1) ----------------------------------
# Raport powstaje PRZED kosztownymi cechami 04-07, takze przy run_split=FALSE.
# Zapisujemy go od razu, nawet jezeli nastepna kontrola zatrzyma analize.
# Przy separate bez GTF raport opisuje FASTA-only; w experimental_only nie
# dolaczamy referencji bez FASTA (sluzy wylacznie jako tlo genomowe).
fasta_coverage <- do.call(rbind,lapply(names(fasta_meta),function(src)
  calc_fasta_coverage(fasta_meta[[src]],tx_meta[[src]],src)))
rownames(fasta_coverage) <- NULL
print(fasta_coverage,row.names=FALSE)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if(write_outputs) write_tsv(fasta_coverage,file.path(out_dir,"fasta_gtf_coverage.tsv"))
# Brak klasy jest bledem TYLKO przy przygotowaniu kohort i w trybie ensembl.
# Tryb separate zachowuje wlasna logike etykiet (rowniez bez GTF).
if(run_split && input_mode=="ensembl")
  check_reference_classes(fasta_coverage[fasta_coverage$source=="reference",],extract_lnc)
# Przy bledzie zmien pliki, wybierz nowy/pusty out_dir i powtorz od 00.
# Nie wczytuj checkpointow obliczonych dla innego zbioru FASTA.

# 04. CECHY SEKWENCYJNE -------------------------------------------------------
# Szybki benchmark przed cala analiza (opcjonalnie, nie wykonuje sie sam):
# i <- seq_len(min(1000L,length(seqs$reference)))
# system.time(test <- calc_sequence_features(seqs$reference[i],fasta_meta$reference$transcript_id[i],pas_window,min_orf_nt,calc_kmers))
seq_features <- list()
for(src in names(seqs)) {
  s <- seqs[[src]]; ids <- fasta_meta[[src]]$transcript_id
  chunks <- sequence_chunks(s,sequence_chunk_size,sequence_chunk_nt)
  results <- vector("list",length(chunks)); t0 <- proc.time()[["elapsed"]]
  cat("\n04: ",src," | ",length(s)," transkryptow | ",sum(as.double(width(s)))," nt | ",length(chunks)," porcji\n",sep="")
  for(i in seq_along(chunks)) {
    idx <- chunks[[i]]; t1 <- proc.time()[["elapsed"]]
    results[[i]] <- calc_sequence_features(s[idx],ids[idx],pas_window,min_orf_nt,calc_kmers)
    if(progress_enabled) {
      cat(sprintf("  %d/%d | %d/%d (%.1f%%) | blok %.2f s | razem %.1f s\n",
        i,length(chunks),max(idx),length(s),100*max(idx)/length(s),proc.time()[["elapsed"]]-t1,proc.time()[["elapsed"]]-t0))
      flush.console()
    }
  }
  seq_features[[src]] <- do.call(rbind,results); rownames(seq_features[[src]]) <- NULL
  stopifnot(identical(seq_features[[src]]$transcript_id,ids))
  print(head(seq_features[[src]][,1:8]))
  # Checkpoint po calym zrodle. W razie przerwania w tej samej sesji lista
  # results zachowuje ukonczone porcje; nie kasuj jej przed ewentualnym wznowieniem.
  # Reczny odczyt RDS tylko po sprawdzeniu zgodnosci FASTA, ID i parametrow!
  if(write_outputs && save_sequence_checkpoint)
    saveRDS(list(features=seq_features[[src]],ids=ids,parameters=list(pas_window=pas_window,min_orf_nt=min_orf_nt,calc_kmers=calc_kmers)),
      file.path(out_dir,paste0("checkpoint_sequence_",src,".rds")),compress=FALSE)
  rm(results,chunks,s,ids)
}

# 05. ARCHITEKTURA EKSON-INTRON ------------------------------------------------
struct_features <- list()
for(src in names(exons_by_tx)) {
  e <- exons_by_tx[[src]]; tx <- tx_ranges[[src]]
  x <- timed_step(paste(src,"05: statystyki eksonow/intronow"),calc_exon_features(e,introns_by_tx[[src]]))
  loc <- data.frame(transcript_id=names(tx),seqname=as.character(seqnames(tx)),strand=as.character(strand(tx)),
                    tx_start=start(tx),tx_end=end(tx))
  x <- join_by_tx(x,loc)
  m <- tx_meta[[src]]; g <- m$gene_id
  tab <- table(g[!is.na(g)]); m$isoforms_per_gene <- as.integer(tab[match(g,names(tab))])
  struct_features[[src]] <- join_by_tx(x,m[,c("transcript_id","isoforms_per_gene")])
  print(head(struct_features[[src]],2))
}

# 06. KONTEKST GENOMOWY -------------------------------------------------------
pc_ref <- NULL
if("reference" %in% names(tx_meta)) pc_ref <- timed_step("06: referencja PC",build_pc_reference(tx_meta$reference,exons_by_tx$reference,tx_ranges$reference,introns_by_tx$reference))
if(is.null(pc_ref)) message("Brak wiarygodnej referencji PC: kontekst genomowy = NA.")
genomic_features <- list()
for(src in names(tx_ranges)) {
  genomic_features[[src]] <- timed_step(paste(src,"06: kontekst genomowy"),calc_genomic_features(tx_ranges[[src]],exons_by_tx[[src]],tx_meta[[src]],pc_ref,src))
  print(table(genomic_features[[src]]$genomic_class,useNA="ifany"))
}

# 07. OPCJONALNE CECHY Z GENOMU DNA -------------------------------------------
# Do motywow intronowych potrzebujemy genomu, NIE FASTA cDNA. W V4
# biostrings czyta oryginalny plik .fa.gz bez przygotowania indeksow.
# Rsamtools pozostaje opcja: raz przygotowany FAI (i GZI dla BGZF).
# Nie uruchamiamy tutaj indexFa(), gunzip() ani bgzip(); nie modyfikujemy wejsc.
# W trybie none nie czytamy genomu i nie wywolujemy scanFa(). Nadal zostaja
# wszystkie kolumny splicingu (NA). GTF i kontekst PC nie sa wylaczane.
# Jezeli wystapil problem z pamiecia, ustaw genome_engine="none" i wykonaj
# ponownie 07-16 w tej samej sesji TYLKO gdy obiekty 01-06 sa nadal dostepne.
# Jezeli zmieniasz out_dir w trakcie sesji, utworz nowy katalog przed eksportem:
# dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
# W takim przypadku checkpointy sekwencji z 04 pozostaja w starym katalogu.
# Dla nowego przebiegu wybierz nowy/pusty out_dir przed zapisem.
# Jezeli sesja R ulegla awarii, potrzebny jest restart i odtworzenie 01-06.
# Dalsze moduly (repeats, BAM, MFE, unikalne junctions) pozostaja poza V4 robocza.
if(is.null(genome_fasta)) genome_engine <- "none"
cat("07: silnik genomowy = ",genome_engine,"\n",sep="")
splice_features <- run_splice_analysis(introns_by_tx,genome_fasta,genome_engine,splice_chunk_size)
# Po wyjsciu z funkcji genom Biostrings nie jest potrzebny w dalszych etapach.
# gc() odzyskuje nieosiagalne obiekty; system nie musi oddac calej pamieci od razu.
print(gc())
splice_summary <- NULL
if(length(splice_features)) {
  splice_summary <- do.call(rbind,lapply(names(splice_features),function(src)
    summarize_splice_features(splice_features[[src]],introns_by_tx[[src]],src,genome_engine)))
  rownames(splice_summary) <- NULL
  print(splice_summary,row.names=FALSE)
}

# 08. INTEGRACJA CECH ---------------------------------------------------------
# Pelna unia ID z FASTA i GTF: GTF-only pozostaje opisowo, ale bez sekwencji
# nie moze wejsc do zbioru treningowego. FASTA-only nie znika z tabeli.
# Szablony kolumn dla braku GTF (NA, nigdy nie uznajemy braku za zero).
meta_fields <- c("gtf_transcript_id","gene_id","gene_name","gene_biotype","transcript_biotype","transcript_version")
struct_fields <- c("exon_count","tx_length","mean_exon_len","median_exon_len","min_exon_len","max_exon_len","sd_exon_len",
  "first_exon_len","last_exon_len","genomic_span","exon_density","first_exon_fraction","last_exon_fraction","intron_count",
  "total_intron_len","mean_intron_len","median_intron_len","min_intron_len","max_intron_len","sd_intron_len","exon_intron_ratio",
  "tx_start","tx_end","isoforms_per_gene")
genome_fields <- c("dist_to_pc","dist_to_nearest_tss","exonic_overlap_bp","exonic_overlap_fraction","overlapping_pc_genes")
splice_fields <- c("gt_ag","gc_ag","at_ac","canonical_junctions","noncanonical_junctions","unknown_junctions","noncanonical_fraction")
features <- list()
for(src in names(seqs)) {
  f <- fasta_meta[[src]]; m <- tx_meta[[src]]
  ids <- union(f$transcript_id,if(!is.null(m))m$transcript_id else character())
  x <- data.frame(transcript_id=ids,source=src,record_id=paste(src,ids,sep="::"))
  x <- join_by_tx(x,f); x$has_fasta <- x$has_fasta %in% TRUE
  x$has_gtf <- ids %in% (if(!is.null(m)) m$transcript_id else character())
  for(y in list(m,seq_features[[src]],struct_features[[src]],genomic_features[[src]],splice_features[[src]]))
    if(!is.null(y)) x <- join_by_tx(x,y)
  for(nm in setdiff(c(meta_fields,"seqname","strand","nearest_pc_gene_id","genomic_class","neighbor_orientation"),names(x))) x[[nm]] <- NA_character_
  for(nm in setdiff(c(struct_fields,genome_fields,splice_fields),names(x))) x[[nm]] <- NA_real_
  for(nm in setdiff(c("overlap_pc","overlap_same_strand","overlap_opposite_strand"),names(x))) x[[nm]] <- NA
  x$annotation_status <- ifelse(!x$has_gtf,"annotation_missing",ifelse(!x$has_fasta,"fasta_missing","matched"))
  x$length_difference <- x$seq_length-x$tx_length
  x$length_match <- ifelse(is.na(x$length_difference),NA,x$length_difference==0)
  # Version checking after Ensembl suffix normalization. A mismatch is not
  # resolved by stripping the version; it must be investigated.
  fv <- ifelse(grepl("^ENS[A-Z]*T[0-9]+\\.[0-9]+$",x$fasta_id),sub("^.*\\.","",x$fasta_id),NA_character_)
  gv <- x$transcript_version
  graw <- x$gtf_transcript_id
  suffix <- ifelse(grepl("^ENS[A-Z]*T[0-9]+\\.[0-9]+$",graw),sub("^.*\\.","",graw),NA_character_)
  gv[is.na(gv)] <- suffix[is.na(gv)]
  x$version_match <- ifelse(is.na(fv)|is.na(gv),NA,fv==gv)
  if(any(x$version_match %in% FALSE)) warning("Niezgodne wersje ID w ",src,"; rekordy nie trafia do treningu.")
  # Identyczne sekwencje moga miec rozne locus; zachowujemy kazdy rekord.
  a <- as.character(seqs[[src]])
  x$sequence_group <- match(a,a)[match(x$transcript_id,f$transcript_id)]
  features[[src]] <- x
  print(table(x$annotation_status,useNA="ifany")); print(table(x$length_match,useNA="ifany"))
}
# Rbind z uzupelnieniem tych samych kolumn w identycznej kolejnosci.
all_cols <- unique(unlist(lapply(features,names)))
features <- lapply(features,function(x){for(nm in setdiff(all_cols,names(x))) x[[nm]]<-NA; x[,all_cols,drop=FALSE]})
final_features <- do.call(rbind,unname(features)); rownames(final_features)<-NULL
stopifnot(!anyDuplicated(final_features$record_id))

# 09. KLASYFIKACJA BIOTYPOW ---------------------------------------------------
# Biotyp transkryptu, a NIE biotyp genu, decyduje o klasie referencyjnej.
# Brak adnotacji i niejasna adnotacja to rozne przypadki.
# W starym wariancie osobne FASTA moga byc jedynym zrodlem etykiet; jest to
# zalozenie uzytkownika, zapisane w label_source, nigdy nie wnioskujemy z ORF.
b <- final_features$transcript_biotype; ref <- final_features$source=="reference"
final_features$main_class <- "unknown"
final_features$is_lncrna <- NA
final_features$label_source <- NA_character_
annotated <- ref & final_features$has_gtf & !is.na(b) & nzchar(b)
final_features$main_class[annotated] <- "ambiguous"
final_features$main_class[annotated & b %in% pc_biotypes] <- "protein_coding"
final_features$main_class[annotated & b %in% all_ncrna_biotypes] <- "non_coding"
final_features$is_lncrna[annotated] <- b[annotated] %in% lncrna_biotypes
final_features$label_source[annotated] <- "transcript_biotype"
final_features$main_class[ref & !annotated] <- "annotation_missing"
hint <- final_features$group_hint
fallback <- ref & !annotated & !is.na(hint) & input_mode=="separate"
final_features$main_class[fallback] <- ifelse(hint[fallback]=="lncRNA","non_coding",hint[fallback])
# Gdy mamy NC bez GTF/lnc_file, dopuszczamy jawny, konfigurowalny regex naglowka.
legacy_lnc_pattern <- "lncRNA|lnc_rna|lincRNA|long_non_coding"
final_features$is_lncrna[fallback] <- hint[fallback]=="lncRNA" |
  (hint[fallback]=="non_coding" & grepl(legacy_lnc_pattern,final_features$original_header[fallback],ignore.case=TRUE))
final_features$label_source[fallback] <- "user_fasta_group_or_header"
# Rozbieznosci GTF vs plik wejsciowy sa raportowane. GTF ma pierwszenstwo.
final_features$label_conflict <- annotated & !is.na(hint) &
  (ifelse(hint=="lncRNA","non_coding",hint)!=final_features$main_class |
     (hint=="lncRNA" & !final_features$is_lncrna))
final_features$label_conflict[is.na(final_features$label_conflict)]<-FALSE
if(any(final_features$label_conflict)) warning("Niezgodne etykiety GTF / FASTA: wykluczone z treningu.")
final_classified <- final_features
final_features_extended <- final_features # alias zgodny z poprzednim skryptem
print(table(final_features$main_class,useNA="ifany")); print(table(b[ref & final_features$main_class=="ambiguous"],useNA="ifany"))

# 10. STATYSTYKI OPISOWE I WYKRESY ---------------------------------------------
# lncRNA jest podzbiorem non_coding. Nie dodajemy jej liczebnosci drugi raz.
rf <- final_features[ref,,drop=FALSE]
group_masks <- list(protein_coding=rf$main_class=="protein_coding",non_coding=rf$main_class=="non_coding",
                    lncRNA=rf$is_lncrna %in% TRUE,ambiguous=rf$main_class=="ambiguous")
summary_table <- do.call(rbind,lapply(names(group_masks),function(g)
  calc_group_stats(rf[group_masks[[g]],,drop=FALSE],g,sum(rf$has_fasta))))
summary_table$relation <- c("disjoint","contains lncRNA","subset of non_coding","disjoint")
print(summary_table,row.names=FALSE)
experimental_summary <- NULL
if (any(!ref)) {
  ef <- final_features[!ref,,drop=FALSE]
  experimental_summary <- calc_group_stats(ef,"experimental_unknown",sum(ef$has_fasta))
  print(experimental_summary,row.names=FALSE)
}
# Dodatkowa kontrola pokrycia: cztery grupy nie obejmuja annotation_missing.
print(table(rf$main_class,rf$has_fasta,useNA="ifany"))
if(make_plots && nrow(rf)) {
  pdat <- rf[rf$has_fasta & rf$main_class %in% c("protein_coding","non_coding","ambiguous"),]
  pdat$plot_class <- ifelse(pdat$is_lncrna %in% TRUE,"lncRNA",ifelse(pdat$main_class=="non_coding","nc_other",pdat$main_class))
  # Rozlaczne klasy na wykresach, w odroznieniu od hierarchicznych tabel.
  plot_vars <- c("seq_length","exon_count","gc_content","orf_coverage")
  for(v in plot_vars) {
    d <- pdat[is.finite(pdat[[v]]),,drop=FALSE]; if(!nrow(d)) next
    d$value <- if(v=="seq_length") log10(d[[v]]) else d[[v]]
    p <- ggplot2::ggplot(d,ggplot2::aes(x=plot_class,y=value,fill=plot_class))+
      ggplot2::geom_boxplot(outlier.alpha=.15)+ggplot2::theme_bw()+
      ggplot2::labs(x=NULL,y=if(v=="seq_length") "log10(seq_length)" else v)+ggplot2::theme(legend.position="none")
    print(p); if(write_outputs) ggplot2::ggsave(file.path(out_dir,paste0(v,".png")),p,width=8,height=5,dpi=160)
  }
  d <- pdat[!is.na(pdat$genomic_class),,drop=FALSE]
  if(nrow(d)) {
    p <- ggplot2::ggplot(d,ggplot2::aes(x=plot_class,fill=genomic_class))+ggplot2::geom_bar(position="fill")+
      ggplot2::theme_bw()+ggplot2::labs(x=NULL,y="Fraction",fill="Genomic class")
    print(p); if(write_outputs) ggplot2::ggsave(file.path(out_dir,"genomic_class.png"),p,width=9,height=5,dpi=160)
  }
}

# 11. EKSPORT GRUP BIOLOGICZNYCH ----------------------------------------------
# Zachowujemy oryginalne naglowki FASTA. Niezalezne od pozniejszego splitu.
if("reference" %in% names(seqs)) {
  cDNA <- seqs$reference; cdna_tx_ids <- fasta_meta$reference$transcript_id
  groups <- lapply(group_masks,function(m) cDNA[cdna_tx_ids %in% rf$transcript_id[m]])
  fasta_pc <- groups$protein_coding; fasta_noncoding <- groups$non_coding
  fasta_lncrna <- groups$lncRNA; fasta_ambiguous <- groups$ambiguous
  stopifnot(all(names(fasta_lncrna) %in% names(fasta_noncoding)),
            !length(intersect(names(fasta_pc),names(fasta_noncoding))),
            !length(intersect(names(fasta_ambiguous),c(names(fasta_pc),names(fasta_noncoding)))))
  if(write_outputs) for(g in names(groups)) write_fasta(groups[[g]],file.path(out_dir,paste0(g,".fa.gz")))
  missing_ids <- rf$transcript_id[rf$main_class=="annotation_missing"]
  if(write_outputs) write_fasta(cDNA[cdna_tx_ids %in% missing_ids],file.path(out_dir,"annotation_missing.fa.gz"))
}

# 12. DEDUPLIKACJA I KONTROLA JAKOSCI ------------------------------------------
# Oryginalna pelna tabela nadal zawiera WSZYSTKIE rekordy. Tylko do splitu
# przygotowujemy mniejsza tabele zawierajaca sekwencje oraz klucze do cech.
# Brak cechy GTF nie wyklucza automatycznie sekwencji: Core moze miec NA.
if(run_split) {
  raw <- list()
  for(src in names(seqs)) {
    ids <- fasta_meta[[src]]$transcript_id
    x <- final_features[match(paste(src,ids,sep="::"),final_features$record_id),]
    stopifnot(all(x$has_fasta),!anyNA(x$record_id))
    raw_group <- if(src=="experimental") rep("ex",nrow(x)) else ifelse(x$main_class=="protein_coding","cds",
        ifelse(x$is_lncrna %in% TRUE,"lnc_ext",ifelse(x$main_class=="non_coding","nc",x$main_class)))
    bad_ref <- src=="reference" & (x$label_conflict | x$version_match %in% FALSE |
       (exclude_length_mismatch & x$length_match %in% FALSE))
    raw_group[bad_ref] <- "annotation_missing"
    raw[[src]] <- data.frame(record_id=x$record_id,source=src,old_name=x$original_header,raw_group=raw_group,
                             is_lncrna=x$is_lncrna,seq_char=as.character(seqs[[src]]))
  }
  input_records <- do.call(rbind,unname(raw))
  dane_surowe <- load_and_dedup(df=input_records,extract_lnc=extract_lnc,dedup_mode=dedup_mode)
  dedup_audit <- attr(dane_surowe,"dedup_audit")
  inspect_counts(dane_surowe)
  # Zmiana A/B: zmien extract_lnc i uruchom ponownie sekcje 12-16.
}

# 13. PLANOWANIE LICZEBNOSCI --------------------------------------------------
# EPV/EPP = orientacyjny scenariusz dla klasy ograniczajacej, NIE prog decyzji.
# Riley et al. BMJ 2020, doi:10.1136/bmj.m441: nie ma jednej uniwersalnej reguly.
epv_plan <- estimate_sample_size(planned_ml_parameters,epv_values,cp_stage_prop,ml_tt_prop)
print(epv_plan)
# Przyblizona liczba niezaleznych pozytywnych w tescie dla CI czulosci;
# analogicznie osobno negatywne dla specyficznosci. Wald = plan, nie finalny CI.
test_n_approx <- ceiling(qnorm(.975)^2*test_target_sensitivity*(1-test_target_sensitivity)/test_ci_half_width^2)
cat("Orientacyjna liczebnosc jednej klasy testowej (Wald):",test_n_approx,"\n")
test_size_plan <- data.frame(anticipated_proportion=test_target_sensitivity,half_width=test_ci_half_width,
                             approximate_n_per_class=test_n_approx,method="Wald planning approximation")
# Krzywa uczenia wymaga faktycznego treningu; nie rysujemy fikcyjnego wyniku.
# Plan do eksperymentu w osobnym learning_curves.R: ponizsze liczebnosci i seed.
learning_curve_plan <- expand.grid(n_per_class=c(100,250,500,1000,2000,5000),repeat_seed=seed+0:4)
if(run_split) {
  limiting <- if(extract_lnc) min(sum(dane_surowe$class=="lnc"),sum(dane_surowe$class=="cds")) else
    min(sum(dane_surowe$class=="nc"),sum(dane_surowe$class=="cds"))
  cat("Pula ograniczajaca przed split:",limiting,"\n")
  cat("Orientacyjnie ML train przed limitowaniem:",round(limiting*(1-cp_stage_prop)*ml_tt_prop),"\n")
}

# 14. PODZIAL CP / ML / EX ----------------------------------------------------
if(run_split) {
  wynik <- split_datasets(dane_surowe,prefix,seed,cp_stage_prop,cp_tt_prop,ml_tt_prop,max_n_cp,max_n_ml,
                          balance_train,ml_neg_mode,extract_lnc)
  cp_train <- extract_subset(wynik,"cp","train"); cp_test <- extract_subset(wynik,"cp","test")
  ml_train <- extract_subset(wynik,"ml","train"); ml_test <- extract_subset(wynik,"ml","test")
  ex_pred <- extract_subset(wynik,"inference")
  sets <- list(cp_train=cp_train,cp_test=cp_test,ml_train=ml_train,ml_test=ml_test,ex_pred=ex_pred)
  df_ok <- wynik$legend[wynik$legend$keep,]
  print(table(paste(df_ok$stage,df_ok$tt_set,sep="_"),df_ok$role))
  print(table(paste(df_ok$stage,df_ok$tt_set,sep="_"),df_ok$class))
}

# 15. WALIDACJA I SELEKCJA CECH ML --------------------------------------------
# Jawna biala lista: metadane, nazwy genow, biotypy, odleglosc do PC i adnotowana
# dlugosc CDS nie trafia do domyslnej macierzy. Nie wykonujemy skalowania/imputacji
# przed podzialem; parametry preprocessing nalezy uczyc TYLKO na ml_train.
core_cols <- c("seq_length","gc_content","max_orf_len","orf_coverage","orf_count","sequence_entropy",
               "exon_count","mean_exon_len","median_exon_len","min_exon_len","max_exon_len","sd_exon_len",
               "first_exon_fraction","last_exon_fraction","intron_count","total_intron_len","mean_intron_len",
               "median_intron_len","min_intron_len","max_intron_len","sd_intron_len","genomic_span","exon_density",
               "pas_present","pas_distance_from_end","terminal_A_run")
# Ten sam zestaw cech Core w REF i EX. Gdy wylaczono genom, usuwamy ceche
# z panelu ML (pelna tabela nadal ma kolumny NA). Przy porownywaniu modeli
# z roznych przebiegow nalezy jawnie uzgodnic wspolny panel cech.
if(genome_engine!="none") core_cols <- c(core_cols,"noncanonical_fraction")
# Wariant rozszerzony jest WYLACZNIE kandydatem do eksperymentu: wlaczenie
# genomowego kontekstu wymaga niezaleznej referencji, bez self-match/leakage.
extended_candidates <- c("isoforms_per_gene","dist_to_pc","dist_to_nearest_tss","neighbor_orientation",
                         "overlap_same_strand","overlap_opposite_strand","exonic_overlap_fraction","overlapping_pc_genes")
ml_features_core <- final_features[,c("record_id",core_cols),drop=FALSE]
ml_features_extended <- final_features[,c("record_id",core_cols,extended_candidates),drop=FALSE]
# Zostawiamy oryginalne liczby w final_features. W panelach ML nie ufamy
# strukturze GTF przypisanej do innej wersji lub dlugosci sekwencji.
bad_structure <- final_features$length_match %in% FALSE | final_features$version_match %in% FALSE
struct_ml <- intersect(core_cols,c(struct_fields,splice_fields))
ml_features_core[bad_structure,struct_ml] <- NA
ml_features_extended[bad_structure,c(struct_ml,extended_candidates)] <- NA
# Slownik: kazda kolumna ma status; definicje podstawowych parametrow w MD.
feature_dictionary <- data.frame(feature=names(final_features),
  use=ifelse(names(final_features) %in% core_cols,"core_candidate",ifelse(names(final_features) %in% extended_candidates,"extended_requires_validation","description_or_qc")),
  stringsAsFactors=FALSE)
# Uzupelnienie jednostek i zrodla dla cech liczbowych oraz kategorii.
feature_dictionary$source <- ifelse(feature_dictionary$feature %in% c("transcript_id","record_id","source","has_fasta","has_gtf","annotation_status","fasta_files"),"integration",
  ifelse(feature_dictionary$feature %in% c("seq_length","gc_content","a_content","c_content","g_content","t_content","n_fraction","ambiguous_fraction","sequence_entropy","max_orf_len","orf_count","orf_start","orf_end","orf_coverage","pas_present","pas_distance_from_end","terminal_A_run") | grepl("^k[23]_",feature_dictionary$feature),"FASTA",
  ifelse(feature_dictionary$feature %in% splice_fields,"GTF+genome",ifelse(feature_dictionary$feature %in% c(genome_fields,"genomic_class","nearest_pc_gene_id","neighbor_orientation","overlap_pc","overlap_same_strand","overlap_opposite_strand"),"GTF+PC reference",
  ifelse(feature_dictionary$feature %in% c(struct_fields,meta_fields,"seqname","strand"),"GTF","integration/classification")))))
feature_dictionary$unit <- ifelse(grepl("fraction|coverage|density|_content$|^k[23]_",feature_dictionary$feature),"fraction_0_1",
  ifelse(grepl("length|_len$|_bp$|^dist_|^tx_start$|^tx_end$|^orf_start$|^orf_end$|genomic_span|terminal_A_run|length_difference",feature_dictionary$feature),"nt_or_bp",
  ifelse(grepl("count|junctions|^gt_ag$|^gc_ag$|^at_ac$|isoforms_per_gene|overlapping_pc_genes",feature_dictionary$feature),"count",ifelse(feature_dictionary$feature=="sequence_entropy","bits","category_or_ratio"))))
feature_dictionary$unit[feature_dictionary$feature=="pas_distance_from_end"] <- "nt"
# Jednoznaczne definicje najwazniejszych i potencjalnie mylacych kolumn.
definitions <- c(seq_length="FASTA sequence width",tx_length="Sum of annotated exon widths",
  gc_content="(G+C)/(A+C+G+T)",n_fraction="N / FASTA length",ambiguous_fraction="Non-ACGT / FASTA length",
  sequence_entropy="Shannon entropy of A,C,G,T frequencies (bits)",max_orf_len="Longest complete ATG-to-first-in-frame-stop ORF, including stop",
  orf_count="Number of complete start-to-stop ORFs; nested starts included",orf_start="Longest ORF start, 1-based cDNA",
  orf_end="Longest ORF end inclusive, 1-based cDNA",orf_coverage="max_orf_len / seq_length",
  pas_present="AATAAA or ATTAAA within terminal window",pas_distance_from_end="Bases after last PAS motif to sequence end",
  terminal_A_run="Consecutive terminal A bases; not proof of polyadenylation",exon_count="Number of exons in transcript",
  intron_count="Number of gaps between adjacent exons",total_intron_len="Sum of intron widths",genomic_span="Last genomic coordinate minus first plus one",
  exon_density="tx_length / genomic_span",exon_intron_ratio="tx_length / total_intron_len; NA without introns",
  first_exon_fraction="First exon in transcription direction / tx_length",last_exon_fraction="Last exon in transcription direction / tx_length",
  isoforms_per_gene="Number of transcript IDs in input annotation per gene_id",dist_to_pc="Gap in bases to nearest reference protein-coding gene span; adjacency also 0",
  dist_to_nearest_tss="Absolute coordinate difference to nearest reference PC transcript TSS",
  overlap_pc="True genomic span overlap with PC gene, irrespective of strand",exonic_overlap_bp="Union of intersecting query and PC exonic bases, either strand",
  exonic_overlap_fraction="exonic_overlap_bp / tx_length",overlapping_pc_genes="Number of overlapping reference PC gene spans",
  genomic_class="Descriptive position category; own PC locus may encode annotation label",
  canonical_junctions="GT-AG + GC-AG + AT-AC introns (project operational definition)",
  noncanonical_junctions="Other resolved splice motifs",unknown_junctions="Introns without resolved motif",
  noncanonical_fraction="noncanonical_junctions / resolved introns; NA if none",length_difference="seq_length - tx_length",
  fasta_files="Input FASTA path(s); identical merged IDs retain all paths separated by |; NA for GTF-only records",
  sequence_group="Index of first identical FASTA sequence within source; descriptive only")
feature_dictionary$definition <- unname(definitions[feature_dictionary$feature])
# Regular names have systematic meaning; provide a fallback to MD for metadata.
for(i in which(is.na(feature_dictionary$definition))) {
  nm <- feature_dictionary$feature[i]
  feature_dictionary$definition[i] <- if(grepl("^(mean|median|min|max|sd)_(exon|intron)_len$",nm)) gsub("_"," ",nm) else
    if(grepl("^k[23]_",nm)) paste("Frequency among valid windows:",sub("^k[23]_","",nm)) else paste("Metadata or QC field:",nm)
}
feature_dictionary$notes <- ifelse(feature_dictionary$use=="extended_requires_validation","Check self-match, annotation bias and identical computation for EX",ifelse(feature_dictionary$use=="core_candidate","Check missingness, collinearity and domain shift; fit preprocessing on training only","Not automatically an ML predictor; see companion Markdown"))
if(run_split) {
  validation <- validate_subsets(sets); print(validation)
  if(any(validation$shared_ids>0 | validation$shared_sequences>0)) stop("Data leakage: wspolne ID/sekwencje. Zmien dedup_mode lub dane.")
  # Kontrola rowniez wewnatrz zbiorow, w tym trybow none i header_id.
  if(any(vapply(sets,function(s) anyDuplicated(as.character(s))>0,logical(1)))) stop("Duplikaty wewnatrz zbioru: uzyj exact_seq.")
  stopifnot(!anyDuplicated(wynik$legend$ID),all(wynik$legend$record_id %in% final_features$record_id))
  # Cechy dla CP/ML/EX laczymy przez record_id, a wyniki CP pozniej przez ID.
  # Rola pos/neg jest ETYKIETA, nie cecha. ID i record_id tez nie sa cechami.
  dataset_features <- lapply(sets,function(s) {
    l <- wynik$legend[match(names(s),wynik$legend$ID),c("ID","record_id","role"),drop=FALSE]
    x <- ml_features_core[match(l$record_id,ml_features_core$record_id),core_cols,drop=FALSE]
    cbind(l,x,row.names=NULL)
  })
  # Gene-level leakage pozostaje diagnostyka, bo glowna analiza jest na poziomie transkryptu.
  genes <- final_features$gene_id[match(wynik$legend$record_id,final_features$record_id)]
  gs <- lapply(sets,function(s) unique(na.omit(genes[match(names(s),wynik$legend$ID)])))
  # EX gene_id ma inna przestrzen nazw, wiec ten raport tylko dla czterech referencyjnych kohort.
  gene_overlap <- do.call(rbind,lapply(combn(names(gs)[1:4],2,simplify=FALSE),function(p)
    data.frame(set1=p[1],set2=p[2],shared_gene_ids=length(intersect(gs[[p[1]]],gs[[p[2]]])))))
  print(gene_overlap)
  # Porownanie docelowych liczebnosci z faktycznymi, nie gwarancja statystyczna.
  ml_n <- table(factor(wynik$legend$role[wynik$legend$stage=="ml" & wynik$legend$tt_set=="train" & wynik$legend$keep],levels=c("pos","neg")))
  print(epv_plan); cat("Rzeczywista mniejsza klasa ML train:",min(ml_n),"\n")
}

# 16. EKSPORT I REPRODUKOWALNOSC -----------------------------------------------
# Pliki wyjsciowe nie zawieraja wynikow CP ani wytrenowanego ML. Te powstana
# w kolejnych modulach, na nieprzecinajacych sie pulach danych.
if(write_outputs) {
  write_tsv(final_features,file.path(out_dir,"transcript_features.tsv"))
  write_tsv(summary_table,file.path(out_dir,"transcriptome_summary.tsv"))
  if(!is.null(experimental_summary)) write_tsv(experimental_summary,file.path(out_dir,"experimental_summary.tsv"))
  if(!is.null(splice_summary)) write_tsv(splice_summary,file.path(out_dir,"splice_summary.tsv"))
  write_tsv(feature_dictionary,file.path(out_dir,"feature_dictionary.tsv"))
  write_tsv(ml_features_core,file.path(out_dir,"ml_features_core.tsv"))
  write_tsv(ml_features_extended,file.path(out_dir,"ml_features_extended_candidates.tsv"))
  write_tsv(epv_plan,file.path(out_dir,"epv_plan.tsv"))
  write_tsv(test_size_plan,file.path(out_dir,"test_size_plan.tsv"))
  write_tsv(learning_curve_plan,file.path(out_dir,"learning_curve_plan.tsv"))
  write_tsv(data.frame(feature=names(final_features),missing_n=vapply(final_features,function(x) sum(is.na(x)),integer(1))),file.path(out_dir,"missingness.tsv"))
  qc_cols <- c("record_id","fasta_files","annotation_status","has_fasta","has_gtf","length_match","length_difference","version_match","label_conflict","main_class")
  write_tsv(final_features[,qc_cols],file.path(out_dir,"annotation_qc.tsv"))
  if(run_split) {
    write_tsv(dedup_audit,file.path(out_dir,"dedup_audit.tsv"))
    write_tsv(wynik$legend,file.path(out_dir,"dataset_legend.tsv"))
    write_tsv(validation,file.path(out_dir,"subset_validation.tsv"))
    write_tsv(gene_overlap,file.path(out_dir,"gene_overlap.tsv"))
    for(nm in names(sets)) {
      write_fasta(sets[[nm]],file.path(out_dir,paste0(nm,".fa.gz")))
      write_tsv(dataset_features[[nm]],file.path(out_dir,paste0(nm,"_features.tsv")))
    }
    saveRDS(wynik,file.path(out_dir,"dataset_split.rds"))
  }
  saveRDS(final_features,file.path(out_dir,"transcript_features.rds"))
  config <- list(script_version="4.1_draft",input_mode=input_mode,ref_fasta=ref_fasta,ref_gtf=ref_gtf,pc_file=pc_file,nc_file=nc_file,lnc_file=lnc_file,
    ambiguous_file=ambiguous_file,separate_gtf=separate_gtf,separate_representation=separate_representation,include_ex=include_ex,
    ex_file=ex_file,ex_gtf=ex_gtf,genome_fasta=genome_fasta,genome_engine=genome_engine,
    splice_chunk_size=splice_chunk_size,assembly=assembly,reference_release=reference_release,
    ref_id_mode=ref_id_mode,ex_id_mode=ex_id_mode,structure_engine=structure_engine,calc_kmers=calc_kmers,pas_window=pas_window,
    min_orf_nt=min_orf_nt,sequence_chunk_size=sequence_chunk_size,sequence_chunk_nt=sequence_chunk_nt,progress_enabled=progress_enabled,save_sequence_checkpoint=save_sequence_checkpoint,run_split=run_split,extract_lnc=extract_lnc,dedup_mode=dedup_mode,exclude_length_mismatch=exclude_length_mismatch,seed=seed,prefix=prefix,
    cp_stage_prop=cp_stage_prop,cp_tt_prop=cp_tt_prop,ml_tt_prop=ml_tt_prop,max_n_cp=max_n_cp,max_n_ml=max_n_ml,
    balance_train=balance_train,ml_neg_mode=ml_neg_mode,planned_ml_parameters=planned_ml_parameters,epv_values=epv_values,
    test_target_sensitivity=test_target_sensitivity,test_ci_half_width=test_ci_half_width,
    out_dir=out_dir,write_outputs=write_outputs,make_plots=make_plots,run_self_tests=run_self_tests,
    self_tests_only=self_tests_only,self_test_rsamtools=self_test_rsamtools,save_workspace=save_workspace,
    pc_biotypes=pc_biotypes,lncrna_biotypes=lncrna_biotypes,all_ncrna_biotypes=all_ncrna_biotypes,core_cols=core_cols)
  dput(config,file=file.path(out_dir,"run_config.R"))
  active_files <- c(if(input_mode=="ensembl") c(ref_fasta,ref_gtf) else if(input_mode=="separate") c(pc_file,nc_file,lnc_file,ambiguous_file,separate_gtf) else ref_gtf,
                    if(include_ex || input_mode=="experimental_only") c(ex_file,ex_gtf),
                    if(genome_engine!="none") genome_fasta,
                    if(genome_engine=="rsamtools") c(paste0(genome_fasta,".fai"),paste0(genome_fasta,".gzi")))
  local <- unique(active_files[file.exists(active_files)])
  write_tsv(data.frame(path=local,md5=unname(tools::md5sum(local))),file.path(out_dir,"input_checksums.tsv"))
  capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
  if(save_workspace) save.image(file=file.path(out_dir,"dataset_v4.1_workspace.RData"))
}
message("Koniec przygotowania danych. Sprawdz raporty QC, brakujace cechy i liczebnosci przed treningiem.")


# ///////////////////////////////////////////////////////////////////////////////
# NOTATKI DO PRZYSZLEGO TUTORIALU V4.1 (przyklady, NIE kod wykonywany)
# ///////////////////////////////////////////////////////////////////////////////
# TEST 0: male testy syntetyczne bez danych. Ustaw self_tests_only <- TRUE
# w 00 i uruchom source("prepare_dataset_v4.1.R"). Kontrolowany stop() po OK
# jest celowy. Przy pracy sekcjami w RStudio nie wykonuj 01 po tym komunikacie.
# self_test_rsamtools=TRUE wlacza dodatkowa kontrole obu silnikow na malym FASTA.
# Po tescie ustaw self_tests_only <- FALSE, run_self_tests <- FALSE.
#
# TEST 1: nowy Ensembl, tylko baza. Ustaw ref_fasta/ref_gtf i wydanie/assembly,
# include_ex=FALSE, genome_engine="none" (odtworzenie przebiegu V3 bez DNA),
# run_split=TRUE, nowy out_dir. Wykonaj 00-16. Jesli masz genom z tego assembly,
# mozesz powtorzyc test z biostrings. Nie korzystamy z CDS-only do treningu cDNA.
#
# TEST 2: Mm115 REF+EX, aktywny przyklad w 00. Dwa FASTA: cdna.all + ncrna.
# Po 03A sprawdz fasta_gtf_coverage.tsv; oczekujemy 275134 ID w REF FASTA
# oraz 152133 dopasowanych lncRNA (przed QC/deduplikacja). Oryginalny genom .fa.gz,
# genome_engine="biostrings", include_ex=TRUE, run_split=TRUE. Zadbaj o RAM.
# Po 07 sprawdz splice_summary (reference + experimental), zwlaszcza unknown,
# resolved_percent oraz nazwy kontigow. Po 08 sprawdz zgodnosc ID/dlugosci.
# Po 15 sprawdz subset_validation i gene_overlap. EX ma role unknown.
#
# TEST 3: REF+EX bez genomu. genome_engine="none". Plik DNA nie jest potrzebny;
# pozostale cechy liczone jak w TEST 2, pelna tabela ma NA dla splicingu.
# W panelu Core nie ma noncanonical_fraction. Nie porownuj modeli z roznymi
# panelami bez jawnego uwzglednienia tej roznicy w metodologii.
#
# TEST 4 (opcjonalny, osobna sesja R): Rsamtools. Wymaga poprawnego odczytu
# duzego FASTA na danym systemie. Przyklad testu GRCm39 release 115 (nie dla
# innych assembly):
# p <- "Mm115/Mus_musculus.GRCm39.dna.primary_assembly.fa.bgz"
# q <- GRanges("7",IRanges::IRanges(100000000,100000020))
# stopifnot(as.character(Rsamtools::scanFa(p,param=q))=="TTAGAGCATATCTCTGAATAC")
# Nie uruchamiaj na Windows w sesji z niezapisanymi wynikami. Test malego FASTA
# nie gwarantuje prawidlowego dostepu do duzych plikow. Brak fallbacku po awarii.
#
# Zakres V4 roboczej: nie zmieniamy klasyfikacji biotypow, splitu, ORF/PAS,
# obliczen kontekstu PC, definicji motywow ani nie dodajemy nowych cech ML.
# Dalsza optymalizacja pamieci, unikalne junctions, porownawcze wykresy,
# walidacja gene-aware i pelny Markdown pozostaja do nastepnego etapu.
