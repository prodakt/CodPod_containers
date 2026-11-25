# ==============================================================================
# prepare_dataset_v2.R | CodPod_containers / lncRNA
# Charakterystyka transkryptow i przygotowanie zbiorow CP -> ML -> EX
# Dokumentacja rownolegla: prepare_dataset_v2.md (sekcje 00-16).
# Uruchamiaj kolejne sekcje w RStudio; ustaw sciezki w sekcji 00.
# Nie instaluje pakietow, nie pobiera danych i nie zapisuje nic na GitHubie.
# Przyklad korzysta z pelnego cDNA, nie z wycietych sekwencji CDS.
# ==============================================================================

# 00. PAKIETY I KONFIGURACJA ----------------------------------------------------
# Instalacja jednorazowa (odkomentuj w razie potrzeby):
# install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "rtracklayer", "GenomicRanges",
#                        "GenomicFeatures", "txdbmaker", "Rsamtools"))
# install.packages("ggplot2")

# Ustawienia wejscia: "ensembl", "separate" lub "experimental_only".
input_mode <- "ensembl"
ref_fasta <- "cdna.fa.bgz"
ref_gtf <- "genes.gtf.bgz"           # mozna zmienic na genes.gtf.gz
ref_id_mode <- "ensembl"             # "exact" dla innych baz; patrz 02
reference_release <- "Ensembl geneset 2026_04"  # zmien zgodnie z plikami
assembly <- "GRCm39"
prefix <- "Mus"

# Alternatywa: osobne FASTA. Najlepiej pelne transkrypty kodujace, NIE CDS.
pc_file <- "protein_coding.fa.gz"
nc_file <- "non_coding.fa.gz"
lnc_file <- NULL                     # opcjonalny, moze byc podzbiorem nc_file
ambiguous_file <- NULL               # opcjonalny
separate_gtf <- NULL                 # np. "genes.gtf.gz"; NULL = cechy FASTA
separate_representation <- "full_transcript" # lub "CDS_only" (nie do pelnego ML)

# EX mozna dolaczyc teraz lub pozniej. Wylaczone w domyslnym przykladzie.
include_ex <- FALSE
ex_file <- "merged_transcriptome.fasta"
ex_gtf <- "merged.gtf"
ex_id_mode <- "exact"                # NIE usuwamy kropek z MSTRG.123.1
# Przy experimental_only mozna podac ref_gtf jako niezalezna adnotacje PC.
# Dla braku takiej referencji ustaw ref_gtf <- NULL.

# GTF bezposrednio = zachowuje metadane i obsluguje GTF bez wierszy transcript.
# "txdb" = alternatywa zgodna z poprzednim kodem, kontroluje zgodnosc eksonow.
structure_engine <- "gtf"            # "gtf" lub "txdb"
genome_fasta <- NULL                  # opcjonalny lokalny, indeksowany genome.fa
# FASTA genomu i GTF musza miec zgodny assembly, nazwy chromosomow i orientacje.
calc_kmers <- FALSE                   # 16 dinukleotydow zawsze opisowo; TRUE dodaje 64 trimery
pas_window <- 100L                    # ostatnie 100 nt; motywy w sekwencji DNA
min_orf_nt <- 6L                     # ATG + STOP, wliczamy STOP

# Ustawienia do pozniejszego podzialu (sekcje 12-14).
run_split <- TRUE                    # FALSE: tylko charakterystyka i eksport
extract_lnc <- TRUE                  # TRUE = A, FALSE = B; decyzja reczna!
dedup_mode <- "exact_seq"            # "exact_seq", "header_id", "none"
exclude_length_mismatch <- TRUE       # wyklucz niezgodne FASTA/GTF z treningu ref
seed <- 123L
cp_stage_prop <- 0.4
cp_tt_prop <- 0.8
ml_tt_prop <- 0.75
max_n_cp <- 2000                     # limit NA ROLE, train + test lacznie
max_n_ml <- 3000
balance_train <- TRUE
ml_neg_mode <- "cds_only"            # lub "cds_and_nc_other" (A)
# EPV jest jedynie heurystyka, nie gwarancja mocy. Uwzglednij przyszle wyniki CP.
planned_ml_parameters <- 25L
epv_values <- c(10, 15, 20)
test_target_sensitivity <- 0.90
test_ci_half_width <- 0.03

out_dir <- "dataset_v2_output"       # najlepiej nowy katalog dla kazdego przebiegu
write_outputs <- TRUE
make_plots <- TRUE                   # ggplot2; wykresy opisowe
run_self_tests <- FALSE              # TRUE uruchamia male testy przed importem

# Podstawowe pakiety; stosujemy jawne namespace, by uniknac konfliktu select().
required <- c("Biostrings", "rtracklayer", "GenomicRanges", "IRanges", "S4Vectors")
if (structure_engine == "txdb") required <- c(required, "GenomicFeatures", "txdbmaker")
if (!is.null(genome_fasta)) required <- c(required, "Rsamtools")
if (make_plots) required <- c(required, "ggplot2")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Brak pakietow: ", paste(missing_packages, collapse = ", "))
# S4 metody width/start/end/strand itd. dostepne po zaladowaniu pakietow.
library(Biostrings)
library(GenomicRanges)
stopifnot(input_mode %in% c("ensembl", "separate", "experimental_only"),
          structure_engine %in% c("gtf", "txdb"),
          separate_representation %in% c("full_transcript", "CDS_only"))
if (input_mode == "separate" && separate_representation == "CDS_only" && run_split)
  stop("CDS_only: ustaw run_split=FALSE lub dostarcz pelne transkrypty kodujace.")
if (input_mode == "experimental_only") run_split <- FALSE
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

# ==============================================================================
# MALE FUNKCJE POMOCNICZE (definicje przed kodem glownym 01-16)
# ==============================================================================

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
  if (anyDuplicated(unique(data.frame(id=ids, raw=raw))$id))
    stop("Rozne wersje ID GTF po normalizacji. Sprawdz release i identyfikatory.")
  # Pierwszy niepusty atrybut na ID, ze sprawdzeniem konfliktow.
  fields <- c("gene_id", "gene_name", "gene_biotype", "transcript_biotype", "transcript_version")
  vals <- lapply(fields, function(f) gtf_field(q, f)); names(vals) <- fields
  # Czesci baz uzywaja gene_type i transcript_type zamiast *_biotype.
  for (pair in list(c("gene_biotype", "gene_type"), c("transcript_biotype", "transcript_type"))) {
    idx <- is.na(vals[[pair[1]]]); vals[[pair[1]]][idx] <- gtf_field(q, pair[2])[idx]
  }
  # split raz; nie uzywamy wielokrotnego filtrowania calego GTF na kazdy ID.
  grp <- split(seq_along(ids), ids)
  pick <- function(x, field) vapply(grp, function(i) {
    z <- unique(x[i][!is.na(x[i]) & nzchar(x[i])])
    if (length(z) > 1L) stop("Konflikt ", field, " dla ", ids[i[1]])
    if (length(z)) z[1] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  out <- data.frame(transcript_id = names(grp),
                    gtf_transcript_id = pick(raw, "transcript_id"), stringsAsFactors = FALSE)
  for (f in fields) out[[f]] <- pick(vals[[f]], f)
  out
}

# 03: Eksony bezposrednio z GTF. Usuwamy tylko powtorzone zakresy tego samego
# transkryptu, nie eksony wspolne dla roznych izoform.
get_exons_by_tx <- function(gr, mode, engine = "gtf") {
  if (engine == "txdb") {
    txdb <- txdbmaker::makeTxDbFromGRanges(gr)
    ex <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
    names(ex) <- extract_transcript_ids(names(ex), mode)
    if (anyNA(names(ex)) || anyDuplicated(names(ex))) stop("TxDb: nieunikalne ID")
    return(ex)
  }
  e <- gr[which(gtf_field(gr, "type") == "exon")]
  ids <- extract_transcript_ids(gtf_field(e, "transcript_id"), mode)
  ok <- !is.na(ids) & nzchar(ids)
  e <- e[ok]; ids <- ids[ok]
  if (!length(e)) stop("Brak eksonow z transcript_id w GTF.")
  S4Vectors::mcols(e) <- NULL
  e <- split(e, ids)
  # Uklad 5' -> 3' rowniez na nici minus; nieznana nic pozostaje oznaczona '*'.
  endoapply(e, function(z) {
    z <- unique(z)
    if (length(unique(as.character(seqnames(z)))) != 1L || length(unique(as.character(strand(z)))) != 1L)
      stop("Transkrypt obejmuje rozne chromosomy lub nici.")
    z[order(start(z), decreasing = as.character(strand(z)[1]) == "-")]
  })
}

# 04: Najdluzszy kompletny ORF: ATG do PIERWSZEGO STOP w tej samej ramce.
# Wszystkie starty liczone osobno (ORF zagniezdzone tez); 3 ramki na nici cDNA.
# Standardowy kod jadrowy. Nie jest to dowod translacji; nie liczymy ORF bez STOP.
calc_longest_orf_bp <- function(dna_set, min_nt = 6L) {
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

# 04: Cechy FASTA obliczane niezaleznie od etykiety. Sklad jako ulamek 0-1.
# GC = (G+C)/(A+C+G+T); niejednoznaczne symbole raportujemy osobno.
calc_sequence_features <- function(s, ids, pas_window=100L, min_orf_nt=6L, kmers=FALSE) {
  n <- width(s); acgt <- letterFrequency(s, letters=c("A","C","G","T"))
  denom <- rowSums(acgt); p <- acgt / ifelse(denom == 0, NA_real_, denom)
  ent <- -rowSums(ifelse(p > 0, p * log2(p), 0))
  ent[denom == 0] <- NA_real_
  out <- data.frame(transcript_id=ids, seq_length=n, gc_content=p[,"G"]+p[,"C"],
                    a_content=p[,"A"], c_content=p[,"C"], g_content=p[,"G"], t_content=p[,"T"],
                    n_fraction=as.numeric(letterFrequency(s, "N"))/n,
                    ambiguous_fraction=1-denom/n, sequence_entropy=ent)
  orf <- calc_longest_orf_bp(s, min_orf_nt)
  out <- cbind(out, orf); out$orf_coverage <- out$max_orf_len/n
  # PAS: tylko dwa klasyczne motywy DNA AATAAA/ATTAAA w koncowym oknie.
  # Odleglosc liczona od KONCA motywu do konca sekwencji. Brak motywu = NA.
  tail_seq <- as.character(subseq(s, start=pmax(1L, n-pas_window+1L)))
  pas <- vapply(tail_seq, function(x) {
    # Lookahead pozwala znalezc rowniez nakladajace sie motywy.
    hits <- gregexpr("(?=AATAAA|ATTAAA)", x, perl=TRUE)[[1]]
    if (hits[1] < 0) NA_real_ else nchar(x)-max(hits)-5L
  }, numeric(1))
  out$pas_present <- !is.na(pas); out$pas_distance_from_end <- pas
  # Terminalny ciag A NIE jest dowodem obecnosci ogona poly(A) in vivo.
  tail_a <- regexpr("A+$", as.character(s), perl=TRUE)
  out$terminal_A_run <- pmax(0L, attr(tail_a, "match.length"))
  # Częstosc dinukleotydow na wsrod prawidlowych okien (bez niejednoznacznych).
  for (k in if (kmers) c(2L,3L) else 2L) {
    freq <- oligonucleotideFrequency(s, width=k, as.prob=FALSE)
    freq <- freq/ifelse(rowSums(freq)==0, NA_real_, rowSums(freq))
    colnames(freq) <- paste0("k",k,"_",colnames(freq))
    out <- cbind(out, as.data.frame(freq))
  }
  out
}

# 05: Statystyki eksonow i intronow z rzeczywistych zakresow.
# Brak intronow: liczba i suma = 0, srednia/min/max/SD/ratio = NA (nie dotyczy).
# Nakladajace sie eksony w jednym transkrypcie sa bledem struktury, nie sumujemy
# ich po cichu. Nieznana nic: nie interpretujemy pierwszego/ostatniego eksonu.
calc_exon_features <- function(ex) {
  rows <- lapply(seq_along(ex), function(i) {
    z <- ex[[i]]; w <- width(z); o <- order(start(z)); z <- z[o]
    gap <- if (length(z)>1L) start(z)[-1L]-end(z)[-length(z)]-1L else integer()
    if (any(gap <= 0L)) stop("Nakladajace lub przylegajace eksony: ", names(ex)[i])
    span <- max(end(z))-min(start(z))+1L
    minus <- as.character(strand(z)[1]) == "-"; known <- as.character(strand(z)[1]) != "*"
    first <- if (known) if (minus) tail(w[o],1) else w[o][1] else NA_real_
    last <- if (known) if (minus) w[o][1] else tail(w[o],1) else NA_real_
    c(exon_count=length(w), tx_length=sum(w), mean_exon_len=mean(w), median_exon_len=median(w),
      min_exon_len=min(w), max_exon_len=max(w), sd_exon_len=if(length(w)>1) sd(w) else NA_real_,
      first_exon_len=first, last_exon_len=last, genomic_span=span, exon_density=sum(w)/span,
      first_exon_fraction=first/sum(w), last_exon_fraction=last/sum(w),
      intron_count=length(gap), total_intron_len=sum(gap),
      mean_intron_len=if(length(gap)) mean(gap) else NA_real_,
      median_intron_len=if(length(gap)) median(gap) else NA_real_,
      min_intron_len=if(length(gap)) min(gap) else NA_real_,
      max_intron_len=if(length(gap)) max(gap) else NA_real_,
      sd_intron_len=if(length(gap)>1) sd(gap) else NA_real_,
      exon_intron_ratio=if(length(gap)) sum(w)/sum(gap) else NA_real_)
  })
  data.frame(transcript_id=names(ex), as.data.frame(do.call(rbind, rows)), row.names=NULL)
}

# 05: Intron jako przerwa miedzy kolejnymi eksonami tej samej izoformy.
get_introns_by_tx <- function(ex) {
  endoapply(ex, function(z) {
    z <- z[order(start(z))]; n <- length(z)
    if (n < 2L) return(z[FALSE])
    a <- end(z)[-n]+1L; b <- start(z)[-1L]-1L
    if (any(a>b)) stop("Niepoprawny intron.")
    GRanges(seqnames=as.character(seqnames(z)[-n]), IRanges::IRanges(a,b),
            strand=as.character(strand(z)[-n]))
  })
}

# 05: Zakresy transkryptow od pierwszej do ostatniej pozycji genomowej.
get_tx_ranges <- function(ex) {
  z <- unlist(range(ex), use.names=FALSE)
  if (length(z)!=length(ex)) stop("Niepoprawny zakres transkryptu.")
  names(z) <- names(ex); z
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
  gene_tx <- split(tx[idx], gene[idx]); gene_ranges <- unlist(range(gene_tx), use.names=FALSE)
  if (length(gene_ranges)!=length(gene_tx)) stop("Gen PC na wielu kontigach/ niciach; sprawdz adnotacje.")
  names(gene_ranges) <- names(gene_tx)
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
    inter <- IRanges::IRanges(pmax(start(qe)[eq],start(pc_ref$exons)[es]),
                             pmin(end(qe)[eq],end(pc_ref$exons)[es]))
    by_tx <- split(inter, owner[eq]); bp <- vapply(by_tx,function(z) sum(width(reduce(z))),numeric(1))
    out$exonic_overlap_bp[as.integer(names(bp))] <- bp
  }
  txlen <- vapply(ex[ids],function(z) sum(width(z)),numeric(1))
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

# 07: Splice motifs z genomu. Pobieramy tylko po 2 nt z obu koncow intronu.
# Dla nici minus zamieniamy donor/akceptor i stosujemy reverse complement.
# GT-AG (major) oraz GC-AG i AT-AC (rzadsze) raportujemy oddzielnie.
# Nieznana nic, intron <4 nt, brak kontigu lub niejednoznaczne zasady = unknown.
calc_splice_features <- function(introns, genome_file) {
  n <- length(introns); counts <- lengths(introns)
  out <- data.frame(transcript_id=names(introns), gt_ag=rep(NA_integer_,n), gc_ag=rep(NA_integer_,n),
    at_ac=rep(NA_integer_,n), canonical_junctions=rep(NA_integer_,n), noncanonical_junctions=rep(NA_integer_,n),
    unknown_junctions=rep(NA_integer_,n), noncanonical_fraction=rep(NA_real_,n))
  if(is.null(genome_file)) return(out)
  if (!file.exists(genome_file)) stop("Brak genomu: ",genome_file)
  if (!file.exists(paste0(genome_file,".fai"))) stop("Najpierw Rsamtools::indexFa(genome_fasta).")
  fa <- Rsamtools::FaFile(genome_file); idx <- Rsamtools::scanFaIndex(fa)
  lens <- setNames(width(idx), as.character(seqnames(idx)))
  z <- unlist(introns,use.names=FALSE); owner <- rep(seq_len(n),counts)
  motif <- rep("unknown",length(z)); chr <- as.character(seqnames(z)); st <- as.character(strand(z))
  valid <- which(st %in% c("+","-") & width(z)>=4L & chr %in% names(lens) &
                 start(z)>=1L & end(z)<=lens[chr])
  if(length(valid)) {
    v <- z[valid]
    left <- GRanges(as.character(seqnames(v)), IRanges::IRanges(start(v),start(v)+1L))
    right <- GRanges(as.character(seqnames(v)), IRanges::IRanges(end(v)-1L,end(v)))
    a <- Rsamtools::scanFa(fa,param=left); b <- Rsamtools::scanFa(fa,param=right)
    minus <- st[valid]=="-"; d <- as.character(a); ac <- as.character(b)
    d[minus] <- as.character(reverseComplement(b[minus])); ac[minus] <- as.character(reverseComplement(a[minus]))
    good <- grepl("^[ACGT]{2}$",d) & grepl("^[ACGT]{2}$",ac)
    motif[valid[good]] <- paste(d[good],ac[good],sep="-")
  }
  count <- function(m) tabulate(owner[m],nbins=n)
  out$gt_ag <- count(motif=="GT-AG"); out$gc_ag <- count(motif=="GC-AG"); out$at_ac <- count(motif=="AT-AC")
  out$canonical_junctions <- out$gt_ag+out$gc_ag+out$at_ac
  out$unknown_junctions <- count(motif=="unknown")
  out$noncanonical_junctions <- counts-out$canonical_junctions-out$unknown_junctions
  tested <- counts-out$unknown_junctions
  out$noncanonical_fraction <- ifelse(tested>0,out$noncanonical_junctions/tested,NA_real_)
  out
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
  message("Self-tests OK (ID, ORF, struktura, A/B, deduplikacja i rozlacznosc).")
}

# ==============================================================================
# KOD GLOWNY — wykonuj kolejne sekcje w RStudio
# ==============================================================================

# 01. DANE WEJSCIOWE -----------------------------------------------------------
# Nowy Ensembl: przypiety przyklad z plikow uzytkownika. Pobieranie opcjonalne:
# base_url <- "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/635/9/ensembl/2026_04/geneset/"
# options(timeout=1200)
# download.file(paste0(base_url,"cdna.fa.bgz"), "cdna.fa.bgz", mode="wb")
# download.file(paste0(base_url,"genes.gtf.bgz"), "genes.gtf.bgz", mode="wb")
# Nie uzywaj niepelnego adresu https://ensembl.org jako adresu FASTA/GTF.

seqs <- list(); group_hint <- list(); anno <- list(); id_mode <- list()
if(input_mode=="ensembl") {
  seqs$reference <- read_fasta(ref_fasta); group_hint$reference <- rep(NA_character_,length(seqs$reference))
  anno$reference <- rtracklayer::import(ref_gtf,format="gtf"); id_mode$reference <- ref_id_mode
} else if(input_mode=="separate") {
  # LncRNA moze byc celowo obecne w obu plikach NC i lnc. Po identyfikatorach
  # scalimy te rekordy w 02, a nie potraktujemy jako dwa rozne transkrypty.
  files <- list(protein_coding=pc_file,non_coding=nc_file,lncRNA=lnc_file,ambiguous=ambiguous_file)
  ss <- lapply(files,read_fasta)
  seqs$reference <- do.call(c,unname(ss)); group_hint$reference <- rep(names(ss),lengths(ss))
  if(!is.null(separate_gtf)) anno$reference <- rtracklayer::import(separate_gtf,format="gtf")
  id_mode$reference <- ref_id_mode
} else if(!is.null(ref_gtf)) {
  # Tylko EX: referencyjny GTF sluzy jako tlo genomowe, nie do nadawania klas EX.
  anno$reference <- rtracklayer::import(ref_gtf,format="gtf"); id_mode$reference <- ref_id_mode
}
if(include_ex || input_mode=="experimental_only") {
  seqs$experimental <- read_fasta(ex_file); group_hint$experimental <- rep(NA_character_,length(seqs$experimental))
  anno$experimental <- rtracklayer::import(ex_gtf,format="gtf"); id_mode$experimental <- ex_id_mode
}
print(vapply(seqs,length,integer(1)))

# 02. IDENTYFIKATORY ----------------------------------------------------------
# Nie zmieniamy oryginalnych naglowkow FASTA. record_id jest kluczem globalnym
# (source + transcript_id), a transcript_id jest kluczem wewnatrz danego zrodla.
fasta_meta <- list()
for(src in names(seqs)) {
  s <- seqs[[src]]; h <- names(s); ids <- extract_transcript_ids(h,id_mode[[src]])
  if(anyNA(ids)||any(!nzchar(ids))) stop("Brak ID FASTA: ",src)
  hints <- group_hint[[src]]
  # Zachowanie wersji/oryginalnego ID. Nie gubimy naglowka w eksportach biologicznych.
  if(anyDuplicated(ids)) {
    if(input_mode!="separate" || src!="reference") stop("Nieunikalne ID FASTA po normalizacji: ",src)
    groups <- split(which(duplicated(ids)|duplicated(ids,fromLast=TRUE)), ids[duplicated(ids)|duplicated(ids,fromLast=TRUE)])
    for(i in groups) {
      if(length(unique(as.character(s[i])))!=1L) stop("Ten sam ID, rozne sekwencje: ",ids[i[1]])
      # Tolerujemy tylko identyczne NC/lnc i duplikat tej samej klasy.
      if(length(unique(hints[i]))>1L && !all(hints[i] %in% c("non_coding","lncRNA"))) stop("Konflikt plikow dla ID: ",ids[i[1]])
      if(any(hints[i]=="lncRNA")) hints[i] <- "lncRNA"
    }
    keep <- !duplicated(ids); s <- s[keep]; h <- h[keep]; ids <- ids[keep]; hints <- hints[keep]
  }
  seqs[[src]] <- s
  fasta_meta[[src]] <- data.frame(transcript_id=ids,fasta_id=extract_transcript_ids(h,"exact"),original_header=h,
                                  group_hint=hints,has_fasta=TRUE,stringsAsFactors=FALSE)
  print(head(fasta_meta[[src]],2))
}

# 03. ADNOTACJA I METADANE -----------------------------------------------------
# Ten sam kod dla Ensembl i StringTie. GFF3 ma ID/Parent; najpierw przekonwertuj
# na GTF (np. gffread -T -o genes.gtf genes.gff3) i sprawdz zachowanie biotypow.
tx_meta <- list(); exons_by_tx <- list(); tx_ranges <- list(); introns_by_tx <- list()
for(src in names(anno)) {
  tx_meta[[src]] <- get_tx_meta(anno[[src]],id_mode[[src]])
  exons_by_tx[[src]] <- get_exons_by_tx(anno[[src]],id_mode[[src]],structure_engine)
  tx_ranges[[src]] <- get_tx_ranges(exons_by_tx[[src]])
  introns_by_tx[[src]] <- get_introns_by_tx(exons_by_tx[[src]])
  # Kontrola TxDb / GTF: nie pozwalamy na ciche zgubienie eksonow przez adapter.
  if(!all(names(exons_by_tx[[src]]) %in% tx_meta[[src]]$transcript_id)) stop("Eksony bez metadanych: ",src)
  print(head(tx_meta[[src]]))
}

# 04. CECHY SEKWENCYJNE -------------------------------------------------------
seq_features <- list()
for(src in names(seqs)) {
  seq_features[[src]] <- calc_sequence_features(seqs[[src]],fasta_meta[[src]]$transcript_id,pas_window,min_orf_nt,calc_kmers)
  print(head(seq_features[[src]][,1:8]))
}

# 05. ARCHITEKTURA EKSON-INTRON ------------------------------------------------
struct_features <- list()
for(src in names(exons_by_tx)) {
  e <- exons_by_tx[[src]]; tx <- tx_ranges[[src]]
  x <- calc_exon_features(e)
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
if("reference" %in% names(tx_meta)) pc_ref <- build_pc_reference(tx_meta$reference,exons_by_tx$reference,tx_ranges$reference,introns_by_tx$reference)
if(is.null(pc_ref)) message("Brak wiarygodnej referencji PC: kontekst genomowy = NA.")
genomic_features <- list()
for(src in names(tx_ranges)) {
  genomic_features[[src]] <- calc_genomic_features(tx_ranges[[src]],exons_by_tx[[src]],tx_meta[[src]],pc_ref,src)
  print(table(genomic_features[[src]]$genomic_class,useNA="ifany"))
}

# 07. OPCJONALNE CECHY --------------------------------------------------------
# Do analizy motywow intronowych potrzebujemy genomu, NIE FASTA cDNA.
# Rsamtools::indexFa(genome_fasta) # raz dla lokalnego, nieskompresowanego lub BGZF FASTA
# Opcjonalne dalsze moduly: repeats, konserwacja, BAM/junction reads, ekspresja,
# UTR/CDS z adnotacji (tylko opis), MFE. Nie udajemy, ze zostaly juz obliczone.
splice_features <- list()
for(src in names(introns_by_tx)) splice_features[[src]] <- calc_splice_features(introns_by_tx[[src]],genome_fasta)

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
if(!is.null(genome_fasta)) core_cols <- c(core_cols,"noncanonical_fraction")
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
feature_dictionary$source <- ifelse(feature_dictionary$feature %in% c("transcript_id","record_id","source","has_fasta","has_gtf","annotation_status"),"integration",
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
  write_tsv(feature_dictionary,file.path(out_dir,"feature_dictionary.tsv"))
  write_tsv(ml_features_core,file.path(out_dir,"ml_features_core.tsv"))
  write_tsv(ml_features_extended,file.path(out_dir,"ml_features_extended_candidates.tsv"))
  write_tsv(epv_plan,file.path(out_dir,"epv_plan.tsv"))
  write_tsv(test_size_plan,file.path(out_dir,"test_size_plan.tsv"))
  write_tsv(learning_curve_plan,file.path(out_dir,"learning_curve_plan.tsv"))
  write_tsv(data.frame(feature=names(final_features),missing_n=vapply(final_features,function(x) sum(is.na(x)),integer(1))),file.path(out_dir,"missingness.tsv"))
  qc_cols <- c("record_id","annotation_status","has_fasta","has_gtf","length_match","length_difference","version_match","label_conflict","main_class")
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
  config <- list(input_mode=input_mode,ref_fasta=ref_fasta,ref_gtf=ref_gtf,pc_file=pc_file,nc_file=nc_file,lnc_file=lnc_file,
    ambiguous_file=ambiguous_file,separate_gtf=separate_gtf,separate_representation=separate_representation,include_ex=include_ex,
    ex_file=ex_file,ex_gtf=ex_gtf,genome_fasta=genome_fasta,assembly=assembly,reference_release=reference_release,
    ref_id_mode=ref_id_mode,ex_id_mode=ex_id_mode,structure_engine=structure_engine,calc_kmers=calc_kmers,pas_window=pas_window,
    min_orf_nt=min_orf_nt,run_split=run_split,extract_lnc=extract_lnc,dedup_mode=dedup_mode,exclude_length_mismatch=exclude_length_mismatch,seed=seed,prefix=prefix,
    cp_stage_prop=cp_stage_prop,cp_tt_prop=cp_tt_prop,ml_tt_prop=ml_tt_prop,max_n_cp=max_n_cp,max_n_ml=max_n_ml,
    balance_train=balance_train,ml_neg_mode=ml_neg_mode,planned_ml_parameters=planned_ml_parameters,epv_values=epv_values,
    pc_biotypes=pc_biotypes,lncrna_biotypes=lncrna_biotypes,all_ncrna_biotypes=all_ncrna_biotypes,core_cols=core_cols)
  dput(config,file=file.path(out_dir,"run_config.R"))
  active_files <- c(if(input_mode=="ensembl") c(ref_fasta,ref_gtf) else if(input_mode=="separate") c(pc_file,nc_file,lnc_file,ambiguous_file,separate_gtf) else ref_gtf,
                    if(include_ex || input_mode=="experimental_only") c(ex_file,ex_gtf),genome_fasta)
  local <- unique(active_files[file.exists(active_files)])
  write_tsv(data.frame(path=local,md5=unname(tools::md5sum(local))),file.path(out_dir,"input_checksums.tsv"))
  capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
}
message("Koniec przygotowania danych. Sprawdz raporty QC, brakujace cechy i liczebnosci przed treningiem.")
