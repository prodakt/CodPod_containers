# TODO — przygotowanie zbiorów danych i identyfikacja lncRNA (`CP → ML → EX`)

> **Notatnik roboczy / backlog projektu**  
> Aktualizacja: 2026-11-22  
> Repozytorium: [CodPod_containers / scripts](https://github.com/prodakt/CodPod_containers/tree/main/scripts)  
> **Zasada:** ten plik jest przygotowany lokalnie. Bez wyraźnego polecenia **nie zmieniać niczego w repozytorium GitHub**.

## Jak korzystać z notatnika

- `[ ]` — do zrobienia; `[x]` — zakończone.
- **P0** — poprawność biologiczna/metodologiczna, możliwy błąd lub *data leakage*.
- **P1** — ważne eksperymenty i rozwój metody.
- **P2** — dokumentacja, wygoda użycia i rozwój repozytorium.
- Przy zamykaniu zadania dopisać **wynik, decyzję, datę** i ewentualnie odsyłacz do testu/publikacji/commitu.

## 1. Poprawność przygotowania danych — P0

### Referencje i reprezentacja transkryptów

- [ ] **P0** Rozstrzygnąć, czy klasa kodująca (`cds_file`) powinna zawierać wyłącznie CDS, czy pełne transkrypty kodujące (cDNA z UTR). **Problem:** eksperymentalne transkrypty StringTie (`ex_file`) są zwykle pełnymi transkryptami, a sekwencja CDS odpowiada jedynie ich fragmentowi; deduplikacja `exact_seq` nie wykryje takiej zgodności.
  - [ ] Porównać warianty CDS vs pełne cDNA pod względem uczciwości uczenia CP, długości sekwencji i usuwania znanych transkryptów z `ex`.
  - [ ] Oddzielić ewentualnie referencję do **treningu CP** od referencji do **identyfikacji znanych transkryptów w EX**.
- [ ] **P0** Sprawdzić, jak `ex` jest odfiltrowywany względem referencji: identyczne sekwencje, różne długości/UTR, alternatywne izoformy, fragmenty, różne wersje adnotacji.
- [ ] **P0** Zdefiniować politykę rekordów konfliktowych: ta sama sekwencja/ID pojawia się jako `cds`, `lnc`, `nc_other` lub `ex`. Obecny priorytet deduplikacji to `cds > lnc_ext > nc > ex`; ocenić jego konsekwencje biologiczne.
- [ ] **P0** Zdecydować, jak traktować identyczne i niemal identyczne izoformy tego samego genu. Podstawową jednostką analizy pozostaje **transkrypt**; grupowanie po genie rozważyć jako dodatkowy, bardziej restrykcyjny benchmark (nie jako automatyczne zastąpienie obecnej metody).

### Funkcje R: zabezpieczenia

- [ ] **P0** W `split_datasets()` wybierać Ścieżkę A/B jawnie (obecnie skrypt wnioskuje po `"lnc" %in% df_ref$class`). Gdy wybrano Ścieżkę A bez wykrytych `lnc`, zwrócić czytelny błąd zamiast cichej zmiany logiki.
- [ ] **P0** Walidować `ml_neg_mode` (`"cds_only"` albo `"cds_and_nc_other"`) i pozostałe parametry: proporcje `0–1`, limity, brakujące pliki, puste klasy, niepoprawne etykiety.
- [ ] **P0** Weryfikować, że po limicie i balansowaniu każda wymagana kohorta ma próbki obu klas (w szczególności `ml_train`, `ml_test`, `cp_train`, `cp_test`).
- [ ] **P0** Zdefiniować dopuszczalne użycie `dedup_mode = "none"` i `"header_id"`; udokumentować ryzyko pozostawienia identycznych sekwencji pod różnymi identyfikatorami.
- [ ] **P0** Sprawdzić efekty `max_n_cp`, `max_n_ml` oraz `balance_train` na skład podklas, zwłaszcza udział `lnc` i `nc_other` we wspólnej klasie pozytywnej CP.
- [ ] **P0** Dodać kompleksową kontrolę rozłączności **wszystkich pięciu zbiorów**: `cp_train`, `cp_test`, `ml_train`, `ml_test`, `ex_pred` (identyfikatory, identyczne sekwencje; później opcjonalnie podobieństwo sekwencji i wspólne geny).
- [ ] **P0** Zapewnić, że do trenowania ML trafiają predykcje CP wykonane na danych **niewidzianych przez modele CP podczas treningu**, a `ex_pred` jest używane tylko do inferencji.
- [ ] **P0** Oddzielić zbiór służący do wyboru cech/hiperparametrów od ostatecznego testu; nie stroić modeli na zbiorze raportowanym jako niezależny test.

**Warunek zakończenia sekcji:** testy na małych sztucznych FASTA obejmują Ścieżkę A i B, duplikaty, braki klas, limity i brak przenikania rekordów między zbiorami.

## 2. Jak duże zbiory są wystarczające? — P1 / eksperyment metodologiczny

**Cel:** wyznaczyć praktycznie wystarczające/optymalne liczebności dla CP i ML na podstawie wyników, a nie sztywnej reguły.

- [ ] **P1** Przygotować eksperyment *learning curves* osobno dla etapu CP oraz meta-klasyfikatora ML.
- [ ] **P1** Sprawdzić kilka rosnących liczebności próbek **na klasę** (np. 100, 250, 500, 1 000, 2 000, 5 000 — wyłącznie gdy tyle danych jest dostępnych); każdy wariant powtórzyć dla kilku ziaren losowania.
- [ ] **P1** Testować nie tylko wielkość `ml_train`, ale również `cp_train` i proporcję puli referencyjnej między `cp` i `ml` (`cp_stage_prop`).
- [ ] **P1** Zbadać wpływ podziałów `cp_tt_prop`, `ml_tt_prop`, balansu klas oraz limitów `max_n_cp` i `max_n_ml`.
- [ ] **P1** Porównać Ścieżkę A (`lnc` osobno) i B (całe `nc`), szczególnie przy słabszej adnotacji gatunków niemodelowych.
- [ ] **P1** Utrzymywać niezależną pulę końcową do ewaluacji; wielokrotne losowania i wybór parametrów ograniczyć do części rozwojowej.
- [ ] **P1** Raportować: ROC-AUC, PR-AUC, czułość, specyficzność, precision, F1/MCC (zależnie od celu), kalibrację, przedziały ufności oraz koszt obliczeń.
- [ ] **P1** Zbudować wykresy wyniku względem liczby próbek wraz ze zmiennością między powtórzeniami; zdefiniować **praktyczne kryterium plateau** (minimalna poprawa względem kosztu, bez zakładania z góry konkretnej granicy).
- [ ] **P1** Dodatkowo ocenić wpływ liczebności **niezależnego zbioru testowego** na szerokość przedziałów ufności i stabilność wyników.
- [ ] **P1** Zrewidować progi `<600`, `600–1500`, `>1500` oraz sugestię `1500–2500 na klasę` w obecnym README: na razie traktować je jako **heurystyki do sprawdzenia**, nie udowodnione progi uniwersalne.

**Produkt końcowy:** tabela rekomendowanych liczebności i proporcji z uzasadnieniem, krzywe uczenia, przedziały ufności i opis ograniczeń zależnych od gatunku/danych.

## 3. Katalog cech dla meta-klasyfikatora ML — P1

**Założenie:** podstawą są wyniki predyktorów CP; do nich dodajemy informację biologiczną i strukturalną. Cechy podzielić na **podstawowe** (dostępne szeroko) i **opcjonalne** (zależne od adnotacji i technologii).

### A. Wyniki narzędzi coding-potential (wejście obowiązkowe)

- [ ] Zebranie i normalizacja wyników CPC2, CPAT, CNCI, PLEK i innych docelowych narzędzi: wynik ciągły, klasa, skala i interpretacja.
- [ ] Zbadanie zgodności/rozbieżności predyktorów, brakujących wyników i korelacji między ich cechami.

### B. Cechy sekwencyjne — zwykle z FASTA

- [ ] Długość transkryptu; zawartość GC; częstość k-merów / skład nukleotydowy.
- [ ] Długość najdłuższego ORF, udział ORF w transkrypcie, liczba potencjalnych ORF.
- [ ] Skład kodonów / hexamery / miary podobieństwa do sekwencji kodujących (sprawdzić, czy nie dublują predyktorów CP).
- [ ] Motywy poliadenylacji (*polyadenylation signals*, np. PAS). **Odróżnić przewidywany sygnał od rzeczywistego, zaobserwowanego ogona poly(A)**; sam FASTA nie zawsze pozwala potwierdzić poliadenylację.
- [ ] Zawartość sekwencji powtórzonych i elementów transpozycyjnych (jeśli dostępne adnotacje).

### C. Struktura transkryptu — wymaga GTF/GFF i często genomu

- [ ] Liczba eksonów; długości eksonów i intronów; długość całkowita i rozkład długości intronów.
- [ ] Liczba połączeń ekson–ekson, warianty alternatywnego splicingu, *retained introns*.
- [ ] Kanoniczne (`GT–AG`, `GC–AG`, `AT–AC`) i niekanoniczne miejsca splicingowe; udział nietypowych złączy i jakość ich poparcia odczytami.
- [ ] Liczba nowych / niezaadnotowanych *splice junctions* względem referencji.
- [ ] Ocena, które cechy da się wyliczyć bezpiecznie dla referencji i dla transkryptów StringTie (bez asymetrii jakości danych).

### D. Położenie względem genów i genomu — wymaga adnotacji genomowej

- [ ] Klasy lokalizacji: `intergenic`, `intronic`, `exonic-overlapping`, `sense`, `antisense`, `bidirectional/divergent`, `convergent`.
- [ ] Nakładanie na znane geny kodujące / niekodujące i kierunek transkrypcji; odległość od najbliższego genu/TSS.
- [ ] Liczba nakładających się transkryptów/genów, położenie w obrębie genu i relacja do granic znanych eksonów.

### E. Ekspresja i wiarygodność transkryptu — wymaga danych RNA-seq

- [ ] Poziom ekspresji i wykrywalność w próbkach; spójność między replikatami.
- [ ] Specyficzność tkankowa / warunkowa (gdy dostępne odpowiednie dane).
- [ ] Poparcie pokryciem odczytów, junction reads, wiarygodność końców transkryptu.

### F. Zaawansowane / opcjonalne

- [ ] Konserwacja sekwencji / syntenia między gatunkami.
- [ ] Homologia do białek / domen białkowych; dodatkowe dowody translacji (np. Ribo-seq), jeśli dostępne.
- [ ] Przewidywana struktura wtórna RNA i cechy stabilności — ocenić koszt i użyteczność.
- [ ] Dane o lokalizacji subkomórkowej, jeśli wiarygodne i porównywalne.

### Zasady wyboru cech

- [ ] **P1** Utworzyć tabelę cech: `feature`, definicja, biologiczna interpretacja, typ danych, źródło (FASTA/GTF/genom/BAM/...), narzędzie R, dostępność dla `ex`, potencjalny błąd/bias, koszt obliczeń, priorytet.
- [ ] **P1** Oddzielić panel **Core** (transkrypt FASTA + wyniki CP) od panelu **Extended** (GTF/genom/RNA-seq/dodatkowe eksperymenty).
- [ ] **P1** Unikać cech, które ujawniają etykietę referencyjną zamiast właściwości biologicznej (np. bezpośredni biotyp z adnotacji), i cech dostępnych tylko dla danych treningowych.
- [ ] **P1** Przeprowadzić analizę braków danych, korelacji, ważności cech i *ablation study* (Core vs Core + struktura vs Core + genom vs pełny model).
- [ ] **P1** Uwzględnić przenośność modelu na gatunki słabo zaadnotowane; ocenić warianty bez cech zależnych od rozbudowanej adnotacji.

**Produkt końcowy:** katalog cech z wymaganiami wejściowymi, moduły ich obliczania oraz wyniki porównania zestawów cech.

## 4. Literatura i uzasadnienie statystyczne — P1

- [ ] **P1** Zgromadzić **zweryfikowane publikacje** (pełna bibliografia, DOI/URL, kluczowy wniosek, zakres stosowalności) dotyczące minimalnej liczebności prób dla klasyfikatorów ML.
- [ ] **P1** Sprawdzić, do jakich modeli i sytuacji odnosi się reguła **EPV** (Peduzzi/Harrell i rozwinięcia). **Nie przenosić automatycznie „15–20 obserwacji na zmienną” na RF/XGBoost/LightGBM.**
- [ ] **P1** Zebrać publikacje o *learning curves*, efektywności danych, doborze próby oraz analizie stabilności modeli.
- [ ] **P1** Udokumentować metody wyznaczania przedziałów ufności dla ROC-AUC, PR-AUC, sensitivity/specificity itp. oraz doboru liczebności zbioru testowego (np. bootstrap, właściwe metody dla konkretnych miar).
- [ ] **P1** Zebrać literaturę o zależności próbek, *data leakage*, transkryptach z tego samego genu, podobnych sekwencjach, walidacji międzygatunkowej i zewnętrznej.
- [ ] **P1** Wyszukać publikacje uzasadniające poszczególne grupy cech lncRNA oraz dowody ich przydatności; odróżniać dowód biologiczny od przewagi predykcyjnej w konkretnym modelu.
- [ ] **P1** Przejrzeć istniejący `scripts/readme.md` i oznaczyć twierdzenia zbyt kategoryczne lub pozbawione źródeł (zwłaszcza EPV, „plateau 1500–2500”, „<0,2%”, sztywne progi przejścia A/B).
- [ ] **P1** Do każdej istotnej decyzji metodologicznej dopisać: **co wiemy z literatury / co jest hipotezą / co sprawdziliśmy we własnych testach**.

**Produkt końcowy:** krótka bibliografia komentowana + tabela `decyzja → źródło → ograniczenia → wynik własnej walidacji`.

## 5. Reprodukowalność, testy i wyniki — P1/P2

- [ ] **P1** Dodać automatyczne testy jednostkowe dla `load_and_dedup()`, `inspect_counts()`, `split_datasets()`, `extract_subset()`.
- [ ] **P1** Przygotować małe, sztuczne zbiory FASTA do testów (duplikaty, różne nagłówki, kilka izoform, brak `lnc`, brak `ex`, niezbalansowane klasy).
- [ ] **P1** Raportować liczbę rekordów **przed i po** każdym kroku: wczytanie, deduplikacja, oznaczanie klas, etap CP/ML, train/test, odrzucenia przez limity i `ex_pred`.
- [ ] **P1** Eksportować tabelę pochodzenia rekordów (`ID`, `old_name`, `raw_group`, `class`, `stage`, `tt_set`, `role`, `keep`, przyczyna odrzucenia) do TSV/CSV; zachować oryginalne identyfikatory.
- [ ] **P1** Rejestrować `seed`, źródła i wersje referencji (np. Ensembl/release/assembly), parametry wywołań, `sessionInfo()` i wersje oprogramowania.
- [ ] **P2** Przygotować przykład uruchomienia od wejściowych FASTA do wynikowych FASTA oraz tabel kontrolnych.
- [ ] **P2** Docelowo umożliwić niezależne odtwarzanie testów/analiz dla myszy i wybranych gatunków niemodelowych.

## 6. Dokumentacja i organizacja repozytorium — P2

**Obecny stan katalogu `scripts/` (sprawdzony):** `prepare_dataset.R`, `usage_prepare_dataset.R`, `readme.md`, `IAChat_DatasetPrep.txt`.

- [ ] **P2** Na początku `usage_prepare_dataset.R` dodać `source("prepare_dataset.R")` lub opisać wymagany sposób załadowania funkcji; przykład powinien być samodzielny.
- [ ] **P2** Rozdzielić krótki **Quick Start / README** od szczegółowego opisu metodologii; unikać drugiej pełnej kopii kodu R w README (łatwo o rozjazd wersji).
- [ ] **P2** Dodać tabelę parametrów (`extract_lnc`, `lnc_pattern`, `dedup_mode`, `cp_stage_prop`, `cp_tt_prop`, `ml_tt_prop`, `max_n_cp`, `max_n_ml`, `balance_train`, `ml_neg_mode`) z domyślnymi wartościami i skutkami.
- [ ] **P2** Opisać kontrakt wejścia/wyjścia: formaty FASTA/GTF, wymagane pakiety (R, Biostrings), format `legend`, lista generowanych zbiorów, sposób zapisania FASTA.
- [ ] **P2** Poprawić niespójność nazwy skryptu w obecnym README (`dataset_preparation.R`) względem pliku w repo (`prepare_dataset.R`).
- [ ] **P2** Włączyć `scripts/` do głównego README projektu `CodPod_containers`, jeśli opis repo ma obejmować również moduł przygotowania danych.
- [ ] **P2** Rozważyć przeniesienie `IAChat_DatasetPrep.txt` do katalogu archiwalnego (`docs/archive/`) i oddzielenie historii dyskusji od dokumentacji użytkowej.
- [ ] **P2** **Utrzymać prostą, modularną strukturę**: kod funkcji, przykład, metodologia, testy, ewentualnie przykładowe dane. Nie mnożyć katalogów bez realnej potrzeby.

**Docelowy kierunek (orientacyjny; nie oznacza polecenia reorganizacji teraz):**

```text
scripts/
├── README.md                         # krótka instrukcja i nawigacja
├── prepare_dataset.R                 # funkcje
├── usage_prepare_dataset.R           # przykład użycia
├── docs/
│   ├── dataset_preparation.md        # szczegółowa metodologia
│   └── archive/                      # materiały historyczne
├── tests/                            # testy + małe dane testowe
└── todo.md                           # backlog (opcjonalnie)
```

## 7. Decyzje otwarte / do omówienia

- [ ] Czy **CP** porównuje CDS z ncRNA, czy pełne transkrypty kodujące z ncRNA? Czy do czyszczenia `ex` potrzebna jest osobna referencja?
- [ ] Czy i kiedy stosować *gene-aware split* (test dodatkowy vs główna analiza na poziomie transkryptów)?
- [ ] Jakie kryteria biologiczne/ilościowe powinny decydować o wyborze Ścieżki A/B?
- [ ] Czy klasa negatywna ML ma obejmować tylko `cds`, czy też `nc_other` — i jak wpływa to na interpretację wyniku „lncRNA”?
- [ ] Jak traktować transkrypty częściowe i niepełne oraz różną jakość adnotacji gatunków?
- [ ] Które cechy mają być obowiązkowe w modelu przenośnym na gatunki niemodelowe?
- [ ] Jaki protokół końcowej walidacji przyjąć: niezależny test, split po genach, walidacja międzygatunkowa, dane z innego eksperymentu?

## 8. Dziennik wyników / notatki

| Data | Zadanie / eksperyment | Wynik / obserwacja | Decyzja / kolejny krok |
|---|---|---|---|
| 2026-11-22 | Utworzenie backlogu | Zebrano tematy z przeglądu kodu, dokumentacji i kolejnych uwag | Najpierw P0 (poprawność i leakage), potem eksperymenty P1 |
| — | — | — | — |

---

**Status:** lista propozycji i hipotez do przetestowania; zadania są niezrealizowane, dopóki nie zostaną oznaczone jako ukończone. **Żadne zmiany w repozytorium nie zostały wykonane przy przygotowaniu tego notatnika.**
