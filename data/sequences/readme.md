Data description
This repository contains full coding (cds) and non‑coding (nc) sequence datasets for three species: Homo sapiens, Mus musculus and Sus scrofa.
For each species, coding and non‑coding sequences were split into training and test subsets in a 60/40 ratio (60% training, 40% test).

Sequences are organized under:

data/sequences/Homo_sapiens/

data/sequences/Mus_musculus/

data/sequences/Sus_scrofa/

Each species directory contains four FASTA files:

nc_<Species>_train_balanced.fasta
Non‑coding training set. The dataset is balanced with respect to the smaller class, i.e. the number of nc sequences is matched to the number of pc sequences in the corresponding training set.

pc_<Species>_train_balanced.fasta
Coding training set. The dataset is balanced relative to the less frequent class, so that the training data for coding and non‑coding sequences have the same size. This is intended for training classifiers under class‑balanced conditions.

seqs2predict_<Species>_full.fasta
Full test set for prediction. It contains all sequences from the test subset (both coding and non‑coding), with the original class proportions preserved. No balancing is applied here, to better reflect “natural” prediction conditions where class frequencies may be skewed.

seqs2predict_<Species>_small.fasta
Small validation set for quick evaluation. It consists of 500 coding (pc) and 500 non‑coding (nc) sequences sampled from the test data. This file is intended for rapid sanity‑checks, method prototyping and lightweight benchmarking without loading the full test set.

Train/test split and balancing
Initial cds and nc datasets were split into training and test sets using a 60/40 ratio (60% of sequences assigned to training, 40% to test).

For each species, balanced training sets were constructed by downsampling the majority class to the size of the minority class. As a result, pc_*_train_balanced.fasta and nc_*_train_balanced.fasta contain the same number of sequences and can be used to train models under balanced class conditions.

Test data are kept in their original, imbalanced form in seqs2predict_*_full.fasta. This allows evaluation of classifiers in a more realistic scenario, where the underlying coding/non‑coding ratio is not artificially adjusted.

The seqs2predict_*_small.fasta files provide a compact, balanced validation subset (500 pc + 500 nc) that can be used for rapid method comparison or debugging.

Using the datasets with Git LFS
Full FASTA files are stored via Git Large File Storage (LFS). To clone the repository and download all data:

bash
git clone https://github.com/prodakt/CodPod_containers.git
cd CodPod_containers
git lfs install
git lfs pull
After these steps, all FASTA files in data/sequences/ will be available locally in their full size and can be directly used as input for classification tools and pipelines.
