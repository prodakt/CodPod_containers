# CodPod_containers

*A unified collection of Docker containers, benchmark datasets, and reproducible workflows for transcript coding potential prediction.*

---

## Overview

**CodPod_containers** is an open-source companion project developed to support the **lncRna** framework. Its primary objective is to provide a standardized and reproducible environment for running coding potential prediction software without the need for complicated installation or dependency management.

The repository gathers Docker containers, benchmark datasets, example analyses, and standardized execution procedures for the most commonly used coding potential prediction tools.

Instead of spending hours installing individual programs with incompatible dependencies, users should be able to execute any supported predictor using a single Docker command.

Although the repository has been created as part of the **lncRna** ecosystem, it can be used independently in any transcriptomics or lncRNA analysis project.

---

# Motivation

Coding potential prediction is one of the fundamental steps during long non-coding RNA (lncRNA) identification.

Over the last decade, numerous prediction tools have been developed. They employ diverse computational strategies, including analyses based on

- Open Reading Frame (ORF) characteristics
- k-mer composition
- codon usage bias
- sequence intrinsic features
- GC content
- protein homology
- evolutionary conservation
- classical machine learning
- deep learning
- ensemble learning

Each program has its own installation procedure, software requirements, dependencies, input formats, output formats, and execution syntax. Consequently, comparing multiple predictors or combining them into a single workflow is often difficult and time-consuming.

**CodPod_containers** aims to eliminate these obstacles by providing a unified infrastructure in which all supported tools can be executed in a standardized manner.

---

# Project goals

The project has several objectives:

- collect Docker images for major coding potential prediction tools;
- prepare custom Docker containers whenever official images are unavailable;
- provide standardized benchmark datasets;
- prepare reproducible execution procedures;
- document input and output formats;
- collect example prediction results;
- simplify comparative benchmarking;
- facilitate integration with the **lncRna** package.

Ultimately, users should be able to choose any subset of prediction tools and execute them using nearly identical workflows regardless of their original implementation.

---

# Repository organization

```text
CodPod_containers/

├── containers/
│   ├── CPC2/
│   ├── CPAT/
│   ├── CNCI/
│   ├── FEELnc/
│   ├── LGC/
│   ├── LncADeep/
│   ├── LncFinder/
│   ├── PLEK/
│   └── ...
│
├── datasets/
│   ├── training/
│   ├── testing/
│   ├── benchmark/
│   └── examples/
│
├── sequences/
│   ├── fasta/
│   ├── annotation/
│   └── reference/
│
├── results/
│   ├── CPC2/
│   ├── CPAT/
│   ├── CNCI/
│   └── ...
│
├── workflows/
│
├── documentation/
│
├── metadata/
│
└── README.md
```

---

# Supported predictors

The repository aims to include Docker containers for the most widely used transcript coding potential prediction software.

Current target software includes:

| Predictor | Primary methodology |
|------------|---------------------|
| CPC2 | Support Vector Machine |
| CPAT | Logistic Regression |
| CNCI | Intrinsic sequence features |
| PLEK | Improved k-mer frequencies |
| LGC | ORF length and GC content |
| CPPred | Multiple sequence-derived features |
| LongDist | Distance-based classifier |
| FEELnc | Random Forest |
| LncFinder | Machine learning |
| LncADeep | Deep learning |
| PredLnc-GFStack | Ensemble learning |

Additional predictors are welcome.

---

# What is provided for each predictor?

Whenever possible, each software package will include:

- Docker image or Dockerfile
- links to official software repositories
- software version information
- installation notes
- execution examples
- input data specification
- output description
- benchmark dataset
- expected output files
- algorithm summary
- publication reference

---

# Standardized workflow

Each predictor is prepared according to the same general procedure.

1. Obtain the Docker image.
2. Run the prediction using benchmark sequences.
3. Export prediction results.
4. Validate output.
5. Compare results with other predictors.

This standardized workflow enables direct comparison between different coding potential prediction methods.

---

# Benchmark datasets

The repository contains standardized datasets intended for

- model training
- testing
- benchmarking
- reproducibility studies
- software validation

Whenever licensing permits, publicly available reference datasets will be included. Additional benchmark datasets generated during the development of the **lncRna** package may also be provided.

---

# Docker philosophy

Whenever an official Docker image exists, it will be preferred.

Priority order:

1. Official Docker images maintained by software authors.
2. Community-maintained Docker images (e.g. BioContainers).
3. Custom Docker images developed within this repository.

This strategy maximizes software reproducibility while minimizing maintenance effort.

---

# Contribution guidelines

Contributors are encouraged to add new predictors.

Each contribution should include:

- Docker image or Dockerfile
- software version
- execution instructions
- benchmark example
- expected output
- software citation
- short description of the prediction algorithm

---

# Long-term vision

The long-term objective of **CodPod_containers** is to become a comprehensive collection of reproducible coding potential prediction environments.

By integrating standardized Docker containers with benchmark datasets and unified execution procedures, the repository aims to make comparative evaluation of transcript coding potential predictors straightforward, reproducible, and platform-independent.

Ultimately, users of the **lncRna** package—or any other transcriptomics workflow—should be able to select one or more prediction tools and execute them with minimal effort:

```bash
docker run ...
```

No manual installation.

No dependency conflicts.

No operating system issues.

Just reproducible analyses.
