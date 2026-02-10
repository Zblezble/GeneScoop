## 🧬 GeneScoop

Detect and extract multiple gene sequences from multiple GenBank files.

---

## 📂 Contents

- GeneScoop.py – The code itself
- GeneScoop.ipynb – The code formatted for Jupyter Lab (how I used it)
- LICENSE – Licensing information  
- README.md – This file
- ncbi_queries.md – Example of NCBI queries used to download genome files (using [ncbi-acc-download](https://github.com/kblin/ncbi-acc-download), or to download sequences using [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt)) which can then be subjected to GeneScoop, same keywords can be then used in the code at `product_and_note_keywords`

---

## 🔍 Features

- Extracts **CDS** features based on:
  - Exact gene name match (e.g., `mcrA`)
  - Flexible product or note keyword matching
- Outputs sequences in **FASTA** format
- Offers parallel processing

---

## 🛠️ Usage

The code is formatted for use with Jupyter Lab

GeneScoop is released under a BSD-3-Clause license. See LICENSE for more details.
