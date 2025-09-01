Markdown

# ADMET and Drug-Likeness Analysis Project

This document contains a complete guide to the project, including a brief README and installation dependencies.

## 1. README

This project analyzes drug candidates from an Excel file. It runs them through drug-likeness filters (Lipinski, Veber, Ghose) and calculates a comprehensive 5-part ADMET score (Absorption, Distribution, Metabolism, Excretion, Toxicity) to rank their potential. The final output is a summary Excel file.

**How to Run:**

1. Set up the environment and install dependencies (see section 2).
2. Place your `Validation_dataset.xlsx` file in a `data/` folder.
3. Place the Python script (from section 3) in a `src/` folder.
4. From the `src` folder, run `python src\analyze_drugs.py`. The `output.xlsx` will be generated in the parent directory.

## 2. Requirements (`requirements.txt`)

Save the following content as `requirements.txt` and install using `pip install -r requirements.txt`.

**Note:** RDKit is best installed first via conda/mamba: `conda install -c conda-forge rdkit`.
