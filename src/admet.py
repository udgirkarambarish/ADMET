import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen

# --- Helper Functions for ADMET Scoring ---

def score_optimal_range(value, r_low, r_opt_low, r_opt_high, r_high):
    """Calculates a desirability score (0-1) for a value within an optimal range."""
    if pd.isna(value): return 0
    if r_opt_low <= value <= r_opt_high:
        return 1.0
    elif value < r_opt_low and value > r_low:
        return (value - r_low) / (r_opt_low - r_low)
    elif value > r_opt_high and value < r_high:
        return (r_high - value) / (r_high - r_opt_high)
    else:
        return 0.0

def score_increasing(value, r_low, r_high):
    """Calculates a desirability score (0-1) where higher is better."""
    if pd.isna(value): return 0
    if value <= r_low: return 0.0
    if value >= r_high: return 1.0
    return (value - r_low) / (r_high - r_low)

def score_decreasing(value, r_low, r_high):
    """Calculates a desirability score (0-1) where lower is better."""
    if pd.isna(value): return 0
    if value <= r_low: return 1.0
    if value >= r_high: return 0.0
    return (r_high - value) / (r_high - r_low)

# --- ADMET Sub-score Calculation Functions ---

def calculate_absorption_score(row):
    """Calculates the Absorption sub-score."""
    scores = {
        'logS': score_optimal_range(row['logS'], -6, -4, -1, 0),
        'caco2': score_increasing(row['caco2'], 0, 20),
        'hia': 1.0 if row['hia'] == 1 else 0.0,
        'pgp_sub': 1.0 if row['pgp_sub'] == 0 else 0.0,
    }
    return np.mean(list(scores.values()))

def calculate_distribution_score(row):
    """Calculates the Distribution sub-score."""
    scores = {
        'PPB': score_optimal_range(row['PPB'], 0, 50, 95, 100),
        'BBB': score_decreasing(row['BBB'], -1, 1), # Lower is generally safer
        'logVDss': score_optimal_range(row['logVDss'], -2, -1, 1, 2)
    }
    return np.mean(list(scores.values()))

def calculate_metabolism_score(row):
    """Calculates the Metabolism sub-score."""
    cyp_cols = ['CYP1A2-inh', 'CYP2C9-inh', 'CYP2C19-inh', 'CYP2D6-inh', 'CYP3A4-inh']
    scores = [1.0 if row[col] == 0 else 0.0 for col in cyp_cols]
    return np.mean(scores)

def calculate_excretion_score(row):
    """Calculates the Excretion sub-score."""
    scores = {
        'cl-plasma': score_optimal_range(row['cl-plasma'], 0, 5, 30, 60),
        'BCRP': 1.0 if row['BCRP'] == 0 else 0.0 # Non-interaction is desirable
    }
    return np.mean(list(scores.values()))

def calculate_toxicity_score(row):
    """Calculates the Toxicity sub-score."""
    tox_cols = ['hERG', 'DILI', 'Ames', 'SkinSen', 'Carcinogenicity']
    scores = [1.0 if row[col] == 0 else 0.0 for col in tox_cols]
    return np.mean(scores)

# --- RDKit and Drug-Likeness Functions ---

def get_rdkit_desc(smiles_string, descriptor_func):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return descriptor_func(mol) if mol else None
    except:
        return None

def apply_lipinski_rules(row):
    violations = 0
    if row['MW'] > 500: violations += 1
    if row['nHD'] > 5: violations += 1
    if row['nHA'] > 10: violations += 1
    if not (-2 <= row['logP'] <= 5): violations += 1
    if row['nRot'] > 5: violations += 1
    return violations

def apply_veber_rules(row):
    violations = 0
    if row['MW'] >= 500: violations += 1
    if row['nHD'] > 5: violations += 1
    if row['nHA'] > 10: violations += 1
    if row['TPSA'] >= 140: violations += 1
    if row['Atom_Count'] is None or not (20 <= row['Atom_Count'] <= 70): violations += 1
    return violations

def apply_ghose_filter(row):
    violations = 0
    if not (160 <= row['MW'] <= 480): violations += 1
    if not (-0.4 <= row['logP'] <= 5.6): violations += 1
    if row['Atom_Count'] is None or not (20 <= row['Atom_Count'] <= 70): violations += 1
    if row['Molar_Refractivity'] is None or not (40 <= row['Molar_Refractivity'] <= 130): violations += 1
    return violations

# --- Main Analysis Function ---

def analyze_drug_properties(input_file, output_file):
    """
    Reads a molecule dataset, applies comprehensive drug-likeness filters and
    a 5-part weighted ADMET score, and saves a full summary to an Excel file.
    """
    try:
        df = pd.read_excel(input_file, engine='openpyxl')

        # --- Step 1: Check for all required columns ---
        required_cols = [
            'Protein sequence', 'smiles', 'MW', 'logP', 'nHD', 'nHA', 'nRot', 'TPSA',
            'logS', 'caco2', 'hia', 'pgp_sub', 'PPB', 'BBB', 'logVDss',
            'CYP1A2-inh', 'CYP2C9-inh', 'CYP2C19-inh', 'CYP2D6-inh', 'CYP3A4-inh',
            'cl-plasma', 'BCRP', 'hERG', 'DILI', 'Ames', 'SkinSen', 'Carcinogenicity'
        ]
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            return f"Error: Missing required columns: {', '.join(missing)}"

        # --- Step 2: Calculate RDKit Descriptors ---
        print("Calculating RDKit descriptors...")
        df['Atom_Count'] = df['smiles'].apply(lambda smi: get_rdkit_desc(smi, lambda m: m.GetNumAtoms()))
        df['Molar_Refractivity'] = df['smiles'].apply(lambda smi: get_rdkit_desc(smi, Crippen.MolMR))
        print("RDKit calculations complete.")

        # --- Step 3: Run Drug-Likeness Analysis ---
        df['Lipinski_Violations'] = df.apply(apply_lipinski_rules, axis=1)
        df['Passes_Lipinski'] = df['Lipinski_Violations'] < 2
        df['Veber_Violations'] = df.apply(apply_veber_rules, axis=1)
        df['Passes_Veber'] = df['Veber_Violations'] == 0
        df['Ghose_Violations'] = df.apply(apply_ghose_filter, axis=1)
        df['Passes_Ghose'] = df['Ghose_Violations'] == 0

        # --- Step 4: Calculate 5-Part ADMET Score ---
        print("Calculating 5-part ADMET scores...")
        df['Score_Absorption'] = df.apply(calculate_absorption_score, axis=1)
        df['Score_Distribution'] = df.apply(calculate_distribution_score, axis=1)
        df['Score_Metabolism'] = df.apply(calculate_metabolism_score, axis=1)
        df['Score_Excretion'] = df.apply(calculate_excretion_score, axis=1)
        df['Score_Toxicity'] = df.apply(calculate_toxicity_score, axis=1)

        df['ADMET_Score'] = 0.2 * (df['Score_Absorption'] + df['Score_Distribution'] +
                                   df['Score_Metabolism'] + df['Score_Excretion'] +
                                   df['Score_Toxicity'])
        print("ADMET scoring complete.")

        # --- Step 5: Create and Save Final Summary Output ---
        summary_cols = [
            'Protein sequence', 'smiles', 'ADMET_Score', 
            'Score_Absorption', 'Score_Distribution', 'Score_Metabolism', 
            'Score_Excretion', 'Score_Toxicity',
            'Passes_Lipinski', 'Passes_Veber', 'Passes_Ghose',
            'Lipinski_Violations', 'Veber_Violations', 'Ghose_Violations',
            'Atom_Count', 'Molar_Refractivity'
        ]
        summary_df = df[summary_cols]
        summary_df.to_excel(output_file, index=False, engine='openpyxl')
        
        return f"Analysis complete. Full summary saved to '{output_file}'"

    except FileNotFoundError:
        return f"Error: The file was not found at {input_file}"
    except Exception as e:
        return f"An unexpected error occurred: {e}"

# --- Main execution block ---
if __name__ == "__main__":
    input_excel_file = 'data/Validation_dataset.xlsx'
    output_excel_file = 'output.xlsx'

    result_message = analyze_drug_properties(input_excel_file, output_excel_file)
    print(result_message)
