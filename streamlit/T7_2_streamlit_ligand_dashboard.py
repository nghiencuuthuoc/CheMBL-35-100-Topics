# File: notebook/T7_2_streamlit_ligand_dashboard.py
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import warnings

# Disable RDKit warnings
warnings.filterwarnings("ignore", category=UserWarning)

# === Load Data ===
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
data_path = os.path.join(base_path, 'data', 'T7_egfr_ligands.csv')

st.set_page_config(page_title="Ligand-Based Drug Design", layout="wide")
st.title("🧪 Ligand-Based Drug Design (EGFR Targets)")
st.caption("Powered by ChEMBL + RDKit")

# Load CSV
df = pd.read_csv(data_path)

# === Compute Descriptors ===
@st.cache_data
def compute_descriptors(smiles_list):
    result = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            result.append({
                "MolWt": Descriptors.MolWt(mol),
                "LogP": Descriptors.MolLogP(mol),
                "TPSA": Descriptors.TPSA(mol)
            })
        else:
            result.append({"MolWt": None, "LogP": None, "TPSA": None})
    return pd.DataFrame(result)

desc_df = compute_descriptors(df['canonical_smiles'])
df_final = pd.concat([df, desc_df], axis=1)

# === Show Data ===
st.subheader("📄 Ligand Table (Top 10)")
st.dataframe(df_final.head(10))

# === Plotting ===
st.subheader("📊 Descriptor Distributions")

cols = st.multiselect("Choose descriptors to visualize:", ["MolWt", "LogP", "TPSA"], default=["MolWt", "LogP"])

for col in cols:
    fig, ax = plt.subplots()
    sns.histplot(df_final[col].dropna(), kde=True, ax=ax, bins=30)
    ax.set_title(f"{col} Distribution")
    st.pyplot(fig)

# === Download ===
st.subheader("⬇ Export Data")
csv_export = df_final.to_csv(index=False)
st.download_button("Download as CSV", csv_export, file_name="T7_egfr_ligands_with_descriptors.csv")

