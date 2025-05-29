# T5_5_streamlit_druglikeness_app.py
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(page_title="Drug-likeness Checker", layout="wide")
st.title("🧪 Drug-likeness Checker with RDKit")

st.markdown("""
This app evaluates molecular drug-likeness based on **Lipinski's Rule of Five** and other physicochemical descriptors using RDKit.

- Input SMILES string(s) (one per line)
- View computed descriptors and Lipinski rule check
""")

# Input box
smiles_input = st.text_area("Enter one SMILES per line:", height=200)

# Parse input
def analyze_smiles(smiles_list):
    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            molwt = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hdonor = Lipinski.NumHDonors(mol)
            hacceptor = Lipinski.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            rotb = Descriptors.NumRotatableBonds(mol)
            lipinski_pass = (
                molwt <= 500 and logp <= 5 and hdonor <= 5 and hacceptor <= 10
            )
            results.append({
                "SMILES": smi,
                "MolWt": molwt,
                "LogP": logp,
                "HDonors": hdonor,
                "HAcceptors": hacceptor,
                "TPSA": tpsa,
                "RotBonds": rotb,
                "Lipinski_Pass": lipinski_pass
            })
    return pd.DataFrame(results)

if st.button("Analyze"):
    smiles_list = [s.strip() for s in smiles_input.splitlines() if s.strip()]
    if smiles_list:
        df_result = analyze_smiles(smiles_list)
        st.success(f"✅ Processed {len(df_result)} molecule(s)")
        st.dataframe(df_result)

        # Filter for Lipinski passed molecules
        passed = df_result[df_result['Lipinski_Pass']]
        failed = df_result[~df_result['Lipinski_Pass']]

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Lipinski Passed", len(passed))
        with col2:
            st.metric("Lipinski Violated", len(failed))

        # Visualization
        st.subheader("Descriptor Distribution")
        if len(df_result) > 1:
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.histplot(df_result['MolWt'], bins=20, kde=True, ax=ax)
            ax.set_title("Molecular Weight Distribution")
            st.pyplot(fig)

            fig2, ax2 = plt.subplots(figsize=(10, 5))
            sns.histplot(df_result['LogP'], bins=20, kde=True, ax=ax2)
            ax2.set_title("LogP Distribution")
            st.pyplot(fig2)
    else:
        st.warning("⚠️ Please enter at least one SMILES.")
