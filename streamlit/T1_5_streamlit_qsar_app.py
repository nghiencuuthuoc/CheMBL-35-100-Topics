import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import joblib
import os

# Load model
model_path = os.path.join("..", "model", "qsar_model.pkl")
model = joblib.load(model_path)

# App title
st.set_page_config(page_title="QSAR IC50 Predictor")
st.title("🔬 QSAR Web App – Predict IC50 from SMILES")
st.markdown("Enter a molecule **SMILES** below to predict its IC50 (µM) value.")

# Input
smiles = st.text_input("👉 SMILES input:", "CCOc1ccc2nc(S(N)(=O)=O)sc2c1")

# Function to compute descriptors
def calc_desc(smiles_str):
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None
    return {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'RingCount': Descriptors.RingCount(mol)
    }

# Predict button
if st.button("🚀 Predict IC50"):
    desc = calc_desc(smiles)
    if desc is None:
        st.error("❌ Invalid SMILES. Please try again.")
    else:
        df = pd.DataFrame([desc])
        prediction = model.predict(df)[0]
        st.success(f"✅ Predicted IC50: **{prediction:.2f} µM**")
        st.dataframe(df.T.rename(columns={0: "Value"}))
