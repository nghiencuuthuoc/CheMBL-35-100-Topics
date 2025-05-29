# app/T1_3_qsar_streamlit.py
import streamlit as st
import pandas as pd
import joblib
import os
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from qsar_utils import calc_rdkit_descriptors


st.set_page_config(page_title="QSAR IC50 Explorer", layout="wide")

st.title("🧪 QSAR Model on ChEMBL IC50 Data")

# Load model
model_path = os.path.join("..", "model", "qsar_model.pkl")
model = joblib.load(model_path)

# Load data
# data_path = os.path.join("..", "data", "T1_ic50_data.csv")
data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "T1_ic50_data.csv"))

df = pd.read_csv(data_path)
st.subheader("📄 Raw ChEMBL Data (IC50 in nM)")
st.dataframe(df.head(10))

# Choose molecule
st.subheader("🔍 Predict pIC50 for new SMILES")
input_smiles = st.text_input("Enter SMILES", value="CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
desc = calc_rdkit_descriptors(input_smiles)

if desc:
    desc_df = pd.DataFrame([desc])
    prediction = model.predict(desc_df)[0]
    st.success(f"✅ Predicted pIC50: **{prediction:.3f}**")
    st.write("Descriptors used:", desc)
else:
    st.warning("Invalid SMILES or descriptor calculation failed.")

# Optional: visualize pIC50 distribution
st.subheader("📊 pIC50 Distribution (calculated)")
df = df.dropna(subset=["canonical_smiles", "ic50_nM"])
df["pIC50"] = -df["ic50_nM"].apply(lambda x: float(x) * 1e-9).apply(lambda x: pd.np.log10(x)) * -1

st.line_chart(df["pIC50"])
