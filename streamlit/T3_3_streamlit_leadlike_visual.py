# File: T3_3_streamlit_leadlike_visual.py
# Author: Nghiên Cứu Thuốc | www.nghiencuuthuoc.com

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 📁 Load data
base_path = '..'
path_all = os.path.join(base_path, 'data', 'T3_leadlike_with_descriptors.csv')
path_filtered = os.path.join(base_path, 'data', 'T3_leadlike_filtered_descriptors.csv')

st.set_page_config(page_title="Lead-like Compound Visualization", layout="wide")

# 🧬 Title
st.title("🧪 Lead-like Compound Descriptor Visualization")
st.markdown("Visual comparison of molecular properties of compounds from **ChEMBL 35**.")

# 📊 Load data
@st.cache_data
def load_data():
    df_all = pd.read_csv(path_all)
    df_filtered = pd.read_csv(path_filtered)
    return df_all, df_filtered

df_all, df_filtered = load_data()

# 🔘 Select descriptor
descriptor = st.selectbox("Select descriptor to visualize:", ['MolWt', 'LogP', 'HBD', 'HBA', 'RotBonds'])

# 📈 Plot histogram
fig, ax = plt.subplots(figsize=(8, 4))
sns.histplot(df_all[descriptor], color='blue', label='All Compounds', kde=True, stat='density', bins=30, ax=ax)
sns.histplot(df_filtered[descriptor], color='green', label='Lead-like', kde=True, stat='density', bins=30, ax=ax)
ax.set_title(f'Distribution of {descriptor}')
ax.set_xlabel(descriptor)
ax.set_ylabel("Density")
ax.legend()
st.pyplot(fig)

# 📑 Show summary stats
with st.expander("📊 Show summary statistics"):
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("All Compounds")
        st.dataframe(df_all[descriptor].describe())
    with col2:
        st.subheader("Lead-like Only")
        st.dataframe(df_filtered[descriptor].describe())

# 🔁 Optional: pairplot
if st.checkbox("📌 Show pairplot (lead-like only)", value=False):
    with st.spinner("Generating pairplot..."):
        fig2 = sns.pairplot(df_filtered[['MolWt', 'LogP', 'HBD', 'HBA', 'RotBonds']])
        st.pyplot(fig2)

# 📬 Footer
st.markdown("---")
st.markdown("""
| Copyright 2025 | 🧠 Nghiên Cứu Thuốc | PharmApp |
| www.nghiencuuthuoc.com | Zalo: +84888999311 |
""")
