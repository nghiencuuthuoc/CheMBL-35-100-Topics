import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings

# Tắt cảnh báo UserWarning
warnings.filterwarnings("ignore", category=UserWarning)

# Cấu hình thư mục
base_path = ".."
raw_path = os.path.join(base_path, "data", "T4_lipinski_raw_with_descriptors.csv")
filtered_path = os.path.join(base_path, "data", "T4_lipinski_filtered.csv")

# Load dữ liệu
@st.cache_data
def load_data():
    raw_df = pd.read_csv(raw_path)
    filtered_df = pd.read_csv(filtered_path)
    return raw_df, filtered_df

raw_df, filtered_df = load_data()

# Giao diện
st.title("🔬 Lipinski Rule Dashboard – ChEMBL 35")
st.markdown("Visualize and compare descriptor distributions before and after Lipinski filtering.")

# Tùy chọn đặc trưng
lipinski_columns = ['MolWt', 'LogP', 'HDonors', 'HAcceptors']
selected_col = st.selectbox("📌 Select descriptor:", lipinski_columns)

# Overlay histogram
fig, ax = plt.subplots()
sns.histplot(raw_df[selected_col], bins=30, stat="density", kde=True, label="Before", color='gray', alpha=0.5)
sns.histplot(filtered_df[selected_col], bins=30, stat="density", kde=True, label="After", color='blue', alpha=0.6)

plt.title(f"Distribution of {selected_col}")
plt.xlabel(selected_col)
plt.ylabel("Density")
plt.legend()
st.pyplot(fig)

# Xem dữ liệu
with st.expander("📋 Show raw data"):
    st.dataframe(raw_df.head())

with st.expander("✅ Show filtered data"):
    st.dataframe(filtered_df.head())
