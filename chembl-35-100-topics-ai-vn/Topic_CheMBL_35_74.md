# PharmApp Suite
## 🧠 AI for Drug Discovery and Development 🧪
| Copyright 2025 | Nghiên Cứu Thuốc | www.nghiencuuthuoc.com | Zalo: +84888999311 |

🧩 Topic: Vẽ 3D không gian hóa học phân tử
---
Tuyệt vời! Tôi rất vui được hỗ trợ bạn trong dự án phân tích dữ liệu ChEMBL 35 sử dụng RDKit. Với vai trò là một chuyên gia nghiên cứu và phát triển dược học, tôi hiểu rõ tầm quan trọng của việc khai thác dữ liệu ChEMBL để tìm kiếm các hợp chất tiềm năng. Dưới đây là phân tích chi tiết, hướng dẫn song ngữ, code mẫu SQL và Python, cùng với các ví dụ để bạn bắt đầu.

**1. Phân tích Mô hình Phân tích (Analysis of the Analysis Model)**

Dự án của bạn tập trung vào việc phân tích dữ liệu ChEMBL 35, một cơ sở dữ liệu lớn chứa thông tin về các phân tử hoạt tính sinh học. Mục tiêu là sử dụng RDKit để trích xuất các đặc trưng hóa học từ cấu trúc phân tử và kết hợp chúng với dữ liệu hoạt tính sinh học để xây dựng mô hình dự đoán.

*   **Data Source (Nguồn dữ liệu):** ChEMBL 35 database.
*   **Tools (Công cụ):**
    *   psql (PostgreSQL) for database querying.
    *   RDKit for chemical feature extraction.
    *   Jupyter Notebook for code execution and documentation.
    *   Scikit-learn (or other machine learning libraries) for model building.
*   **Steps (Các bước):**
    1.  **Data Extraction (Trích xuất dữ liệu):** Use SQL queries to extract relevant data from ChEMBL, such as compound structures (SMILES) and activity data (IC50, Ki, etc.).
    2.  **Feature Generation (Tạo đặc trưng):** Use RDKit to generate molecular descriptors (e.g., fingerprints, physicochemical properties) from the SMILES strings.
    3.  **Data Preprocessing (Tiền xử lý dữ liệu):** Clean and prepare the data for modeling. This may involve handling missing values, scaling features, and splitting the data into training and testing sets.
    4.  **Model Building (Xây dựng mô hình):** Train a machine learning model to predict activity based on the molecular descriptors. Common models include Random Forest, Support Vector Machines, and Neural Networks.
    5.  **Model Evaluation (Đánh giá mô hình):** Evaluate the performance of the model using appropriate metrics (e.g., AUC, RMSE, R-squared).
    6.  **Analysis and Interpretation (Phân tích và giải thích):** Analyze the model to identify important features and gain insights into the structure-activity relationship.

**2. Hướng dẫn song ngữ (Bilingual Instructions)**

Dưới đây là hướng dẫn chi tiết bằng cả tiếng Anh và tiếng Việt:

*   **English:**

    1.  **Connect to the ChEMBL database using psql.**
    2.  **Execute SQL queries to extract the desired data (e.g., compound structures and activity data).**
    3.  **Save the extracted data as CSV files in the `../data/` directory.**
    4.  **Load the CSV files into a Jupyter Notebook.**
    5.  **Use RDKit to generate molecular descriptors from the SMILES strings.**
    6.  **Preprocess the data (handle missing values, scale features, etc.).**
    7.  **Split the data into training and testing sets.**
    8.  **Train a machine learning model to predict activity.**
    9.  **Evaluate the model's performance.**
    10. **Analyze the model to identify important features and gain insights into the structure-activity relationship.**
*   **Tiếng Việt:**

    1.  **Kết nối đến cơ sở dữ liệu ChEMBL bằng psql.**
    2.  **Thực hiện các truy vấn SQL để trích xuất dữ liệu mong muốn (ví dụ: cấu trúc hợp chất và dữ liệu hoạt tính).**
    3.  **Lưu dữ liệu đã trích xuất dưới dạng tệp CSV trong thư mục `../data/`.**
    4.  **Tải các tệp CSV vào Jupyter Notebook.**
    5.  **Sử dụng RDKit để tạo các đặc trưng phân tử từ chuỗi SMILES.**
    6.  **Tiền xử lý dữ liệu (xử lý giá trị thiếu, chia tỷ lệ đặc trưng, v.v.).**
    7.  **Chia dữ liệu thành các tập huấn luyện và kiểm tra.**
    8.  **Huấn luyện mô hình học máy để dự đoán hoạt tính.**
    9.  **Đánh giá hiệu suất của mô hình.**
    10. **Phân tích mô hình để xác định các đặc trưng quan trọng và hiểu rõ hơn về mối quan hệ cấu trúc-hoạt tính.**

**3. Code SQL và Python mẫu (Sample SQL and Python Code)**

**SQL (English):**

```sql
-- Extract 100 compounds with IC50 values for a specific target
SELECT
    molregno,
    compound_structures.canonical_smiles,
    activities.standard_value,
    activities.standard_units
FROM
    compound_structures
JOIN
    activities ON compound_structures.molregno = activities.molregno
JOIN
    target_dictionary ON activities.tid = target_dictionary.tid
WHERE
    target_dictionary.pref_name = 'desired_target_name'  -- Replace with your target
    AND activities.standard_type = 'IC50'
    AND activities.standard_units = 'nM'
    AND act.standard_value ~ '^[0-9\.]+$'  -- Only numeric values
LIMIT 100;
```

**SQL (Tiếng Việt):**

```sql
-- Trích xuất 100 hợp chất với giá trị IC50 cho một mục tiêu cụ thể
SELECT
    molregno,
    compound_structures.canonical_smiles,
    activities.standard_value,
    activities.standard_units
FROM
    compound_structures
JOIN
    activities ON compound_structures.molregno = activities.molregno
JOIN
    target_dictionary ON activities.tid = target_dictionary.tid
WHERE
    target_dictionary.pref_name = 'desired_target_name'  -- Thay thế bằng mục tiêu của bạn
    AND activities.standard_type = 'IC50'
    AND activities.standard_units = 'nM'
    AND act.standard_value ~ '^[0-9\.]+$'  -- Chỉ giá trị số
LIMIT 100;
```

**Python (English):**

```python
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Define the base path
base_path = "../data/"

# Load the data from CSV
data = pd.read_csv(os.path.join(base_path, "your_data.csv"))  # Replace with your CSV file

# Function to calculate molecular descriptors
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        descriptors = {}
        descriptors['MolLogP'] = Descriptors.MolLogP(mol)
        descriptors['MolWt'] = Descriptors.MolWt(mol)
        descriptors['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
        descriptors['NumHDonors'] = Descriptors.NumHDonors(mol)
        return descriptors
    else:
        return None

# Apply the function to each SMILES string
data['descriptors'] = data['canonical_smiles'].apply(calculate_descriptors)

# Handle missing descriptors
data = data.dropna(subset=['descriptors'])

# Convert descriptors to DataFrame
descriptors_df = pd.DataFrame(data['descriptors'].tolist())
data = pd.concat([data, descriptors_df], axis=1)

# Prepare data for modeling
X = data[['MolLogP', 'MolWt', 'NumHAcceptors', 'NumHDonors']]  # Use only available descriptors
y = data['standard_value']

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a Random Forest Regressor model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse}")
print(f"R-squared: {r2}")
```

**Python (Tiếng Việt):**

```python
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Xác định đường dẫn cơ sở
base_path = "../data/"

# Tải dữ liệu từ tệp CSV
data = pd.read_csv(os.path.join(base_path, "your_data.csv"))  # Thay thế bằng tệp CSV của bạn

# Hàm tính toán các descriptor phân tử
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        descriptors = {}
        descriptors['MolLogP'] = Descriptors.MolLogP(mol)
        descriptors['MolWt'] = Descriptors.MolWt(mol)
        descriptors['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
        descriptors['NumHDonors'] = Descriptors.NumHDonors(mol)
        return descriptors
    else:
        return None

# Áp dụng hàm cho mỗi chuỗi SMILES
data['descriptors'] = data['canonical_smiles'].apply(calculate_descriptors)

# Xử lý các descriptor bị thiếu
data = data.dropna(subset=['descriptors'])

# Chuyển đổi descriptors thành DataFrame
descriptors_df = pd.DataFrame(data['descriptors'].tolist())
data = pd.concat([data, descriptors_df], axis=1)

# Chuẩn bị dữ liệu cho mô hình hóa
X = data[['MolLogP', 'MolWt', 'NumHAcceptors', 'NumHDonors']]  # Chỉ sử dụng các descriptor có sẵn
y = data['standard_value']

# Chia dữ liệu thành các tập huấn luyện và kiểm tra
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Huấn luyện mô hình Random Forest Regressor
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Dự đoán
y_pred = model.predict(X_test)

# Đánh giá mô hình
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse}")
print(f"R-squared: {r2}")
```

**4. Ví dụ Code SQL và Python (SQL and Python Code Examples)**

Dưới đây là 5 ví dụ về code SQL và Python để bạn tham khảo:

**Example 1: Extract compounds with specific activity range**

*   **SQL (English):**

```sql
SELECT molregno, compound_structures.canonical_smiles, activities.standard_value
FROM compound_structures
JOIN activities ON compound_structures.molregno = activities.molregno
WHERE activities.standard_type = 'IC50'
AND activities.standard_units = 'nM'
AND activities.standard_value >= 100 AND activities.standard_value <= 1000
LIMIT 100;
```

*   **SQL (Tiếng Việt):**

```sql
SELECT molregno, compound_structures.canonical_smiles, activities.standard_value
FROM compound_structures
JOIN activities ON compound_structures.molregno = activities.molregno
WHERE activities.standard_type = 'IC50'
AND activities.standard_units = 'nM'
AND activities.standard_value >= 100 AND activities.standard_value <= 1000
LIMIT 100;
```

*   **Python (English):**

```python
# Read the CSV file
df = pd.read_csv(os.path.join(base_path, "activity_range.csv"))

# Print the number of rows and columns
print(f"Number of rows: {len(df)}")
print(f"Number of columns: {len(df.columns)}")

# Print the first 5 rows of the DataFrame
print(df.head())
```

*   **Python (Tiếng Việt):**

```python
# Đọc file CSV
df = pd.read_csv(os.path.join(base_path, "activity_range.csv"))

# In số lượng hàng và cột
print(f"Số lượng hàng: {len(df)}")
print(f"Số lượng cột: {len(df.columns)}")

# In 5 hàng đầu tiên của DataFrame
print(df.head())
```

**Example 2: Calculate Lipinski's Rule of Five**

*   **Python (English):**

```python
def lipinski_rule_of_five(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        # Check Lipinski's rule of five
        rule_1 = mw <= 500
        rule_2 = logp <= 5
        rule_3 = hbd <= 5
        rule_4 = hba <= 10

        # Count number of rules that fail
        num_failed_rules = sum([not rule_1, not rule_2, not rule_3, not rule_4])
        return num_failed_rules
    else:
        return None

# Apply Lipinski's Rule of Five to each SMILES
df['Lipinski_Failures'] = df['canonical_smiles'].apply(lipinski_rule_of_five)

# Print the results
print(df[['canonical_smiles', 'Lipinski_Failures']].head())
```

*   **Python (Tiếng Việt):**

```python
def lipinski_rule_of_five(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        # Kiểm tra quy tắc 5 của Lipinski
        rule_1 = mw <= 500
        rule_2 = logp <= 5
        rule_3 = hbd <= 5
        rule_4 = hba <= 10

        # Đếm số lượng quy tắc không đạt
        num_failed_rules = sum([not rule_1, not rule_2, not rule_3, not rule_4])
        return num_failed_rules
    else:
        return None

# Áp dụng quy tắc 5 của Lipinski cho mỗi SMILES
df['Lipinski_Failures'] = df['canonical_smiles'].apply(lipinski_rule_of_five)

# In kết quả
print(df[['canonical_smiles', 'Lipinski_Failures']].head())
```

**Example 3: Extract compounds based on substructure**

*   **SQL (English):**

```sql
SELECT molregno, compound_structures.canonical_smiles
FROM compound_structures
WHERE compound_structures.canonical_smiles LIKE '%C=O%'  -- Example: compounds containing a carbonyl group
LIMIT 100;
```

*   **SQL (Tiếng Việt):**

```sql
SELECT molregno, compound_structures.canonical_smiles
FROM compound_structures
WHERE compound_structures.canonical_smiles LIKE '%C=O%'  -- Ví dụ: hợp chất chứa nhóm carbonyl
LIMIT 100;
```

*   **Python (English):**

```python
from rdkit import Chem

def check_substructure(smiles, substructure_smarts):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        substructure = Chem.MolFromSmarts(substructure_smarts)
        if substructure is not None:
            return mol.HasSubstructMatch(substructure)
        else:
            return False
    else:
        return False

# Define substructure SMARTS
substructure_smarts = "C=O"  # Carbonyl group

# Apply the function to each SMILES
df['Has_Carbonyl'] = df['canonical_smiles'].apply(lambda x: check_substructure(x, substructure_smarts))

# Print the results
print(df[['canonical_smiles', 'Has_Carbonyl']].head())
```

*   **Python (Tiếng Việt):**

```python
from rdkit import Chem

def check_substructure(smiles, substructure_smarts):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        substructure = Chem.MolFromSmarts(substructure_smarts)
        if substructure is not None:
            return mol.HasSubstructMatch(substructure)
        else:
            return False
    else:
        return False

# Xác định cấu trúc con SMARTS
substructure_smarts = "C=O"  # Nhóm carbonyl

# Áp dụng hàm cho mỗi SMILES
df['Has_Carbonyl'] = df['canonical_smiles'].apply(lambda x: check_substructure(x, substructure_smarts))

# In kết quả
print(df[['canonical_smiles', 'Has_Carbonyl']].head())
```

**Example 4: Calculate Molecular Weight**

*   **Python (English):**

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None

# Apply the function to each SMILES
df['Molecular_Weight'] = df['canonical_smiles'].apply(calculate_molecular_weight)

# Print the results
print(df[['canonical_smiles', 'Molecular_Weight']].head())
```

*   **Python (Tiếng Việt):**

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None

# Áp dụng hàm cho mỗi SMILES
df['Molecular_Weight'] = df['canonical_smiles'].apply(calculate_molecular_weight)

# In kết quả
print(df[['canonical_smiles', 'Molecular_Weight']].head())
```

**Example 5: Visualize Molecules**

*   **Python (English):**

```python
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Draw.MolToImage(mol)
    else:
        return None

# Visualize the first molecule
img = visualize_molecule(df['canonical_smiles'].iloc[0])
img
```

*   **Python (Tiếng Việt):**

```python
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Draw.MolToImage(mol)
    else:
        return None

# Visualize phân tử đầu tiên
img = visualize_molecule(df['canonical_smiles'].iloc[0])
img
```

**5. Giải quyết lỗi (Error Resolution)**

*   **ERROR: operator does not exist: numeric ~ unknown LINE 12: AND act.standard_value ~ '^[0-9\.]+$'**

    *   **Explanation (Giải thích):** Lỗi này xảy ra vì toán tử `~` (regex match) không được hỗ trợ trực tiếp trên kiểu dữ liệu `numeric` trong PostgreSQL.
    *   **Solution (Giải pháp):** Chuyển đổi cột `standard_value` sang kiểu `text` trước khi so sánh với biểu thức chính quy.
    *   **SQL (English):**

```sql
SELECT
    molregno,
    compound_structures.canonical_smiles,
    activities.standard_value,
    activities.standard_units
FROM
    compound_structures
JOIN
    activities ON compound_structures.molregno = activities.molregno
JOIN
    target_dictionary ON activities.tid = target_dictionary.tid
WHERE
    target_dictionary.pref_name = 'desired_target_name'  -- Replace with your target
    AND activities.standard_type = 'IC50'
    AND activities.standard_units = 'nM'
    AND CAST(activities.standard_value AS TEXT) ~ '^[0-9\.]+$'  -- Convert to TEXT
LIMIT 100;
```

    *   **SQL (Tiếng Việt):**

```sql
SELECT
    molregno,
    compound_structures.canonical_smiles,
    activities.standard_value,
    activities.standard_units
FROM
    compound_structures
JOIN
    activities ON compound_structures.molregno = activities.molregno
JOIN
    target_dictionary ON activities.tid = target_dictionary.tid
WHERE
    target_dictionary.pref_name = 'desired_target_name'  -- Thay thế bằng mục tiêu của bạn
    AND activities.standard_type = 'IC50'
    AND activities.standard_units = 'nM'
    AND CAST(activities.standard_value AS TEXT) ~ '^[0-9\.]+$'  -- Chuyển đổi sang TEXT
LIMIT 100;
```

*   **Phiên bản scikit-learn cũ không hỗ trợ tham số `squared=False` trong hàm `mean_squared_error`**

    *   **Explanation (Giải thích):** Tham số `squared=False` được thêm vào hàm `mean_squared_error` trong các phiên bản scikit-learn mới hơn để trả về Root Mean Squared Error (RMSE) thay vì Mean Squared Error (MSE).
    *   **Solution (Giải pháp):**
        1.  **Update scikit-learn (Cập nhật scikit-learn):** Nâng cấp phiên bản scikit-learn của bạn lên phiên bản mới nhất bằng lệnh `pip install --upgrade scikit-learn`.
        2.  **Calculate RMSE manually (Tính RMSE thủ công):** Nếu bạn không thể nâng cấp scikit-learn, bạn có thể tính RMSE bằng cách lấy căn bậc hai của MSE.

        *   **Python (English):**

```python
from sklearn.metrics import mean_squared_error
import numpy as np

# Calculate MSE
mse = mean_squared_error(y_test, y_pred)

# Calculate RMSE
rmse = np.sqrt(mse)

print(f"Mean Squared Error: {mse}")
print(f"Root Mean Squared Error: {rmse}")
```

        *   **Python (Tiếng Việt):**

```python
from sklearn.metrics import mean_squared_error
import numpy as np

# Tính MSE
mse = mean_squared_error(y_test, y_pred)

# Tính RMSE
rmse = np.sqrt(mse)

print(f"Mean Squared Error: {mse}")
print(f"Root Mean Squared Error: {rmse}")
```

**Lưu ý quan trọng:**

*   Hãy thay thế `"your_data.csv"` và `"desired_target_name"` bằng tên tệp và mục tiêu cụ thể của bạn.
*   Đảm bảo rằng bạn đã cài đặt tất cả các thư viện cần thiết (RDKit, pandas, scikit-learn, etc.) trước khi chạy code Python.
*   Khi làm việc với dữ liệu lớn, hãy cân nhắc sử dụng các kỹ thuật tối ưu hóa để cải thiện hiệu suất.

Chúc bạn thành công với dự án phân tích dữ liệu ChEMBL 35! Nếu bạn có bất kỳ câu hỏi nào khác, đừng ngần ngại hỏi.
