# PharmApp Suite
## 🧠 AI for Drug Discovery and Development 🧪
| Copyright 2025 | Nghiên Cứu Thuốc | www.nghiencuuthuoc.com | Zalo: +84888999311 |

🧩 Topic: Phân tích ligand nhắm GPCR
---
Tuyệt vời! Tôi rất vui được tham gia vào dự án phân tích dữ liệu ChEMBL 35 của bạn. Dựa trên thông tin bạn cung cấp và những yêu cầu cụ thể, tôi sẽ giúp bạn xây dựng quy trình phân tích dữ liệu hiệu quả, tập trung vào việc giải quyết các vấn đề bạn gặp phải và tuân thủ theo template AIMLOps.

**Tổng quan về quy trình phân tích**

Chúng ta sẽ thực hiện các bước sau:

1.  **Data Extraction (Trích xuất dữ liệu):** Sử dụng SQL để truy vấn dữ liệu từ cơ sở dữ liệu ChEMBL 35 và lưu vào file CSV. Chú ý xử lý lỗi `operator does not exist: numeric ~ unknown`.
2.  **Data Preprocessing (Tiền xử lý dữ liệu):** Đọc file CSV, làm sạch và chuẩn hóa dữ liệu bằng RDKit và các thư viện Python khác.
3.  **Feature Engineering (Xây dựng đặc trưng):** Tạo các đặc trưng hóa học (chemical features) từ SMILES strings sử dụng RDKit.
4.  **Model Building (Xây dựng mô hình):** Xây dựng mô hình học máy (ví dụ: regression) để dự đoán hoạt tính của các hợp chất.
5.  **Model Evaluation (Đánh giá mô hình):** Đánh giá hiệu năng của mô hình và tinh chỉnh nếu cần.

**1. Phân tích mô hình (Model Analysis)**

Mô hình phân tích sẽ tập trung vào việc xây dựng một mô hình dự đoán hoạt tính (activity) của các hợp chất dựa trên cấu trúc hóa học của chúng. Cụ thể, chúng ta sẽ sử dụng mô hình hồi quy (regression model) để dự đoán giá trị hoạt tính (ví dụ: IC50, Ki).

*   **Input (Đầu vào):** SMILES strings (biểu diễn cấu trúc hóa học) của các hợp chất.
*   **Process (Quy trình):**
    *   Sử dụng RDKit để chuyển đổi SMILES strings thành các đặc trưng số (numerical features) như fingerprints, descriptors.
    *   Sử dụng các thuật toán hồi quy (ví dụ: Linear Regression, Random Forest, XGBoost) để xây dựng mô hình dự đoán.
*   **Output (Đầu ra):** Mô hình dự đoán giá trị hoạt tính của hợp chất dựa trên cấu trúc hóa học.

**2. Hướng dẫn song ngữ (Bilingual Instructions)**

**2.1. SQL Code (Mã SQL)**

*   **Purpose:** Extract data from ChEMBL database and save to CSV file. (Mục đích: Trích xuất dữ liệu từ cơ sở dữ liệu ChEMBL và lưu vào file CSV.)
*   **Explanation:** The following SQL code retrieves compound information and activity data, addressing the error related to numeric comparison. (Giải thích: Mã SQL sau truy xuất thông tin hợp chất và dữ liệu hoạt tính, giải quyết lỗi liên quan đến so sánh số.)

```sql
-- English
-- Extracting data from ChEMBL database

SELECT
    md.molregno,
    cs.canonical_smiles,
    act.standard_type,
    act.standard_value,
    act.standard_units
FROM
    molecule_dictionary md
JOIN
    compound_structures cs ON md.molregno = cs.molregno
JOIN
    activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
  AND act.standard_units = 'nM'
  AND act.standard_value IS NOT NULL
  AND act.standard_value::text ~ '^[0-9\.]+$' -- Check if the value is numeric
LIMIT 100;

-- Vietnamese
-- Trích xuất dữ liệu từ cơ sở dữ liệu ChEMBL

SELECT
    md.molregno,
    cs.canonical_smiles,
    act.standard_type,
    act.standard_value,
    act.standard_units
FROM
    molecule_dictionary md
JOIN
    compound_structures cs ON md.molregno = cs.molregno
JOIN
    activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
  AND act.standard_units = 'nM'
  AND act.standard_value IS NOT NULL
  AND act.standard_value::text ~ '^[0-9\.]+$' -- Kiểm tra xem giá trị có phải là số không
LIMIT 100;
```

**Giải thích lỗi và cách khắc phục:**

Lỗi `ERROR: operator does not exist: numeric ~ unknown` xảy ra do PostgreSQL không thể so sánh trực tiếp kiểu dữ liệu `numeric` với một biểu thức regular expression (kiểu `unknown`). Để khắc phục, chúng ta cần ép kiểu dữ liệu `standard_value` về kiểu `text` trước khi so sánh bằng `~` (regular expression matching).

**2.2. Python Code (Mã Python)**

*   **Purpose:** Load data from CSV, preprocess data using RDKit, and build a regression model. (Mục đích: Tải dữ liệu từ CSV, tiền xử lý dữ liệu bằng RDKit và xây dựng mô hình hồi quy.)
*   **Explanation:** The following Python code loads data, converts SMILES to features, and builds a linear regression model. It also addresses the scikit-learn version issue. (Giải thích: Mã Python sau tải dữ liệu, chuyển đổi SMILES thành đặc trưng và xây dựng mô hình hồi quy tuyến tính. Nó cũng giải quyết vấn đề phiên bản scikit-learn.)

```python
# English
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

# Define base path
base_path = "../data"
csv_file = "chembl_data.csv" # Replace with your actual file name
csv_path = os.path.join(base_path, csv_file)

# Load data from CSV
try:
    df = pd.read_csv(csv_path)
except FileNotFoundError:
    print(f"Error: The file {csv_path} was not found.")
    exit()

# Data Preprocessing
df = df.dropna(subset=['canonical_smiles', 'standard_value'])

# Convert SMILES to features using RDKit (Morgan Fingerprints)
def smiles_to_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        return list(fp)
    else:
        return None

df['features'] = df['canonical_smiles'].apply(smiles_to_features)
df = df.dropna(subset=['features'])

# Prepare data for modeling
X = list(df['features'])
y = df['standard_value'].astype(float)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Build and train the Linear Regression model
model = LinearRegression()
model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
try:
    mse = mean_squared_error(y_test, y_pred, squared=False) # Try with squared=False first
except TypeError:
    mse = mean_squared_error(y_test, y_pred) # If squared=False is not supported

print(f"Mean Squared Error: {mse}")

# Vietnamese
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

# Định nghĩa đường dẫn cơ sở
base_path = "../data"
csv_file = "chembl_data.csv" # Thay thế bằng tên file thực tế của bạn
csv_path = os.path.join(base_path, csv_file)

# Tải dữ liệu từ CSV
try:
    df = pd.read_csv(csv_path)
except FileNotFoundError:
    print(f"Lỗi: Không tìm thấy file {csv_path}.")
    exit()

# Tiền xử lý dữ liệu
df = df.dropna(subset=['canonical_smiles', 'standard_value'])

# Chuyển đổi SMILES thành đặc trưng sử dụng RDKit (Morgan Fingerprints)
def smiles_to_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        return list(fp)
    else:
        return None

df['features'] = df['canonical_smiles'].apply(smiles_to_features)
df = df.dropna(subset=['features'])

# Chuẩn bị dữ liệu cho mô hình
X = list(df['features'])
y = df['standard_value'].astype(float)

# Chia dữ liệu thành tập huấn luyện và tập kiểm tra
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Xây dựng và huấn luyện mô hình hồi quy tuyến tính
model = LinearRegression()
model.fit(X_train, y_train)

# Dự đoán
y_pred = model.predict(X_test)

# Đánh giá mô hình
try:
    mse = mean_squared_error(y_test, y_pred, squared=False) # Thử với squared=False trước
except TypeError:
    mse = mean_squared_error(y_test, y_pred) # Nếu squared=False không được hỗ trợ

print(f"Sai số bình phương trung bình: {mse}")
```

**Giải thích lỗi và cách khắc phục:**

*   **`squared=False` in `mean_squared_error`:**  Phiên bản cũ của scikit-learn không hỗ trợ tham số `squared=False` trong hàm `mean_squared_error`. Để khắc phục, chúng ta sử dụng một khối `try...except` để thử sử dụng `squared=False` trước. Nếu lỗi `TypeError` xảy ra, chúng ta sẽ sử dụng `mean_squared_error` mà không có tham số `squared=False`.
*   **Xử lý `FileNotFoundError`:** Đoạn code đã được bổ sung khối `try...except` để bắt lỗi `FileNotFoundError` khi file CSV không tồn tại.

**3. Ví dụ Code (Code Examples)**

Dưới đây là 5 ví dụ code SQL và Python mẫu để bạn tham khảo:

**3.1. SQL Examples**

1.  **Extract data for a specific target (Lấy dữ liệu cho một mục tiêu cụ thể):**

```sql
-- English
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
JOIN target_dictionary td ON act.tid = td.tid
WHERE td.target_name = 'CHEMBL205' -- Replace with your target name
AND act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;

-- Vietnamese
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
JOIN target_dictionary td ON act.tid = td.tid
WHERE td.target_name = 'CHEMBL205' -- Thay thế bằng tên mục tiêu của bạn
AND act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;
```

2.  **Extract data for a specific assay (Lấy dữ liệu cho một assay cụ thể):**

```sql
-- English
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
JOIN assays a ON act.assay_id = a.assay_id
WHERE a.description LIKE '%acetylcholinesterase%' -- Replace with your assay description
AND act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;

-- Vietnamese
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
JOIN assays a ON act.assay_id = a.assay_id
WHERE a.description LIKE '%acetylcholinesterase%' -- Thay thế bằng mô tả assay của bạn
AND act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;
```

3.  **Extract data with a specific activity range (Lấy dữ liệu với một khoảng giá trị hoạt tính cụ thể):**

```sql
-- English
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value BETWEEN 10 AND 100 -- Replace with your activity range
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;

-- Vietnamese
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value BETWEEN 10 AND 100 -- Thay thế bằng khoảng giá trị hoạt tính của bạn
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;
```

4.  **Extract data and order by activity (Lấy dữ liệu và sắp xếp theo hoạt tính):**

```sql
-- English
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
ORDER BY act.standard_value ASC -- or DESC for descending order
LIMIT 100;

-- Vietnamese
SELECT md.molregno, cs.canonical_smiles, act.standard_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
ORDER BY act.standard_value ASC -- hoặc DESC để sắp xếp giảm dần
LIMIT 100;
```

5.  **Extract distinct compounds (Lấy các hợp chất khác nhau):**

```sql
-- English
SELECT DISTINCT cs.canonical_smiles
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;

-- Vietnamese
SELECT DISTINCT cs.canonical_smiles
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
AND act.standard_units = 'nM'
AND act.standard_value IS NOT NULL
AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;
```

**3.2. Python Examples**

1.  **Calculate and visualize molecular weight (Tính toán và trực quan hóa khối lượng phân tử):**

```python
# English
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_mw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None

df['mol_weight'] = df['canonical_smiles'].apply(calculate_mw)
df = df.dropna(subset=['mol_weight'])

plt.hist(df['mol_weight'], bins=50)
plt.xlabel('Molecular Weight')
plt.ylabel('Frequency')
plt.title('Distribution of Molecular Weight')
plt.show()

# Vietnamese
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_mw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None

df['mol_weight'] = df['canonical_smiles'].apply(calculate_mw)
df = df.dropna(subset=['mol_weight'])

plt.hist(df['mol_weight'], bins=50)
plt.xlabel('Khối lượng phân tử')
plt.ylabel('Tần số')
plt.title('Phân bố khối lượng phân tử')
plt.show()
```

2.  **Calculate LogP (Tính toán LogP):**

```python
# English
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolLogP(mol)
    else:
        return None

df['logp'] = df['canonical_smiles'].apply(calculate_logp)
df = df.dropna(subset=['logp'])
print(df[['canonical_smiles', 'logp']].head())

# Vietnamese
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolLogP(mol)
    else:
        return None

df['logp'] = df['canonical_smiles'].apply(calculate_logp)
df = df.dropna(subset=['logp'])
print(df[['canonical_smiles', 'logp']].head())
```

3.  **Calculate TPSA (Tính toán TPSA):**

```python
# English
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_tpsa(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.TPSA(mol)
    else:
        return None

df['tpsa'] = df['canonical_smiles'].apply(calculate_tpsa)
df = df.dropna(subset=['tpsa'])
print(df[['canonical_smiles', 'tpsa']].head())

# Vietnamese
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_tpsa(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.TPSA(mol)
    else:
        return None

df['tpsa'] = df['canonical_smiles'].apply(calculate_tpsa)
df = df.dropna(subset=['tpsa'])
print(df[['canonical_smiles', 'tpsa']].head())
```

4.  **Calculate the number of hydrogen bond donors and acceptors (Tính toán số lượng chất cho và nhận liên kết hydro):**

```python
# English
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_hbd(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.NumHDonors(mol)
    else:
        return None

def calculate_hba(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.NumHAcceptors(mol)
    else:
        return None

df['hbd'] = df['canonical_smiles'].apply(calculate_hbd)
df['hba'] = df['canonical_smiles'].apply(calculate_hba)

df = df.dropna(subset=['hbd', 'hba'])
print(df[['canonical_smiles', 'hbd', 'hba']].head())

# Vietnamese
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

base_path = "../data"
csv_file = "chembl_data.csv"
csv_path = os.path.join(base_path, csv_file)

df = pd.read_csv(csv_path)
df = df.dropna(subset=['canonical_smiles'])

def calculate_hbd(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.NumHDonors(mol)
    else:
        return None

def calculate_hba(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.NumHAcceptors(mol)
    else:
        return None

df['hbd'] = df['canonical_smiles'].apply(calculate_hbd)
df['hba'] = df['canonical_smiles'].apply(calculate_hba)

df = df.dropna(subset=['hbd', 'hba'])
print(df[['canonical_smiles', 'hbd', 'hba']].head())
```

5.  **Train a Random Forest Regressor model (Huấn luyện mô hình Random Forest Regressor):**

```python
# English
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# Define base path
base_path = "../data"
csv_file = "chembl_data.csv" # Replace with your actual file name
csv_path = os.path.join(base_path, csv_file)

# Load data from CSV
try:
    df = pd.read_csv(csv_path)
except FileNotFoundError:
    print(f"Error: The file {csv_path} was not found.")
    exit()

# Data Preprocessing
df = df.dropna(subset=['canonical_smiles', 'standard_value'])

# Convert SMILES to features using RDKit (Morgan Fingerprints)
def smiles_to_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        return list(fp)
    else:
        return None

df['features'] = df['canonical_smiles'].apply(smiles_to_features)
df = df.dropna(subset=['features'])

# Prepare data for modeling
X = list(df['features'])
y = df['standard_value'].astype(float)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Build and train the Random Forest Regressor model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)

print(f"Mean Squared Error: {mse}")

# Vietnamese
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# Định nghĩa đường dẫn cơ sở
base_path = "../data"
csv_file = "chembl_data.csv" # Thay thế bằng tên file thực tế của bạn
csv_path = os.path.join(base_path, csv_file)

# Tải dữ liệu từ CSV
try:
    df = pd.read_csv(csv_path)
except FileNotFoundError:
    print(f"Lỗi: Không tìm thấy file {csv_path}.")
    exit()

# Tiền xử lý dữ liệu
df = df.dropna(subset=['canonical_smiles', 'standard_value'])

# Chuyển đổi SMILES thành đặc trưng sử dụng RDKit (Morgan Fingerprints)
def smiles_to_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        return list(fp)
    else:
        return None

df['features'] = df['canonical_smiles'].apply(smiles_to_features)
df = df.dropna(subset=['features'])

# Chuẩn bị dữ liệu cho mô hình
X = list(df['features'])
y = df['standard_value'].astype(float)

# Chia dữ liệu thành tập huấn luyện và tập kiểm tra
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Xây dựng và huấn luyện mô hình Random Forest Regressor
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Dự đoán
y_pred = model.predict(X_test)

# Đánh giá mô hình
mse = mean_squared_error(y_test, y_pred)

print(f"Sai số bình phương trung bình: {mse}")
```

**Lưu ý quan trọng:**

*   Hãy thay thế `"chembl_data.csv"` bằng tên file CSV thực tế của bạn.
*   Đảm bảo bạn đã cài đặt tất cả các thư viện cần thiết (pandas, rdkit, scikit-learn, matplotlib).
*   Điều chỉnh các tham số của mô hình (ví dụ: `n_estimators` trong `RandomForestRegressor`) để đạt được hiệu năng tốt nhất.

Tôi hy vọng những hướng dẫn và ví dụ này sẽ giúp bạn tiến xa hơn trong dự án của mình. Nếu bạn có bất kỳ câu hỏi nào, đừng ngần ngại hỏi!
