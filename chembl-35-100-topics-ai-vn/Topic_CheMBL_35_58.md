# PharmApp Suite
## 🧠 AI for Drug Discovery and Development 🧪
| Copyright 2025 | Nghiên Cứu Thuốc | www.nghiencuuthuoc.com | Zalo: +84888999311 |

🧩 Topic: Phân tích thuốc dành cho trẻ em
---
Tuyệt vời! Dự án phân tích dữ liệu ChEMBL 35 bằng RDKit là một hướng đi rất tiềm năng trong nghiên cứu và phát triển thuốc. Với vai trò là chuyên gia R&D Dược học, tôi sẽ hỗ trợ bạn xây dựng quy trình phân tích hiệu quả, khắc phục các lỗi và đưa ra các ví dụ code cụ thể.

**Tổng quan về chủ đề Topic_CheMBL_35_58**

Chủ đề này có vẻ tập trung vào một khía cạnh cụ thể trong ChEMBL 35, có thể là một nhóm hợp chất, một loại mục tiêu sinh học (target), hoặc một tính chất hóa lý nào đó. Để đưa ra phân tích chính xác hơn, tôi cần biết tên chủ đề đầy đủ hoặc mô tả chi tiết hơn về mục tiêu của bạn. Tuy nhiên, dựa trên thông tin bạn cung cấp, tôi sẽ xây dựng một quy trình phân tích tổng quát, có thể áp dụng và điều chỉnh cho nhiều chủ đề khác nhau trong ChEMBL 35.

**1. Phân tích mô hình phân tích (Analysis Model)**

Mô hình phân tích của chúng ta sẽ bao gồm các bước sau:

1.  **Data Extraction (Trích xuất dữ liệu):** Sử dụng SQL để truy vấn dữ liệu từ cơ sở dữ liệu ChEMBL 35, tập trung vào các bảng quan trọng như `activities`, `molecule_dictionary`, `compound_structures`, và `target_dictionary`.
2.  **Data Cleaning and Preprocessing (Làm sạch và tiền xử lý dữ liệu):**
    *   Loại bỏ các giá trị null hoặc không hợp lệ.
    *   Chuẩn hóa dữ liệu (ví dụ: chuyển đổi IC50, Ki về một đơn vị chung).
    *   Lọc dữ liệu dựa trên các tiêu chí cụ thể (ví dụ: chỉ giữ lại các hoạt chất có IC50 < 10µM).
3.  **Feature Engineering (Xây dựng đặc trưng):** Sử dụng RDKit để tính toán các đặc trưng phân tử (molecular descriptors) từ cấu trúc hóa học của các hợp chất. Các đặc trưng này có thể bao gồm:
    *   Tính chất vật lý (ví dụ: LogP, MW, TPSA).
    *   Hình dạng phân tử (ví dụ: số vòng, số nhánh).
    *   Các đặc trưng dựa trên fingerprint (ví dụ: Morgan fingerprint, MACCS keys).
4.  **Data Analysis and Modeling (Phân tích dữ liệu và xây dựng mô hình):**
    *   **Exploratory Data Analysis (EDA):** Thống kê mô tả, trực quan hóa dữ liệu để hiểu rõ hơn về phân bố và mối quan hệ giữa các biến.
    *   **Machine Learning (ML):** Xây dựng các mô hình dự đoán hoạt tính sinh học (ví dụ: sử dụng Random Forest, Support Vector Machines) hoặc phân loại hợp chất (ví dụ: dựa trên cấu trúc hoặc hoạt tính).
5.  **Model Evaluation and Interpretation (Đánh giá và giải thích mô hình):** Đánh giá hiệu năng của mô hình bằng các metrics phù hợp (ví dụ: AUC, ROC, RMSE) và giải thích các yếu tố quan trọng ảnh hưởng đến hoạt tính.

**2. Hướng dẫn song ngữ (Bilingual Instructions)**

Dưới đây là hướng dẫn chi tiết cho từng bước, kèm theo giải thích bằng cả tiếng Anh và tiếng Việt:

**Step 1: Data Extraction (Trích xuất dữ liệu)**

*   **English:** Use SQL to query data from the ChEMBL 35 database, focusing on tables like `activities`, `molecule_dictionary`, `compound_structures`, and `target_dictionary`.
*   **Tiếng Việt:** Sử dụng SQL để truy vấn dữ liệu từ cơ sở dữ liệu ChEMBL 35, tập trung vào các bảng như `activities`, `molecule_dictionary`, `compound_structures`, và `target_dictionary`.

**Step 2: Data Cleaning and Preprocessing (Làm sạch và tiền xử lý dữ liệu)**

*   **English:**
    *   Remove null or invalid values.
    *   Standardize data (e.g., convert IC50, Ki to a common unit).
    *   Filter data based on specific criteria (e.g., keep only compounds with IC50 < 10µM).
*   **Tiếng Việt:**
    *   Loại bỏ các giá trị null hoặc không hợp lệ.
    *   Chuẩn hóa dữ liệu (ví dụ: chuyển đổi IC50, Ki về một đơn vị chung).
    *   Lọc dữ liệu dựa trên các tiêu chí cụ thể (ví dụ: chỉ giữ lại các hợp chất có IC50 < 10µM).

**Step 3: Feature Engineering (Xây dựng đặc trưng)**

*   **English:** Use RDKit to calculate molecular descriptors from the chemical structures of the compounds. These descriptors can include:
    *   Physicochemical properties (e.g., LogP, MW, TPSA).
    *   Molecular shape (e.g., number of rings, number of branches).
    *   Fingerprint-based features (e.g., Morgan fingerprint, MACCS keys).
*   **Tiếng Việt:** Sử dụng RDKit để tính toán các đặc trưng phân tử từ cấu trúc hóa học của các hợp chất. Các đặc trưng này có thể bao gồm:
    *   Tính chất vật lý (ví dụ: LogP, MW, TPSA).
    *   Hình dạng phân tử (ví dụ: số vòng, số nhánh).
    *   Các đặc trưng dựa trên fingerprint (ví dụ: Morgan fingerprint, MACCS keys).

**Step 4: Data Analysis and Modeling (Phân tích dữ liệu và xây dựng mô hình)**

*   **English:**
    *   **Exploratory Data Analysis (EDA):** Descriptive statistics, data visualization to understand the distribution and relationships between variables.
    *   **Machine Learning (ML):** Build models to predict biological activity (e.g., using Random Forest, Support Vector Machines) or classify compounds (e.g., based on structure or activity).
*   **Tiếng Việt:**
    *   **Exploratory Data Analysis (EDA):** Thống kê mô tả, trực quan hóa dữ liệu để hiểu rõ hơn về phân bố và mối quan hệ giữa các biến.
    *   **Machine Learning (ML):** Xây dựng các mô hình dự đoán hoạt tính sinh học (ví dụ: sử dụng Random Forest, Support Vector Machines) hoặc phân loại hợp chất (ví dụ: dựa trên cấu trúc hoặc hoạt tính).

**Step 5: Model Evaluation and Interpretation (Đánh giá và giải thích mô hình)**

*   **English:** Evaluate the performance of the model using appropriate metrics (e.g., AUC, ROC, RMSE) and interpret the important factors affecting activity.
*   **Tiếng Việt:** Đánh giá hiệu năng của mô hình bằng các metrics phù hợp (ví dụ: AUC, ROC, RMSE) và giải thích các yếu tố quan trọng ảnh hưởng đến hoạt tính.

**3. Code SQL, Python (English)**

Dưới đây là các ví dụ code SQL và Python, tập trung vào việc trích xuất dữ liệu, xây dựng đặc trưng và phân tích cơ bản.

**SQL (PostgreSQL)**

```sql
-- Example 1: Extracting data for a specific target (e.g., CHEMBL205)
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
WHERE
    act.target_chembl_id = 'CHEMBL205'  -- Replace with your target CHEMBL_ID
    AND act.standard_type = 'IC50'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value > 0
    AND act.standard_value < 10000 -- Filter IC50 values
    AND cs.canonical_smiles IS NOT NULL
LIMIT 100;

-- Example 2:  Corrected SQL to filter numeric values for standard_value
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
WHERE
    act.target_chembl_id = 'CHEMBL205'
    AND act.standard_type = 'IC50'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'  -- Ensure standard_value is numeric
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    AND cs.canonical_smiles IS NOT NULL
LIMIT 100;

-- Example 3: Extracting data for a specific activity type (e.g., Ki)
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
WHERE
    act.standard_type = 'Ki'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    AND cs.canonical_smiles IS NOT NULL
LIMIT 100;

-- Example 4: Extracting data with specific units (e.g., nM)
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
WHERE
    act.standard_type = 'IC50'
    AND act.standard_units = 'nM'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    AND cs.canonical_smiles IS NOT NULL
LIMIT 100;

-- Example 5: Extracting data for a specific protein target type
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
JOIN
    target_dictionary td ON act.tid = td.tid
WHERE
    td.target_type = 'PROTEIN'
    AND act.standard_type = 'IC50'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    AND cs.canonical_smiles IS NOT NULL
LIMIT 100;
```

**Python (Jupyter Notebook)**

```python
import os
import pandas as pd
import psycopg2
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Database credentials
DB_HOST = "192.168.206.136"
DB_USER = "rd"
DB_PASS = "rd"
DB_NAME = "chembl_35"
base_path = "." # Assuming the notebook is in the root directory

# Function to connect to the database
def connect_to_db(host, user, password, database):
    conn = psycopg2.connect(host=host, user=user, password=password, database=database)
    return conn

# Function to execute SQL query and return a Pandas DataFrame
def execute_sql_query(conn, query):
    df = pd.read_sql_query(query, conn)
    return df

# Function to calculate molecular descriptors using RDKit
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    descriptors = {}
    descriptors['MW'] = Descriptors.MolWt(mol)
    descriptors['LogP'] = Descriptors.MolLogP(mol)
    # Add more descriptors as needed
    return descriptors

# Main function
def main():
    # Connect to the database
    conn = connect_to_db(DB_HOST, DB_USER, DB_PASS, DB_NAME)

    # SQL query to extract data
    sql_query = """
        SELECT
            act.standard_value,
            act.standard_units,
            md.chembl_id,
            cs.canonical_smiles
        FROM
            activities act
        JOIN
            molecule_dictionary md ON act.molregno = md.molregno
        JOIN
            compound_structures cs ON md.molregno = cs.molregno
        WHERE
            act.target_chembl_id = 'CHEMBL205'  -- Replace with your target CHEMBL_ID
            AND act.standard_type = 'IC50'
            AND act.standard_relation = '='
            AND act.standard_value IS NOT NULL
            AND act.standard_value::TEXT ~ '^[0-9\.]+$'  -- Ensure standard_value is numeric
            AND act.standard_value::NUMERIC > 0
            AND act.standard_value::NUMERIC < 10000
            AND cs.canonical_smiles IS NOT NULL
        LIMIT 100;
    """

    # Execute the SQL query
    df = execute_sql_query(conn, sql_query)
    conn.close()

    # Data Cleaning and Preprocessing
    df = df.dropna(subset=['canonical_smiles', 'standard_value']) # Drop rows with missing SMILES or standard_value
    df['standard_value'] = pd.to_numeric(df['standard_value']) # Ensure standard_value is numeric

    # Feature Engineering
    descriptors_list = []
    for smiles in df['canonical_smiles']:
        descriptors = calculate_descriptors(smiles)
        descriptors_list.append(descriptors)

    descriptors_df = pd.DataFrame(descriptors_list)
    df = pd.concat([df, descriptors_df], axis=1)
    df = df.dropna()  # Drop rows with missing descriptors

    # Data Analysis and Modeling
    X = df[['MW', 'LogP']]  # Features
    y = np.log10(df['standard_value'])  # Target variable (log transform for better distribution)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    # Model Evaluation
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    print(f"Mean Squared Error: {mse}")
    print(f"R-squared: {r2}")

    # Save the DataFrame to a CSV file (optional)
    output_file = os.path.join(base_path, "data", "chembl_data.csv")
    df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")

if __name__ == "__main__":
    main()
```

**Explanation:**

*   **SQL:** The SQL query extracts `standard_value`, `standard_units`, `chembl_id`, and `canonical_smiles` for a specific target from the ChEMBL database.  The `WHERE` clause filters the data based on the target ID, activity type (IC50), relation (=), and ensures that the `standard_value` is numeric, not null, and within a reasonable range. The `LIMIT` clause restricts the output to the first 100 rows.  The corrected SQL uses `::TEXT ~ '^[0-9\.]+$'` and `::NUMERIC` to properly handle numeric filtering.
*   **Python:**
    *   The Python code connects to the ChEMBL database using the provided credentials.
    *   It executes the SQL query and loads the data into a Pandas DataFrame.
    *   It calculates molecular descriptors (MW and LogP) using RDKit.
    *   It builds a Random Forest Regression model to predict the `standard_value` based on the calculated descriptors.
    *   It evaluates the model using Mean Squared Error (MSE) and R-squared.
    *   Finally, it saves the processed DataFrame to a CSV file.

**Fixing the Errors:**

*   **Error a: `ERROR: operator does not exist: numeric ~ unknown LINE 12: AND act.standard_value ~ '^[0-9\.]+$'`**

    This error occurs because PostgreSQL's `~` operator (regular expression match) is not directly applicable to numeric columns. To fix this, you need to cast the `standard_value` to text before using the `~` operator and cast to numeric for comparison:

    ```sql
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'  -- Ensure standard_value is numeric
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    ```

*   **Error b: `mean_squared_error` does not support `squared=False`**

    This error indicates that you are using an older version of scikit-learn.  Either upgrade scikit-learn to a version that supports `squared=False` (which returns the root mean squared error (RMSE) directly) or calculate the RMSE manually:

    ```python
    mse = mean_squared_error(y_test, y_pred)
    rmse = np.sqrt(mse) # Calculate RMSE manually
    print(f"Mean Squared Error: {mse}")
    print(f"Root Mean Squared Error: {rmse}") # Print RMSE
    ```

**4. Five Examples of SQL and Python Code**

Here are five more examples demonstrating different data extraction and analysis techniques:

**Example 1:  Filtering by Molecular Weight (SQL & Python)**

*   **SQL:** Extract compounds with molecular weight between 200 and 500.  This requires joining with a table that contains molecular weight information (which may require creating a custom table or function if not directly available in ChEMBL).  The example assumes a table called `molecule_properties` with a `molregno` and `molecular_weight` column.

```sql
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
LEFT JOIN
    molecule_properties mp ON md.molregno = mp.molregno -- Assuming a table with MW info
WHERE
    act.target_chembl_id = 'CHEMBL205'
    AND act.standard_type = 'IC50'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    AND cs.canonical_smiles IS NOT NULL
    AND mp.molecular_weight BETWEEN 200 AND 500 -- Filter by molecular weight
LIMIT 100;

```

*   **Python:** Add molecular weight filtering after extracting the data using RDKit.

```python
# ... (previous code) ...

    # Feature Engineering (including MW calculation)
    descriptors_list = []
    for smiles in df['canonical_smiles']:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            descriptors_list.append(None)
            continue
        descriptors = {}
        descriptors['MW'] = Descriptors.MolWt(mol)
        descriptors['LogP'] = Descriptors.MolLogP(mol)
        descriptors_list.append(descriptors)

    descriptors_df = pd.DataFrame(descriptors_list)
    df = pd.concat([df, descriptors_df], axis=1)
    df = df.dropna()

    # Filter by molecular weight
    df = df[(df['MW'] >= 200) & (df['MW'] <= 500)]

    # ... (rest of the code) ...
```

**Example 2:  Using Morgan Fingerprints (Python)**

*   **Python:**  Calculate Morgan fingerprints and use them as features in the model.

```python
from rdkit.Chem import AllChem

# Function to calculate Morgan Fingerprint
def calculate_morgan_fingerprint(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    return fp

# ... (previous code) ...

    # Feature Engineering (including Morgan Fingerprint)
    fingerprints = []
    for smiles in df['canonical_smiles']:
        fp = calculate_morgan_fingerprint(smiles)
        if fp is None:
            fingerprints.append(None)
        else:
            fingerprints.append(fp)

    df['fingerprint'] = fingerprints
    df = df.dropna(subset=['fingerprint'])

    # Convert fingerprints to a matrix for machine learning
    X = np.array([list(fp) for fp in df['fingerprint']])
    y = np.log10(df['standard_value'])

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    # ... (rest of the code) ...
```

**Example 3:  Extracting Data for Multiple Target CHEMBL IDs (SQL & Python)**

*   **SQL:** Extract data for several targets at once.

```sql
SELECT
    act.standard_value,
    act.standard_units,
    md.chembl_id,
    cs.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
WHERE
    act.target_chembl_id IN ('CHEMBL205', 'CHEMBL206', 'CHEMBL207') -- Multiple targets
    AND act.standard_type = 'IC50'
    AND act.standard_relation = '='
    AND act.standard_value IS NOT NULL
    AND act.standard_value::TEXT ~ '^[0-9\.]+$'
    AND act.standard_value::NUMERIC > 0
    AND act.standard_value::NUMERIC < 10000
    AND cs.canonical_smiles IS NOT NULL
LIMIT 100;
```

*   **Python:** No significant changes needed in the Python code, as it processes the data returned by the SQL query.

**Example 4:  Calculating TPSA (Topological Polar Surface Area) (Python)**

*   **Python:** Add TPSA calculation to the molecular descriptor calculation.

```python
from rdkit.Chem import rdMolDescriptors

# Function to calculate molecular descriptors including TPSA
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    descriptors = {}
    descriptors['MW'] = Descriptors.MolWt(mol)
    descriptors['LogP'] = Descriptors.MolLogP(mol)
    descriptors['TPSA'] = rdMolDescriptors.CalcTPSA(mol) # Add TPSA
    return descriptors

# ... (rest of the code) ...
```

**Example 5:  Using a Different Machine Learning Model (Python)**

*   **Python:**  Replace the Random Forest Regressor with a Support Vector Regressor (SVR).

```python
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

# ... (previous code) ...

    # Data Analysis and Modeling
    X = df[['MW', 'LogP']]
    y = np.log10(df['standard_value'])

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Feature scaling is often important for SVR
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Create an SVR model
    model = make_pipeline(SVR(kernel='rbf'))  # Use a pipeline for scaling and SVR
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    # Model Evaluation
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    print(f"Mean Squared Error: {mse}")
    print(f"R-squared: {r2}")

    # ... (rest of the code) ...
```

**Important Considerations:**

*   **Database Connection:**  Ensure your database connection details (IP address, username, password, database name) are correct.  It's best practice to store these credentials in environment variables rather than hardcoding them in your script.
*   **Error Handling:** Add more robust error handling to your code, especially when dealing with database connections, file I/O, and RDKit calculations.
*   **Data Scaling:**  For machine learning models like SVR and linear regression, scaling your data (e.g., using `StandardScaler`) can significantly improve performance.
*   **Feature Selection:**  Experiment with different sets of molecular descriptors to find the features that best predict activity.
*   **Model Tuning:**  Tune the hyperparameters of your machine learning models (e.g., `n_estimators` for Random Forest, `C` and `gamma` for SVR) using techniques like cross-validation to optimize performance.
*   **Virtual Environments:**  Use virtual environments (e.g., `venv` or `conda`) to manage your Python dependencies and ensure reproducibility.
*   **Logging:** Implement logging to track the progress of your script and debug any issues.
*   **File Paths:**  Double-check all file paths to ensure they are correct and that you have the necessary permissions to read and write files.

By following these guidelines and adapting the examples to your specific needs, you can build a powerful and effective data analysis pipeline for your ChEMBL 35 project. Let me know if you have any further questions or need more specific assistance!
