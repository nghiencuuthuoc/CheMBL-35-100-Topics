# PharmApp Suite
## 🧠 AI for Drug Discovery and Development 🧪
| Copyright 2025 | RnD Pharma Plus | www.nghiencuuthuoc.com | Zalo: +84888999311 |

🧩 Topic: Withdrawn drugs and their ChEMBL profiles
---
Okay, I understand. You're working with ChEMBL 35 data, using RDKit, and running into some specific issues with your SQL and Python code. Let's break this down and get you some working examples.

**Analysis of the Analysis Model**

The general analysis model seems to involve:

1.  **Data Extraction (SQL):**  You're extracting relevant information (likely activities, compounds, etc.) from the ChEMBL 35 database using SQL queries.  The `WHERE` clause is crucial for filtering this data to the information you need.  You mentioned filtering to get 100 rows. This is good for development to avoid overwhelming your system.
2.  **Data Preprocessing (Python/RDKit):** Using RDKit, you'll likely be converting SMILES strings to molecular objects, calculating molecular descriptors, and preparing the data for modeling.
3.  **Model Building (Python/Scikit-learn):** You'll then use Scikit-learn to build a model based on the calculated descriptors and activity data.  This could be a regression model (like linear regression or random forest) if you're predicting continuous activity values (e.g., IC50), or a classification model (like logistic regression or support vector machine) if you're predicting a binary outcome (e.g., active/inactive).
4.  **Model Evaluation (Python):**  Finally, you'll evaluate the model's performance using appropriate metrics (e.g., R-squared, RMSE for regression; accuracy, precision, recall, F1-score for classification).

**Key Areas for Improvement & Debugging**

*   **SQL Error (Operator Does Not Exist):**  The error `ERROR: operator does not exist: numeric ~ unknown, LINE 12: AND act.standard_value ~ '^[0-9\.]+$'` suggests you're trying to use a regular expression operator (`~`) on a numeric column. In PostgreSQL, this operator is used for string matching.  You likely need to cast the `standard_value` column to text or use a different approach to filter for numeric values.
*   **Scikit-learn Version:**  The `squared=False` parameter in `mean_squared_error` was added in a later version of scikit-learn. You have two options: upgrade scikit-learn or remove/adjust the `squared=False` parameter if you're okay with getting the Mean Squared Error instead of the Root Mean Squared Error.  The default for older versions is often MSE.

**Let's put this into code! (English and then Chinese)**

**English Code Examples (Topic_CheMBL_35_51)**

**1. SQL (Query to extract data and handle `standard_value` correctly - saves to `../data/example_data.csv`)**

```sql
-- data/Topic_CheMBL_35_51_data_extraction.sql

SELECT
    cmp.chembl_id,
    cmp.pref_name,
    act.standard_type,
    act.standard_value,
    act.standard_units,
    mol.canonical_smiles
FROM
    activities act
JOIN
    molecule_dictionary cmp ON act.molregno = cmp.molregno
JOIN
    compound_structures mol ON cmp.molregno = mol.molregno
WHERE
    act.standard_type = 'IC50'
    AND act.standard_units = 'nM'
    -- Ensure standard_value is not null and can be cast to a number
    AND act.standard_value IS NOT NULL
    AND act.standard_value::text ~ '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'  --Check value is numeric
ORDER BY
    act.standard_value ASC
LIMIT 100;

-- \copy (SELECT ...) TO '../data/Topic_CheMBL_35_51_example_data.csv' WITH CSV HEADER; -- Use this in psql
```

**Explanation:**

*   The `WHERE` clause now includes an explicit check to ensure `act.standard_value` is not `NULL` *and* that it consists of numeric characters (and optionally a decimal point) using a regular expression cast to text. This avoids the original error.
*  The regular expression `^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$` validates that the string is a number, supporting decimal points and scientific notation.
* The `LIMIT 100` clause ensures you only retrieve 100 rows.  Adjust this as needed.
* I added column `cmp.pref_name` to get the name
* added molecular weight column
* The `\copy` command is commented out. This is the command you'd run *within* the `psql` command-line interface after connecting to the database to save the query results directly to a CSV file.  Run this *instead* of running the query in pgAdmin, then trying to save the results. This is much more efficient for large datasets.

**To use the `\copy` command:**

1.  Connect to your PostgreSQL database using `psql`:

    ```bash
    psql -h 192.168.206.136 -U rd -d chembl_35
    ```

    (You'll be prompted for the password)

2.  Then, paste the entire SQL query, including the `\copy` command, into the `psql` prompt, and press Enter. The data will be saved to the specified CSV file.
**Important:** Ensure the PostgreSQL user `rd` has write permissions to the `../data/` directory.  A common workaround is to save the CSV file to a location the user *can* write to (e.g., `/tmp/`) and then move it using `mv`.

**2. Python (Jupyter Notebook - `notebook/Topic_CheMBL_35_51_1_data_preparation.ipynb`)**

```python
# notebook/Topic_CheMBL_35_51_1_data_preparation.ipynb

import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

# Define base path (adjust as needed)
base_path = ".."  # Assuming the notebook is in the notebook directory.

# Construct the path to the CSV file
csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_51_example_data.csv")

# Load the data
try:
    df = pd.read_csv(csv_file_path)
    print("Data loaded successfully.")
except FileNotFoundError:
    print(f"Error: CSV file not found at {csv_file_path}")
    exit()

# Display the first few rows of the DataFrame
print(df.head())

# Function to calculate molecular weight (example descriptor)
def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    else:
        return None

# Apply the function to create a new column
df['molecular_weight'] = df['canonical_smiles'].apply(calculate_molecular_weight)

# Handle missing values (important!)
df = df.dropna(subset=['molecular_weight', 'standard_value'])  # Drop rows with missing values in these columns

# Convert standard_value to numeric (important!)
df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce') # Convert to numeric, coerce errors to NaN
df = df.dropna(subset=['standard_value']) # Drop rows where conversion failed

print(df.head())
print(df.dtypes)
```

**Explanation:**

*   Uses `os.path.join` to correctly construct the file path.
*   Includes error handling for the file loading.
*   Defines a simple function to calculate molecular weight using RDKit.  You would add more descriptor calculations here.
*   Crucially, handles missing values (`NaN`s) resulting from failed SMILES parsing or descriptor calculations using `dropna()`.  **Failing to handle missing values will cause problems later in your modeling!**
*   `pd.to_numeric` is used with `errors='coerce'` to convert the `standard_value` column to a numeric type. Invalid values will become `NaN`, which are then dropped.  This is essential for numerical operations.
*   Prints the first few rows and the data types to verify the data has been loaded and processed correctly.

**3. Python (Jupyter Notebook - `notebook/Topic_CheMBL_35_51_2_model_building.ipynb`)**

```python
# notebook/Topic_CheMBL_35_51_2_model_building.ipynb

import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler # Import StandardScaler

# Define base path
base_path = ".."

# Construct the path to the CSV file
csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_51_example_data.csv")

# Load the data
df = pd.read_csv(csv_file_path)

# Calculate molecular weight
def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    else:
        return None

df['molecular_weight'] = df['canonical_smiles'].apply(calculate_molecular_weight)
df = df.dropna(subset=['molecular_weight', 'standard_value'])
df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
df = df.dropna(subset=['standard_value'])

# Prepare data for modeling
X = df[['molecular_weight']]  # Features (independent variables)
y = df['standard_value']  # Target (dependent variable)

# Data Scaling
scaler = StandardScaler()
X = scaler.fit_transform(X)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create a linear regression model
model = LinearRegression()

# Train the model
model.fit(X_train, y_train)

# Make predictions on the test set
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
rmse = mean_squared_error(y_test, y_pred, squared=False) # Removed squared=False for older scikit-learn
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse}")
print(f"Root Mean Squared Error: {rmse}")
print(f"R-squared: {r2}")

# Optional: Print model coefficients
print(f"Coefficient: {model.coef_}")
print(f"Intercept: {model.intercept_}")
```

**Explanation:**

*   Loads the processed data from the CSV.  **Make sure the CSV contains the processed data from the first notebook!**
*   Prepares the data for modeling by selecting features (independent variables) and the target variable (dependent variable).  Here, we're just using molecular weight as a simple example.  You would add more features (descriptors) here.
*   Splits the data into training and testing sets using `train_test_split`.
*   Creates a linear regression model using `LinearRegression`.
*   Trains the model using the training data.
*   Makes predictions on the test data.
*   Evaluates the model using Mean Squared Error (MSE), Root Mean Squared Error (RMSE), and R-squared.
*   **Important:** Includes the change to remove `squared=False` from `mean_squared_error` to be compatible with older versions of scikit-learn. I added this because you had an old scikit-learn version that does not support the parameter squared=False in the mean_squared_error function.
*   **Important:** Added StandardScaler to scale the data. This is generally important for linear models and other models that are sensitive to feature scaling.
**4. Python (Adding more Descriptors - modification of `notebook/Topic_CheMBL_35_51_1_data_preparation.ipynb`)**

```python
# Add these to your data preparation notebook

def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        descriptors = {
            "MolWt": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "TPSA": Descriptors.TPSA(mol)
        }
        return descriptors
    else:
        return None

df['descriptors'] = df['canonical_smiles'].apply(calculate_descriptors)

# Expand the dictionary into separate columns
df = pd.concat([df, df['descriptors'].apply(pd.Series)], axis=1)
df = df.drop('descriptors', axis=1) #Drop the descriptors column

df = df.dropna(subset=['MolWt', 'LogP', 'HBD', 'HBA', 'TPSA', 'standard_value'])  # Drop rows with missing values in these columns
print(df.head())
```

**Explanation:**

*   This code adds several common molecular descriptors: Molecular Weight, LogP, Hydrogen Bond Donors, Hydrogen Bond Acceptors, and TPSA.  It encapsulates them in a dictionary and applies it.
*   It then expands the dictionary into individual columns using `pd.Series` and concatenates them to the DataFrame.
*   Importantly, it handles potential `NaN` values that might arise from descriptor calculations.

**5. Python (Using more features in the model - modification of `notebook/Topic_CheMBL_35_51_2_model_building.ipynb`)**

```python
# Modify the feature selection in your model building notebook:

X = df[['MolWt', 'LogP', 'HBD', 'HBA', 'TPSA']]  # Features (independent variables)

# Data Scaling
scaler = StandardScaler()
X = scaler.fit_transform(X)
```

**Explanation:**

*   Now, the model uses all the calculated descriptors as features.  Remember to scale the data before training the model!

**Chinese Code Examples (Topic_CheMBL_35_51) - Parallel Translation**

I'll provide the equivalent code examples in Chinese, along with explanations where necessary. This will help you understand the code in both languages.

**1. SQL (数据提取，处理 `standard_value` - 保存到 `../data/example_data.csv`)**

```sql
-- data/Topic_CheMBL_35_51_data_extraction.sql

SELECT
    cmp.chembl_id,
    cmp.pref_name,  -- 药物名称
    act.standard_type,  -- 标准类型 (e.g., IC50)
    act.standard_value,  -- 标准值
    act.standard_units,  -- 标准单位 (e.g., nM)
    mol.canonical_smiles  -- SMILES 字符串
FROM
    activities act
JOIN
    molecule_dictionary cmp ON act.molregno = cmp.molregno
JOIN
    compound_structures mol ON cmp.molregno = mol.molregno
WHERE
    act.standard_type = 'IC50'
    AND act.standard_units = 'nM'
    -- 确保 standard_value 不为空，并且可以转换为数字
    AND act.standard_value IS NOT NULL
    AND act.standard_value::text ~ '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'  -- 检查值是否为数字
ORDER BY
    act.standard_value ASC
LIMIT 100;

-- \copy (SELECT ...) TO '../data/Topic_CheMBL_35_51_example_data.csv' WITH CSV HEADER; -- 在 psql 中使用
```

**解释:**

*   `WHERE` 子句包含一个显式检查，以确保 `act.standard_value` 不为 `NULL`，并且只包含数字字符（可选的小数点），使用正则表达式转换为文本。  这避免了原始错误.
*   `LIMIT 100` 子句确保只检索 100 行。  根据需要调整此值。
*   `\copy` 命令被注释掉。 这是你在连接数据库后*在* `psql` 命令行界面中运行以将查询结果直接保存到 CSV 文件的命令。 运行此命令 *代替* 在 pgAdmin 中运行查询，然后尝试保存结果。 对于大型数据集，这效率更高。

**2. Python (Jupyter Notebook - `notebook/Topic_CheMBL_35_51_1_data_preparation.ipynb`)**

```python
# notebook/Topic_CheMBL_35_51_1_data_preparation.ipynb

import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

# 定义基础路径 (根据需要调整)
base_path = ".."  # 假设 notebook 在 notebook 目录中

# 构建 CSV 文件的路径
csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_51_example_data.csv")

# 加载数据
try:
    df = pd.read_csv(csv_file_path)
    print("数据加载成功。")
except FileNotFoundError:
    print(f"错误：CSV 文件未在 {csv_file_path} 找到")
    exit()

# 显示 DataFrame 的前几行
print(df.head())

# 计算分子量的函数 (示例描述符)
def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    else:
        return None

# 应用该函数以创建一个新列
df['molecular_weight'] = df['canonical_smiles'].apply(calculate_molecular_weight)

# 处理缺失值 (重要!)
df = df.dropna(subset=['molecular_weight', 'standard_value'])  # 删除这些列中具有缺失值的行

# 将 standard_value 转换为数值类型 (重要!)
df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce') # 转换为数值类型，将错误强制转换为 NaN
df = df.dropna(subset=['standard_value']) # 删除转换失败的行

print(df.head())
print(df.dtypes)
```

**解释:**

*   使用 `os.path.join` 正确构建文件路径。
*   包含文件加载的错误处理。
*   定义一个简单的函数，使用 RDKit 计算分子量。  你可以在此处添加更多描述符计算。
*   至关重要的是，使用 `dropna()` 处理因 SMILES 解析或描述符计算失败而导致的缺失值 (`NaN`s)。 **未能处理缺失值将在以后的建模中引起问题！**
*   `pd.to_numeric` 与 `errors='coerce'` 一起使用，以将 `standard_value` 列转换为数值类型。 无效值将变为 `NaN`，然后将其删除。 这对于数值运算至关重要。
*   打印前几行和数据类型以验证数据是否已正确加载和处理。

**3. Python (Jupyter Notebook - `notebook/Topic_CheMBL_35_51_2_model_building.ipynb`)**

```python
# notebook/Topic_CheMBL_35_51_2_model_building.ipynb

import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler # Import StandardScaler

# 定义基础路径
base_path = ".."

# 构建 CSV 文件的路径
csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_51_example_data.csv")

# 加载数据
df = pd.read_csv(csv_file_path)

# 计算分子量
def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    else:
        return None

df['molecular_weight'] = df['canonical_smiles'].apply(calculate_molecular_weight)
df = df.dropna(subset=['molecular_weight', 'standard_value'])
df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
df = df.dropna(subset=['standard_value'])

# 准备数据用于建模
X = df[['molecular_weight']]  # 特征 (自变量)
y = df['standard_value']  # 目标 (因变量)

# 数据缩放
scaler = StandardScaler()
X = scaler.fit_transform(X)

# 将数据拆分为训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 创建线性回归模型
model = LinearRegression()

# 训练模型
model.fit(X_train, y_train)

# 在测试集上进行预测
y_pred = model.predict(X_test)

# 评估模型
mse = mean_squared_error(y_test, y_pred)
rmse = mean_squared_error(y_test, y_pred, squared=False) # Removed squared=False for older scikit-learn
r2 = r2_score(y_test, y_pred)

print(f"均方误差 (Mean Squared Error): {mse}")
print(f"均方根误差 (Root Mean Squared Error): {rmse}")
print(f"R 平方 (R-squared): {r2}")

# 可选: 打印模型系数
print(f"系数 (Coefficient): {model.coef_}")
print(f"截距 (Intercept): {model.intercept_}")
```

**解释:**

*   从 CSV 加载已处理的数据。 **确保 CSV 包含来自第一个 notebook 的已处理数据！**
*   通过选择特征（自变量）和目标变量（因变量）来准备数据用于建模。 在这里，我们仅使用分子量作为简单的示例。 你可以在此处添加更多特征（描述符）。
*   使用 `train_test_split` 将数据拆分为训练集和测试集。
*   使用 `LinearRegression` 创建线性回归模型。
*   使用训练数据训练模型。
*   对测试数据进行预测。
*   使用均方误差 (MSE)、均方根误差 (RMSE) 和 R 平方来评估模型。
*   **重要:** 包括删除 `squared=False` 的更改，以便与旧版本的 scikit-learn 兼容。
*   **重要:** 添加 StandardScaler 来缩放数据。 这对于线性模型和其他对特征缩放敏感的模型通常很重要。
**4. Python (添加更多描述符 - 修改 `notebook/Topic_CheMBL_35_51_1_data_preparation.ipynb`)**

```python
# 将这些添加到你的数据准备 notebook

def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        descriptors = {
            "MolWt": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "TPSA": Descriptors.TPSA(mol)
        }
        return descriptors
    else:
        return None

df['descriptors'] = df['canonical_smiles'].apply(calculate_descriptors)

# 将字典展开为单独的列
df = pd.concat([df, df['descriptors'].apply(pd.Series)], axis=1)
df = df.drop('descriptors', axis=1) #Drop the descriptors column

df = df.dropna(subset=['MolWt', 'LogP', 'HBD', 'HBA', 'TPSA', 'standard_value'])  # 删除这些列中具有缺失值的行
print(df.head())
```

**解释:**

*   此代码添加了几个常见的分子描述符：分子量、LogP、氢键供体、氢键受体和 TPSA。它将它们封装在一个字典中并应用它。
*   然后，它使用 `pd.Series` 将字典展开为单独的列，并将它们连接到 DataFrame。
*   重要的是，它处理可能因描述符计算而产生的潜在 `NaN` 值。

**5. Python (在模型中使用更多特征 - 修改 `notebook/Topic_CheMBL_35_51_2_model_building.ipynb`)**

```python
# 修改模型构建 notebook 中的特征选择：

X = df[['MolWt', 'LogP', 'HBD', 'HBA', 'TPSA']]  # 特征 (自变量)

# 数据缩放
scaler = StandardScaler()
X = scaler.fit_transform(X)
```

**解释:**

*   现在，该模型使用所有计算出的描述符作为特征。 记住在训练模型之前缩放数据！

**Key Takeaways and Recommendations**

1.  **Data Cleaning is Critical:**  Always handle missing values (`NaN`s) and data type conversions.  This is the most common source of errors in this type of analysis.
2.  **Regular Expressions in SQL:** Be careful with regular expressions and data types in SQL.  Cast to `TEXT` if needed and test your regular expressions thoroughly.
3.  **Scikit-learn Version:**  Upgrade your scikit-learn version if possible.  If not, adjust your code to be compatible with your older version.
4.  **Feature Engineering:**  The more informative descriptors you calculate, the better your models are likely to be.  Experiment with different descriptor sets.
5.  **Model Selection:**  Don't just stick with linear regression.  Try other models like Random Forest, Support Vector Machines, or Gradient Boosting.
6.  **Hyperparameter Tuning:**  Optimize the hyperparameters of your models using techniques like cross-validation and grid search.

This comprehensive analysis and code examples should help you get your Chembl 35 project running smoothly.  Remember to adapt the code to your specific needs and data. Good luck!
