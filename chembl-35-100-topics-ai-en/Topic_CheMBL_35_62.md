# PharmApp Suite
## 🧠 AI for Drug Discovery and Development 🧪
| Copyright 2025 | RnD Pharma Plus | www.nghiencuuthuoc.com | Zalo: +84888999311 |

🧩 Topic: Deep learning models for activity prediction
---
Okay, I understand. You're working with ChEMBL 35 data, using PostgreSQL, RDKit, and Jupyter Notebooks in an AIMLops-structured project. You want to analyze the data, create some models, and address potential errors, limiting the data to 100 rows for efficiency. Let's break this down.

**1. Analysis Model**

The analysis will focus on exploring the ChEMBL data related to bioactivity and molecular properties. We'll perform the following steps:

*   **Data Retrieval:** Query the ChEMBL database to extract bioactivity data (e.g., IC50 values) and associated compound information (e.g., SMILES strings). We'll focus on a specific target or a set of targets to narrow down the scope.
*   **Data Cleaning and Preprocessing:** Clean the data, handling missing values, converting units, and standardizing activity values. We'll address the `numeric ~ unknown` error in the SQL query.
*   **Molecular Representation:** Use RDKit to generate molecular fingerprints or descriptors from the SMILES strings. These fingerprints represent the molecules numerically and are used as input for machine learning models.
*   **Model Building:** Build a machine learning model to predict bioactivity based on the molecular fingerprints. This could be a regression model (e.g., Linear Regression, Random Forest Regression) or a classification model (e.g., Logistic Regression, Random Forest Classifier) depending on how the bioactivity is represented (continuous or categorical).
*   **Model Evaluation:** Evaluate the performance of the model using appropriate metrics (e.g., Mean Squared Error (MSE) for regression, Accuracy, Precision, Recall, F1-score for classification).
*   **Visualization and Interpretation:** Visualize the data and the model results to gain insights into the relationship between molecular properties and bioactivity.

**2. Code (SQL & Python)**

**SQL (Topic_CheMBL_35_62.sql):**

```sql
-- data/Topic_CheMBL_35_62.csv

SELECT
    act.activity_id,
    act.standard_inchi_key,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.assay_id,
    act.target_id,
    cp.canonical_smiles
FROM
    activities act
JOIN
    compound_structures cp ON act.molregno = cp.molregno
WHERE
    act.standard_type = 'IC50'  -- Filter for IC50 values
    AND act.standard_relation = '='  -- Filter for direct measurements
    AND act.standard_value IS NOT NULL
    AND cp.canonical_smiles IS NOT NULL
    AND act.standard_value::text ~ '^[0-9\.]+$'  -- Corrected line: Cast to text and use regex
LIMIT 100;  -- Limit to 100 rows

-- Explanation:
-- 1. `act.standard_value::text ~ '^[0-9\.]+$'`:  The original error occurs because you're trying to use a regex operator (`~`) directly on a `numeric` type.  The corrected code casts `act.standard_value` to `text` before applying the regex.  This ensures that the regular expression matching can be performed.  The regex `^[0-9\.]+$` checks if the value consists only of digits and periods.
-- 2. `LIMIT 100`: Limits the result set to the first 100 rows.
-- 3. Add filter Standard_inchi_key, cp.canonical_smiles is not null
```

**Steps to run the SQL and save the results:**

1.  **Connect to your PostgreSQL database:**  Open pgAdmin and connect to the database using the IP address, username, password, and database name you provided.
2.  **Create the SQL file:** Copy the SQL code above and save it as `Topic_CheMBL_35_62.sql` in your project directory (likely under a `data` subdirectory).
3.  **Execute the SQL query:** Open the SQL file in pgAdmin and execute it.
4.  **Export the results to CSV:** After the query runs successfully, export the results as a CSV file named `Topic_CheMBL_35_62.csv` and save it in the `data` directory of your project.  Make sure the header row is included.

**Python (Topic_CheMBL_35_62_1_Data_Preprocessing.ipynb):**

```python
# Topic_CheMBL_35_62_1_Data_Preprocessing.ipynb

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Define the base path for your project
base_path = "."  # Assuming the notebook is run from the project root

# Path to the CSV file
csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_62.csv")

# Load the data
try:
    df = pd.read_csv(csv_file_path)
    print("Data loaded successfully.")
except FileNotFoundError:
    print(f"Error: File not found at {csv_file_path}")
    exit()

# Data Cleaning and Preprocessing
print("\nData Cleaning and Preprocessing...")

# Handle missing values (replace with NaN and then drop rows with NaNs)
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.dropna(subset=['canonical_smiles', 'standard_value'], inplace=True)  # Drop rows with missing SMILES or IC50 values

# Convert IC50 values to pIC50 (negative log)
def ic50_to_pic50(ic50_nM):
    """Converts IC50 (nM) to pIC50."""
    pIC50 = -np.log10(ic50_nM / 1e9)  # Convert nM to Molar
    return pIC50

# Only process rows where standard_units == "nM"
df_nM = df[df['standard_units'] == 'nM'].copy()

if not df_nM.empty:
    df_nM['pIC50'] = df_nM['standard_value'].apply(ic50_to_pic50)
    print("pIC50 values calculated for nM data.")

    # Add RDKit Mol objects
    df_nM['ROMol'] = df_nM['canonical_smiles'].apply(lambda x: Chem.MolFromSmiles(x))
    df_nM = df_nM[df_nM['ROMol'].notna()]  # Remove rows where Mol object creation failed
    print("RDKit Mol objects created.")

    # Display the first few rows of the processed DataFrame
    print("\nFirst 5 rows of processed data:")
    print(df_nM.head())
else:
    print("No data with standard_units 'nM' found.")

# Save the processed DataFrame (optional)
processed_csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_62_processed.csv")
if not df_nM.empty:
    df_nM.to_csv(processed_csv_file_path, index=False)
    print(f"\nProcessed data saved to {processed_csv_file_path}")
else:
    print("\nNo processed data to save.")
```

**Python (Topic_CheMBL_35_62_2_Model_Building.ipynb):**

```python
# Topic_CheMBL_35_62_2_Model_Building.ipynb

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import numpy as np

# Define the base path for your project
base_path = "."  # Assuming the notebook is run from the project root

# Path to the processed CSV file
processed_csv_file_path = os.path.join(base_path, "data", "Topic_CheMBL_35_62_processed.csv")

# Load the processed data
try:
    df = pd.read_csv(processed_csv_file_path)
    print("Processed data loaded successfully.")
except FileNotFoundError:
    print(f"Error: File not found at {processed_csv_file_path}")
    exit()

# Molecular Fingerprint Generation
print("\nGenerating Molecular Fingerprints...")

def generate_fingerprint(mol):
    """Generates Morgan fingerprints (ECFP4)."""
    try:
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) #Returns SparseBitVect
    except Exception as e:
        print(f"Error generating fingerprint: {e}")
        return None

df['fingerprint'] = df['ROMol'].apply(generate_fingerprint)
df.dropna(subset=['fingerprint'], inplace=True) # Drop rows where fingerprint generation failed
print("Molecular fingerprints generated.")

# Data Preparation for Modeling
print("\nPreparing data for modeling...")

# Convert fingerprints to numpy arrays
X = np.array([np.array(x) for x in df['fingerprint']]) #Convert SparseBitVect to numpy array
y = df['pIC50']

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize the data (important for linear models)
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Model Building
print("\nBuilding and training the Linear Regression model...")

model = LinearRegression()
model.fit(X_train, y_train)

# Model Evaluation
print("\nEvaluating the model...")

y_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse:.2f}")
print(f"R-squared: {r2:.2f}")

# Interpretation (Example)
print("\nModel Interpretation (Example):")
print("Coefficients for the first 5 features:", model.coef_[:5])

# Visualization (Example - requires matplotlib)
# import matplotlib.pyplot as plt
# plt.scatter(y_test, y_pred)
# plt.xlabel("Actual pIC50")
# plt.ylabel("Predicted pIC50")
# plt.title("Actual vs. Predicted pIC50")
# plt.show()
```

**Explanation of Python Code:**

*   **Data Loading:** Loads the CSV file generated from the SQL query using `pandas`.  Handles potential `FileNotFoundError`.
*   **Data Cleaning:** Handles missing values by dropping rows with missing SMILES or pIC50 values.
*   **pIC50 Conversion:**  Converts IC50 values (in nM) to pIC50 values.  This is a common transformation in drug discovery, as pIC50 values are more normally distributed and easier to work with.
*   **RDKit Mol Object Creation:** Creates RDKit `Mol` objects from the SMILES strings.  These `Mol` objects are required for generating molecular fingerprints.
*   **Molecular Fingerprint Generation:** Uses RDKit to generate Morgan fingerprints (ECFP4).  These fingerprints are numerical representations of the molecules and are used as input features for the machine learning model.
*   **Data Splitting:** Splits the data into training and testing sets.
*   **Model Building:**  Creates and trains a Linear Regression model.
*   **Model Evaluation:** Evaluates the model using Mean Squared Error (MSE) and R-squared.
*   **Model Interpretation:**  Prints the coefficients of the first 5 features of the linear regression model. This provides some insight into which features are most important for predicting pIC50.
*   **Visualization (Commented Out):** Includes example code for visualizing the results using a scatter plot of actual vs. predicted pIC50 values.  This requires `matplotlib`.
*   **Error Handling:** The code includes `try...except` blocks to handle potential errors, such as `FileNotFoundError` when loading the CSV file or errors during fingerprint generation.  If an error occurs, the code prints an informative message and continues to the next row.

**Addressing the Errors:**

*   **SQL Error (`numeric ~ unknown`):** The fix is in the SQL code above.  Cast the `standard_value` to `text` before using the regex operator: `act.standard_value::text ~ '^[0-9\.]+$'`
*   **Scikit-learn Error (`squared=False`):**  The `squared=False` parameter in `mean_squared_error` was introduced in a later version of scikit-learn.  If you're using an older version, remove the `squared=False` parameter.  Alternatively, upgrade your scikit-learn version to the latest one. Since MSE is a popular metric, it is not useful.

**Important Considerations:**

*   **Dependencies:** Make sure you have the necessary Python packages installed: `pandas`, `rdkit`, `scikit-learn`, `numpy`.  You can install them using `pip`:
    ```bash
    pip install pandas rdkit scikit-learn numpy
    ```
*   **RDKit Installation:**  RDKit can be a bit tricky to install.  It's recommended to use conda if possible:
    ```bash
    conda install -c conda-forge rdkit
    ```
*   **File Paths:** Double-check the file paths in the Python code to make sure they are correct for your project structure.
*   **Data Exploration:** Before building the model, it's crucial to explore the data.  Look at the distribution of pIC50 values, the range of fingerprint values, and any potential biases in the data.  This will help you choose the right model and interpret the results.
*   **Model Selection:** Linear Regression is a simple model, but it may not be the best choice for all datasets.  Consider trying other models, such as Random Forest Regression or Support Vector Regression.
*   **Hyperparameter Tuning:** Tune the hyperparameters of the machine learning model to improve its performance.  You can use techniques like cross-validation and grid search to find the optimal hyperparameters.
*   **Data Filtering:** Filter by inchi_key, canonical_smiles to reduce error handling.

**3. Five Examples (Scenario-Based)**

Here are five examples of how you could use this analysis in a pharmaceutical R&D context:

1.  **Target Identification:**  Analyze the bioactivity data for a specific target protein (e.g., a kinase involved in cancer).  The model can help identify compounds that are likely to be active against that target.
2.  **Lead Optimization:**  You have a lead compound and want to improve its potency.  Use the model to predict the pIC50 of various analogs of the lead compound, based on their molecular fingerprints.  This can help you prioritize which analogs to synthesize and test.
3.  **Off-Target Prediction:**  You want to assess the potential for off-target effects of a drug candidate.  Build models for multiple targets and use them to predict the activity of the drug candidate against each target.  This can help identify potential safety concerns.
4.  **Virtual Screening:**  You have a large database of compounds and want to identify potential hits for a target of interest.  Use the model to predict the pIC50 of each compound in the database.  This can help you prioritize which compounds to screen experimentally.
5.  **Structure-Activity Relationship (SAR) Analysis:**  Use the model to understand the relationship between the structure of a molecule and its activity.  Analyze the coefficients of the linear regression model or the feature importance of a random forest model to identify which molecular features are most important for activity.

**In Summary (En résumé):**

You're taking on a well-defined project. The SQL and Python code provides a solid foundation. Remember to adjust file paths, install dependencies, address errors, and explore your data thoroughly. By following these steps, you can build a predictive model and gain valuable insights from the ChEMBL 35 data.

**Analyse (Analyse):**

Vous vous lancez dans un projet bien défini. Le code SQL et Python fournit une base solide. N'oubliez pas d'ajuster les chemins d'accès aux fichiers, d'installer les dépendances, de corriger les erreurs et d'explorer vos données en profondeur. En suivant ces étapes, vous pouvez créer un modèle prédictif et obtenir des informations précieuses à partir des données ChEMBL 35.
