#!/bin/bash

# Prompt user to enter the base filename (without extension)
read -p "Enter the SQL filename (without .sql): " filename

SQL_FILE="${filename}.sql"
TEMP_FILE="temp_result.txt"
OUTPUT_FILE="${filename}.csv"

# Check if SQL file exists
if [ ! -f "$SQL_FILE" ]; then
    echo "❌ File $SQL_FILE not found!"
    exit 1
fi

# Step 1: Run the SQL file and output to temp_result.txt
echo "▶️ Running SQL file: $SQL_FILE ..."
psql -U rd -d chembl_35 -h localhost -f "$SQL_FILE" > "$TEMP_FILE"

# Step 2: Convert the output to CSV
echo "🔄 Converting to CSV: $OUTPUT_FILE ..."
awk '
BEGIN {FS="|"; OFS=","}
/^[-+]/ {next}
/^[ ]*$/ {next}
{
  for (i = 1; i <= NF; i++) gsub(/^ +| +$/, "", $i)
  print $0
}
' "$TEMP_FILE" > "$OUTPUT_FILE"

# Done
echo "✅ CSV saved to $OUTPUT_FILE"
