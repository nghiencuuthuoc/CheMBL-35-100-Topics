#!/bin/bash

# Step 1: Prompt for SQL filename
read -p "Enter the SQL filename (without .sql): " sqlname
SQL_FILE="${sqlname}.sql"

# Step 2: Prompt for output CSV filename
read -p "Enter the desired CSV output filename (without .csv): " csvname
OUTPUT_FILE="${csvname}.csv"
TEMP_FILE="temp_result.txt"

# Step 3: Let user edit SQL content
echo "📝 Opening nano to create/edit: $SQL_FILE"
# echo "# Paste your SQL query below. Save with Ctrl+O, Enter, then exit with Ctrl+X." > "$SQL_FILE"
nano "$SQL_FILE"

# Step 4: Check if SQL file has content
if ! grep -q "[a-zA-Z0-9]" "$SQL_FILE"; then
    echo "❌ No valid SQL content found in $SQL_FILE. Aborting."
    exit 1
fi

# Step 5: Run SQL file and output to temp_result.txt
echo "▶️ Running SQL file: $SQL_FILE"
psql -U rd -d chembl_35 -h localhost -f "$SQL_FILE" > "$TEMP_FILE"

# Step 6: Convert temp_result.txt to CSV
echo "🔄 Converting to CSV: $OUTPUT_FILE"
awk '
BEGIN {FS="|"; OFS=","}
/^[-+]/ {next}
/^[ ]*$/ {next}
{
  for (i = 1; i <= NF; i++) gsub(/^ +| +$/, "", $i)
  print $0
}
' "$TEMP_FILE" > "$OUTPUT_FILE"

# Step 7: Done
echo "✅ CSV saved to $OUTPUT_FILE"
