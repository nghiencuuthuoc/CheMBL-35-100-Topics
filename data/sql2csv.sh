#!/bin/bash

# Step 1: Prompt for filename
read -p "Enter the SQL filename (without .sql): " filename

SQL_FILE="${filename}.sql"
TEMP_FILE="temp_result.txt"
OUTPUT_FILE="${filename}.csv"

# Step 2: Let user paste SQL into nano
echo "📝 Opening nano editor for: $SQL_FILE"
# echo "# Paste your SQL query below. Save with Ctrl+O, Enter, then exit with Ctrl+X." > "$SQL_FILE"
nano "$SQL_FILE"

# Step 3: Check if file has content
if ! grep -q "[a-zA-Z0-9]" "$SQL_FILE"; then
    echo "❌ No valid SQL content detected in $SQL_FILE. Aborting."
    exit 1
fi

# Step 4: Run SQL and redirect output to temp_result.txt
echo "▶️ Running SQL query..."
psql -U rd -d chembl_35 -h localhost -f "$SQL_FILE" > "$TEMP_FILE"

# Step 5: Convert to CSV using awk
echo "🔄 Converting result to CSV: $OUTPUT_FILE"
awk '
BEGIN {FS="|"; OFS=","}
/^[-+]/ {next}
/^[ ]*$/ {next}
{
  for (i = 1; i <= NF; i++) gsub(/^ +| +$/, "", $i)
  print $0
}
' "$TEMP_FILE" > "$OUTPUT_FILE"

# Step 6: Done
echo "✅ CSV saved to $OUTPUT_FILE"
