#!/bin/bash

INPUT_FILE="temp_result.txt"
read -p "Enter the output CSV file name (e.g. result.csv): " OUTPUT_FILE

awk '
BEGIN {FS="|"; OFS=","}
/^[-+]/ {next}
/^[ ]*$/ {next}
{
  for (i = 1; i <= NF; i++) gsub(/^ +| +$/, "", $i)
  print $0
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "✅ CSV saved to $OUTPUT_FILE"
