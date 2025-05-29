@echo off
set INPUT_FILE=temp_result.txt
set OUTPUT_FILE=output.csv

awk "BEGIN {FS=\"|\"; OFS=\",\"} \
/^[-\+]/ {next} \
/^[ ]*$/ {next} \
{for (i = 1; i <= NF; i++) gsub(/^ +| +$/, \"\", $i); print $0}" %INPUT_FILE% > %OUTPUT_FILE%

echo ✅ CSV saved to %OUTPUT_FILE%
