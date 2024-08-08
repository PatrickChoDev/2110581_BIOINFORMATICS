#!/bin/bash

# Define the input and output files
script_filename="$(realpath $1)"

# Check if $2 is a Windows file path
if [[ "$2" == *"\\"* ]]; then
  echo "Windows file path detected."
  # Add code here to handle Windows file path
  # Strip Windows file name to WSL path
  wsl_path="${2//\\//}"
  wsl_path="${realpath $wsl_path}"

  echo "WSL file path: $wsl_path"
  input_filename="$(realpath $wsl_path)"
else
  echo "Linux file path detected."
  input_filename="$(realpath $2)"
fi

if [ -f "$output_filename" ]; then
    rm "$output_filename"
fi

python "$script_filename" < "$input_filename" >> "$output_filename"

cat "$output_filename"