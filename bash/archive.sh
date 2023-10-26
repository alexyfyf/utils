#!/bin/bash

## version 1.01
## 20231026

## used in Biomate vault to tar and archive files to tape
## also generate md5 and check result

# Check if a directory is passed as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 directory_name"
  exit 1
fi

DIRECTORY=$(echo $1 | sed 's/\/$//')

# Step into the directory
cd "$DIRECTORY" || exit

# Create an md5 checksum file for all files in the directory
find . -type f -exec md5sum {} + > ../"${DIRECTORY}.manifest"

# Step back out
cd ..

# Create a tarball of the directory, including the md5 checksum file
tar -czf "${DIRECTORY}.tar.gz" "${DIRECTORY}" "${DIRECTORY}.manifest"

# Extract the tarball to a temporary directory for verification
mkdir -p temp_verification_dir
tar -xzf "${DIRECTORY}.tar.gz" -C temp_verification_dir

# Verify the integrity using the md5 checksum file
cd temp_verification_dir/"$DIRECTORY" || exit
md5sum -c "../${DIRECTORY}.manifest" > ../../${DIRECTORY}_checkresult

# Cleanup: remove the temporary directory and checksum file
cd ../..
rm -r temp_verification_dir

echo "The script has completed. Please check the output above for any discrepancies in the md5 checksum verification."

echo Total file number
cat ${DIRECTORY}_checkresult | wc -l

echo Passed file number
grep ': OK' ${DIRECTORY}_checkresult | wc -l 