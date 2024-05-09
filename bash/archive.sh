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

RAWDIRECTORY=$(echo $1 | sed 's/\/$//')

# cp directory to current location in case no write permission
DIRECTORY=$(basename $RAWDIRECTORY)
cp -r $RAWDIRECTORY .

# Create an md5 checksum file for all files in the directory
find "$DIRECTORY" -type f -exec md5sum {} + > "${DIRECTORY}.manifest"

# Create a tarball of the directory
tar -czf "${DIRECTORY}.tar.gz" "${DIRECTORY}" 

# Extract the tarball to a temporary directory for verification
mkdir -p temp_verification_dir
tar -xzf "${DIRECTORY}.tar.gz" -C temp_verification_dir

# Verify the integrity using the md5 checksum file
cd temp_verification_dir || exit
md5sum -c "../${DIRECTORY}.manifest" > ../${DIRECTORY}.checkresult

# Cleanup: remove the temporary directory and checksum file
cd ..
rm -r temp_verification_dir $DIRECTORY

echo "The script has completed. Please check the output above for any discrepancies in the md5 checksum verification."

echo Total file number
cat ${DIRECTORY}.checkresult | wc -l

echo Passed file number
grep ': OK' ${DIRECTORY}.checkresult | wc -l 

## after generate tar.gz, manifest and checkresults, generate a md5sum for the 3 files for later transfer