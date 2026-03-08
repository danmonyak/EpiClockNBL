#!/bin/bash

echo -n "Permission to remove current installation of TCGAbiolinks (for for data retrieval), if necessary, and reinstall from fork ([y]/n): "
read permissionInput

permissionInput="${permissionInput:0:1}"
permissionCode=$(echo "$permissionInput" | tr '[:upper:]' '[:lower:]')

if [[ "$permissionCode" == "y" || "$permissionCode" == "" ]]; then
  Rscript setup.R
else
  echo "Please provide permission... otherwise, you must manually run:"
  echo "Rscript -e \"remove.packages('TCGAbiolinks')\""
  echo "Rscript setup.R"
fi
