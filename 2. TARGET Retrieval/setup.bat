@echo off
setlocal

set "permissionInput="
set /p permissionInput=Permission to remove current installation of TCGAbiolinks (for data retrieval), if necessary, and reinstall from fork ([y]/n): 

:: Take only the first character
set "permissionCode=%permissionInput:~0,1%"

:: Accept Y or y
if /I "%permissionCode%"=="Y" goto install

:: Accept just pressing Enter
if not defined permissionInput goto install

echo Please provide permission... otherwise, you must manually run:
echo Rscript -e "remove.packages('TCGAbiolinks')"
echo Rscript setup.R
goto end

:install
Rscript setup.R

:end
endlocal