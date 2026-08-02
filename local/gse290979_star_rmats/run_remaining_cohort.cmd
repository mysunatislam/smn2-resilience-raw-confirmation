@echo off
setlocal
set "LOG=E:\smn2_gse290979_local9\star_rmats\logs\remaining_eight_onepass.scheduled.log"

echo [%date% %time%] Starting remaining STAR cohort>>"%LOG%"
wsl.exe -d Ubuntu -- bash -lc "cd /mnt/c/Users/Kotha/OneDrive/Documents/smn2/publication_repo && bash local/gse290979_star_rmats/run_local9_star_rmats.sh align" >>"%LOG%" 2>&1
set "EXIT_CODE=%ERRORLEVEL%"
echo [%date% %time%] Cohort launcher exit code %EXIT_CODE%>>"%LOG%"

exit /b %EXIT_CODE%
