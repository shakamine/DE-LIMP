@echo off
:: ============================================================================
::  Launch_DE-LIMP_WSL.bat — One-click launcher for DE-LIMP on Windows via WSL2
:: ============================================================================
::
::  Runs DE-LIMP natively in a WSL2 Ubuntu distro. No Docker needed.
::  On first run, installs R + all packages inside WSL (20-30 min).
::  Subsequent runs start in ~30 seconds.
::
::  Why WSL over Docker:
::    - No .NET, no 9p msize tuning, no SSH key chmod gymnastics
::    - Bioconductor packages (limpa, ComplexHeatmap, MOFA2) compile cleanly
::    - Faster startup once installed
::
::  Requirements:
::    - Windows 10 2004+ or Windows 11
::    - WSL2 with Ubuntu (installed automatically if missing)
::
:: ============================================================================

setlocal enableextensions
title DE-LIMP (WSL)

echo.
echo  ============================================
echo    DE-LIMP Proteomics (WSL2 launcher)
echo  ============================================
echo.

:: ----------------------------------------------------------------------------
:: 1. Check WSL is installed
:: ----------------------------------------------------------------------------
where wsl >nul 2>&1
if errorlevel 1 (
    echo  ERROR: wsl.exe not found. You need Windows 10 2004+ or Windows 11.
    echo.
    echo  To install WSL2, open PowerShell as Administrator and run:
    echo      wsl --install
    echo  Then restart Windows and re-run this launcher.
    echo.
    pause
    exit /b 1
)

:: ----------------------------------------------------------------------------
:: 2. Check that an Ubuntu distro is installed and reachable
::
::    `wsl --list --quiet` on Windows 11 outputs UTF-16, which findstr can't
::    read. The previous probe `wsl -d Ubuntu -e true && if errorlevel 1`
::    was unreliable: Windows batch's `if errorlevel N` evaluates as
::    "errorlevel >= N", and WSL returns negative-ish exit codes for
::    `WSL_E_DISTRO_NOT_FOUND` that the test misinterprets as success.
::    The launcher then proceeded to copy/run, only to fail downstream
::    with cryptic "There is no distribution with the supplied name."
::
::    v3.10.16 — sentinel-string probe. Run a command that prints a known
::    sentinel; if the sentinel isn't in the output, the distro isn't
::    there. This is exit-code-independent and works regardless of
::    Windows / WSL version quirks.
:: ----------------------------------------------------------------------------
wsl -d Ubuntu -e bash -c "echo __DELIMP_UBUNTU_OK__" 2>&1 | findstr /c:"__DELIMP_UBUNTU_OK__" >nul
if errorlevel 1 (
    echo  No working Ubuntu distro in WSL. Installing Ubuntu now...
    echo  ^(This opens a separate window — follow the prompts, then close it and re-run this launcher.^)
    echo.
    wsl --install -d Ubuntu
    echo.
    echo  Ubuntu install triggered. After it finishes setting up, re-run this launcher.
    pause
    exit /b 0
)

:: ----------------------------------------------------------------------------
:: 3. Copy the setup script into WSL (so it works even if this bat lives
::    on a Windows path WSL can't easily reach, like OneDrive)
:: ----------------------------------------------------------------------------
set "SETUP_SRC=%~dp0delimp_wsl_setup.sh"
if not exist "%SETUP_SRC%" (
    echo  ERROR: delimp_wsl_setup.sh not found next to this .bat file.
    echo  Expected at: %SETUP_SRC%
    echo.
    echo  Make sure you cloned the full DE-LIMP repo, not just this file.
    pause
    exit /b 1
)

echo  Copying setup script into WSL...
:: Convert Windows path to WSL path via wslpath, then cp to ~
wsl -d Ubuntu -e bash -c "cp \"$(wslpath -u '%SETUP_SRC%')\" ~/delimp_wsl_setup.sh && chmod +x ~/delimp_wsl_setup.sh"
if errorlevel 1 (
    echo  ERROR: Failed to copy setup script into WSL.
    pause
    exit /b 1
)

:: ----------------------------------------------------------------------------
:: 4. Run the app (setup is auto-triggered on first run)
:: ----------------------------------------------------------------------------
echo.
echo  ============================================
echo    Starting DE-LIMP...
echo  ============================================
echo.
echo    FIRST INSTALL: 20-30 MINUTES
echo    (R + Bioconductor packages compile from source)
echo.
echo    Subsequent runs: under 1 minute.
echo.
echo    Watch this terminal for progress. When it shows
echo        "Listening on http://0.0.0.0:3838"
echo    the browser opens automatically.
echo.
echo    The browser-open command waits for the app to
echo    actually be listening on port 3838 before opening
echo    -- it will NOT open prematurely on first install.
echo.
echo    You can also open http://localhost:3838 manually
echo    anytime once the "Listening on" line appears.
echo.
echo  ============================================
echo.

:: v3.10.17 — wait for the Shiny port to actually be listening before
:: opening the browser. The previous fixed 90-second delay was way too
:: short for a fresh install (20-30 min of R + Bioconductor compiles),
:: so the browser opened to a "can't connect" error and confused users.
::
:: Polls every 2 seconds for up to 60 minutes. Opens the browser only
:: when port 3838 is listening. Silently exits if the timeout is hit
:: (no false-positive open on a stuck install).
start "" /B powershell -NoProfile -WindowStyle Hidden -Command ^
  "$d=3600; while($d -gt 0){ if (Get-NetTCPConnection -LocalPort 3838 -State Listen -ErrorAction SilentlyContinue) { Start-Process 'http://localhost:3838'; exit 0 }; Start-Sleep -Seconds 2; $d -= 2 }"

:: Run inside WSL — this blocks until the app is stopped with Ctrl+C
wsl -d Ubuntu -e bash -c "bash ~/delimp_wsl_setup.sh"

echo.
echo  DE-LIMP has stopped.
pause
endlocal
