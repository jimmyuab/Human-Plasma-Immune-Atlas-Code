@echo off
REM ============================================================
REM  One-click publish to GitHub
REM  Edit the two variables below, then double-click this file.
REM ============================================================
setlocal enabledelayedexpansion

set REPO_NAME=plasma-immunome-phenome-atlas
set VISIBILITY=public
set RELEASE_TAG=v1.0.0

cd /d "%~dp0"
echo(
echo === Publishing "%REPO_NAME%" (%VISIBILITY%) ===
echo(

REM --- 1. git identity check ---
git rev-parse --is-inside-work-tree >nul 2>&1
if errorlevel 1 (
    echo Initializing git repository...
    git init -b main
)

REM --- 2. stage + commit ---
git add -A
git diff --cached --quiet
if errorlevel 1 (
    git commit -m "Release v1.0.0: open genetics-anchored plasma immunome atlas"
) else (
    echo Nothing new to commit.
)

REM --- 3. create remote + push ---
where gh >nul 2>&1
if %errorlevel%==0 (
    echo Using GitHub CLI...
    gh auth status >nul 2>&1
    if errorlevel 1 (
        echo You are not logged into gh. Run:  gh auth login
        goto :end
    )
    git remote get-url origin >nul 2>&1
    if errorlevel 1 (
        gh repo create %REPO_NAME% --%VISIBILITY% --source=. --remote=origin --push
    ) else (
        git push -u origin main
    )
    echo Creating release %RELEASE_TAG%...
    gh release create %RELEASE_TAG% --title "%REPO_NAME% %RELEASE_TAG%" --notes "Initial public release: four-layer causal-evidence plasma immunome atlas (9 main + 70 supplementary figures, evidence-tiered targets, manuscript)."
    echo(
    echo DONE. Repository and release published.
) else (
    echo GitHub CLI ^(gh^) not found.
    git remote get-url origin >nul 2>&1
    if errorlevel 1 (
        echo(
        echo No 'origin' remote set. To finish, either:
        echo   A^) Install GitHub CLI ^(https://cli.github.com^), then re-run this file, OR
        echo   B^) Create an empty repo on github.com and run:
        echo        git remote add origin https://github.com/^<you^>/%REPO_NAME%.git
        echo        git push -u origin main
    ) else (
        echo Pushing to existing origin remote...
        git push -u origin main
        echo DONE. Pushed to origin.
    )
)

:end
echo(
pause
