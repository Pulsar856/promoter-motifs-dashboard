$ErrorActionPreference = "Stop"
$projectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$venvPython = Join-Path $projectDir ".venv\Scripts\python.exe"

if (-not (Test-Path $venvPython)) {
    Write-Host "No local .venv found. Creating it now..."
    & (Join-Path $projectDir "setup_venv.ps1")
}

& $venvPython (Join-Path $projectDir "app.py")
