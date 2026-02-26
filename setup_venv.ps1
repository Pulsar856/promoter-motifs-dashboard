param(
    [switch]$Recreate
)

$ErrorActionPreference = "Stop"
$projectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$venvDir = Join-Path $projectDir ".venv"
$venvPython = Join-Path $venvDir "Scripts\python.exe"

if ($Recreate -and (Test-Path $venvDir)) {
    Remove-Item -Recurse -Force $venvDir
}

if (-not (Test-Path $venvPython)) {
    python -m venv $venvDir
}

& $venvPython -m pip install --upgrade pip
& $venvPython -m pip install -r (Join-Path $projectDir "requirements.txt")

Write-Host "Local virtual environment is ready:" $venvDir
Write-Host "Run the app with: .\\run.ps1"
