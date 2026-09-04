$ErrorActionPreference = "Stop"

$siteDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = Split-Path -Parent $siteDir

$pythonExe = Join-Path $projectRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $pythonExe)) {
    $pythonExe = "python"
}

$builder = Join-Path $siteDir "build_site.py"
if (-not (Test-Path $builder)) {
    throw "Could not find builder script at $builder"
}

& $pythonExe $builder --project-root "$projectRoot" --output "index.html"

Write-Host "Updated $siteDir\index.html"
