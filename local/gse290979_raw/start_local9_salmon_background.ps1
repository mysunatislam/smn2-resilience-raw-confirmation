[CmdletBinding()]
param(
    [string]$WorkRoot = "E:\smn2_gse290979_local9",
    [ValidateRange(1, 2)]
    [int]$Threads = 2,
    [ValidateRange(1, 8)]
    [int]$DownloadParts = 8,
    [ValidateRange(4, 64)]
    [int]$DownloadChunkMiB = 16,
    [ValidateRange(1, 50)]
    [int]$CurlRetries = 20,
    [string]$WslDistro = "Ubuntu"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$Worker = Join-Path $PSScriptRoot "run_local9_salmon.ps1"
if (-not (Test-Path -LiteralPath $Worker -PathType Leaf)) {
    throw "Local9 Salmon worker not found: $Worker"
}

$Existing = @(
    Get-CimInstance Win32_Process |
        Where-Object {
            $_.Name -eq "powershell.exe" -and
            $_.CommandLine -match "run_local9_salmon[.]ps1.*-Phase all"
        }
)
if ($Existing.Count -gt 0) {
    throw (
        "A local9 Salmon all-phase worker is already active: " +
        (($Existing.ProcessId | Sort-Object) -join ",")
    )
}

$null = New-Item -ItemType Directory -Path $WorkRoot -Force
$ResolvedWorkRoot = (Resolve-Path -LiteralPath $WorkRoot).Path
$LogRoot = Join-Path $ResolvedWorkRoot "logs"
$null = New-Item -ItemType Directory -Path $LogRoot -Force
$RunId = Get-Date -Format "yyyyMMdd_HHmmss"
$StdOut = Join-Path $LogRoot "local9_salmon_$RunId.stdout.log"
$StdErr = Join-Path $LogRoot "local9_salmon_$RunId.stderr.log"
$PidPath = Join-Path $LogRoot "local9_salmon_$RunId.pid.tsv"

$Arguments = @(
    "-NoProfile",
    "-NonInteractive",
    "-ExecutionPolicy", "Bypass",
    "-File", $Worker,
    "-Phase", "all",
    "-WorkRoot", $ResolvedWorkRoot,
    "-Threads", $Threads,
    "-DownloadParts", $DownloadParts,
    "-DownloadChunkMiB", $DownloadChunkMiB,
    "-CurlRetries", $CurlRetries,
    "-WslDistro", $WslDistro
)
$Process = Start-Process `
    -FilePath "powershell.exe" `
    -ArgumentList $Arguments `
    -WindowStyle Hidden `
    -RedirectStandardOutput $StdOut `
    -RedirectStandardError $StdErr `
    -PassThru

[pscustomobject]@{
    status = "STARTED"
    process_id = $Process.Id
    phase = "all"
    work_root = $ResolvedWorkRoot
    threads = $Threads
    download_parts = $DownloadParts
    wsl_distribution = $WslDistro
    stdout_log = $StdOut
    stderr_log = $StdErr
    started_local = (Get-Date).ToString("s")
} | Export-Csv `
    -LiteralPath $PidPath `
    -Delimiter "`t" `
    -NoTypeInformation

Write-Host "Started expedited local9 Salmon workflow: PID $($Process.Id)"
Write-Host "stdout: $StdOut"
Write-Host "stderr: $StdErr"
Write-Host "pid:    $PidPath"
