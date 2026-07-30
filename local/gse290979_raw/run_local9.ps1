[CmdletBinding()]
param(
    [ValidateSet(
        "preflight",
        "reference",
        "download",
        "align",
        "benchmark",
        "evaluate",
        "cleanup-fastq"
    )]
    [string]$Phase = "preflight",
    [string]$SampleId = "",
    [string]$WorkRoot = "E:\smn2_gse290979_local9",
    [ValidateRange(1, 2)]
    [int]$Threads = 1,
    [ValidateRange(2000, 6000)]
    [int]$IndexMemoryMb = 2000,
    [ValidateRange(1, 8)]
    [int]$DownloadParts = 2,
    [ValidateRange(4, 64)]
    [int]$DownloadChunkMiB = 16,
    [ValidateRange(1, 50)]
    [int]$CurlRetries = 20,
    [ValidateSet("ena", "sra")]
    [string]$DownloadMethod = "ena",
    [string]$SraBinDir = ""
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) {
    throw "Rscript was not found at $Rscript"
}

if ([string]::IsNullOrWhiteSpace($SampleId)) {
    $BenchmarkPath = Join-Path $RepoRoot (
        "config\local_raw\GSE290979\local9\benchmark_sample.txt"
    )
    $SampleId = (Get-Content -LiteralPath $BenchmarkPath -TotalCount 1).Trim()
}

$null = New-Item -ItemType Directory -Path $WorkRoot -Force
$ResolvedWorkRoot = (Resolve-Path -LiteralPath $WorkRoot).Path
$DriveRoot = [System.IO.Path]::GetPathRoot($ResolvedWorkRoot)
$DriveName = $DriveRoot.TrimEnd("\").TrimEnd(":")
$Drive = Get-PSDrive -Name $DriveName
$FreeGiB = $Drive.Free / 1GB
$OutputFreeGiB = $FreeGiB
$SplitOutputStorage = $false
$OutputRootPath = Join-Path $ResolvedWorkRoot "output"
if (Test-Path -LiteralPath $OutputRootPath -PathType Container) {
    $OutputItem = Get-Item -LiteralPath $OutputRootPath
    if (
        $OutputItem.LinkType -eq "Junction" -and
        @($OutputItem.Target).Count -gt 0
    ) {
        $OutputTarget = @($OutputItem.Target)[0]
        $OutputDriveRoot = [System.IO.Path]::GetPathRoot($OutputTarget)
        $OutputDriveName = $OutputDriveRoot.TrimEnd("\").TrimEnd(":")
        $OutputDrive = Get-PSDrive -Name $OutputDriveName
        $OutputFreeGiB = $OutputDrive.Free / 1GB
        $SplitOutputStorage = $OutputDriveName -ne $DriveName
    }
}
$FreeGiBInvariant = $FreeGiB.ToString(
    "0.######",
    [System.Globalization.CultureInfo]::InvariantCulture
)
$OutputFreeGiBInvariant = $OutputFreeGiB.ToString(
    "0.######",
    [System.Globalization.CultureInfo]::InvariantCulture
)

$RequiredFreeGiB = switch ($Phase) {
    "preflight" { 15 }
    "reference" { 15 }
    "download" {
        if ($DownloadMethod -eq "sra") { 75 } else { 20 }
    }
    "align" { 30 }
    "benchmark" {
        if ($DownloadMethod -eq "sra") { 90 } else { 30 }
    }
    "evaluate" { 1 }
    "cleanup-fastq" { 1 }
}
if ($FreeGiB -lt $RequiredFreeGiB) {
    throw (
        "Only {0:N1} GiB is free on {1}; phase {2} requires at least {3} GiB." -f
        $FreeGiB, $DriveRoot, $Phase, $RequiredFreeGiB
    )
}
if (
    $SplitOutputStorage -and
    $Phase -in @("align", "benchmark") -and
    $OutputFreeGiB -lt 30
) {
    throw (
        "Only {0:N1} GiB is free on the split output drive; phase {1} " +
        "requires at least 30 GiB." -f $OutputFreeGiB, $Phase
    )
}

$OperatingSystem = Get-CimInstance -ClassName Win32_OperatingSystem
$FreeMemoryGiB = [double]$OperatingSystem.FreePhysicalMemory / 1MB
$MinimumFreeMemoryGiB = [Math]::Max(
    2.5,
    ($IndexMemoryMb / 1024.0) + 0.50
)
if (
    $Phase -in @("reference", "align", "benchmark") -and
    $FreeMemoryGiB -lt $MinimumFreeMemoryGiB
) {
    throw (
        (
            "Only {0:N1} GiB RAM is currently free; this {1} MB split index " +
            "requires at least {2:N1} GiB free."
        ) -f $FreeMemoryGiB, $IndexMemoryMb, $MinimumFreeMemoryGiB
    )
}
if ($FreeMemoryGiB -lt (($IndexMemoryMb / 1024.0) + 2.0)) {
    Write-Warning (
        (
            "Free RAM is {0:N1} GiB. The split index is configured for " +
            "{1} MB; closing other applications is recommended."
        ) -f $FreeMemoryGiB, $IndexMemoryMb
    )
}

Write-Host (
    (
        "Local9 phase={0} sample={1} method={2} work_root={3} " +
        "free_disk={4:N1}GiB output_free={5:N1}GiB " +
        "split_output={6} free_ram={7:N1}GiB"
    ) -f $Phase, $SampleId, $DownloadMethod, $ResolvedWorkRoot,
        $FreeGiB, $OutputFreeGiB, $SplitOutputStorage, $FreeMemoryGiB
)

$Arguments = @(
    (Join-Path $RepoRoot "r\23_gse290979_local_raw_pilot.R"),
    "--phase=$Phase",
    "--work-root=$ResolvedWorkRoot",
    "--sample-id=$SampleId",
    "--threads=$Threads",
    "--index-memory-mb=$IndexMemoryMb",
    "--download-parts=$DownloadParts",
    "--download-chunk-mib=$DownloadChunkMiB",
    "--curl-retries=$CurlRetries",
    "--free-disk-gib=$FreeGiBInvariant",
    "--output-free-disk-gib=$OutputFreeGiBInvariant",
    "--split-output-storage=$($SplitOutputStorage.ToString().ToLowerInvariant())",
    "--download-method=$DownloadMethod"
)
if (-not [string]::IsNullOrWhiteSpace($SraBinDir)) {
    $ResolvedSraBinDir = (Resolve-Path -LiteralPath $SraBinDir).Path
    $Arguments += "--sra-bin-dir=$ResolvedSraBinDir"
}
& $Rscript @Arguments
if ($LASTEXITCODE -ne 0) {
    throw "R local raw pilot failed with exit code $LASTEXITCODE"
}
