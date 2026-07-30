[CmdletBinding()]
param(
    [string]$SampleId = "",
    [string]$WorkRoot = "E:\smn2_gse290979_local9",
    [ValidateRange(4, 64)]
    [int]$DownloadChunkMiB = 16
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$SampleSheetPath = Join-Path $RepoRoot "config\GSE290979_local9_sample_sheet.tsv"
if ([string]::IsNullOrWhiteSpace($SampleId)) {
    $BenchmarkPath = Join-Path $RepoRoot (
        "config\local_raw\GSE290979\local9\benchmark_sample.txt"
    )
    $SampleId = (Get-Content -LiteralPath $BenchmarkPath -TotalCount 1).Trim()
}

$SampleSheet = Import-Csv -LiteralPath $SampleSheetPath -Delimiter "`t"
$Sample = $SampleSheet | Where-Object { $_.sample_id -eq $SampleId }
if ($null -eq $Sample) {
    throw "Sample is not in the frozen local9 profile: $SampleId"
}
if (@($Sample).Count -ne 1) {
    throw "Sample ID is not unique in the frozen local9 profile: $SampleId"
}
if (-not (Test-Path -LiteralPath $WorkRoot -PathType Container)) {
    throw "Work root does not exist: $WorkRoot"
}

$ResolvedWorkRoot = (Resolve-Path -LiteralPath $WorkRoot).Path
$FastqRoot = Join-Path $ResolvedWorkRoot "fastq\$SampleId"
$OutputRoot = Join-Path $ResolvedWorkRoot "output\$SampleId"
$ChunkBytes = [int64]$DownloadChunkMiB * 1MB
$ExpectedBytes = [int64]$Sample.fastq_bytes_total
$FastqUrls = @($Sample.fastq_r1_url, $Sample.fastq_r2_url)
$DurableBytes = [int64]0
$ObservedBytes = [int64]0
$CompletedChunks = 0
$RangeFiles = @()

foreach ($Url in $FastqUrls) {
    $Name = [System.IO.Path]::GetFileName(([uri]$Url).AbsolutePath)
    $FinalPath = Join-Path $FastqRoot $Name
    $PartialPath = "$FinalPath.partial"
    if (Test-Path -LiteralPath $FinalPath -PathType Leaf) {
        $Length = (Get-Item -LiteralPath $FinalPath).Length
        $DurableBytes += $Length
        $ObservedBytes += $Length
        continue
    }
    if (Test-Path -LiteralPath $PartialPath -PathType Leaf) {
        $Length = (Get-Item -LiteralPath $PartialPath).Length
        $DurableBytes += $Length
        $ObservedBytes += $Length
    }
    if (Test-Path -LiteralPath $FastqRoot -PathType Container) {
        $CurrentRanges = @(
            Get-ChildItem `
                -LiteralPath $FastqRoot `
                -Filter "$Name.partial.range_*" `
                -File
        )
        $RangeFiles += $CurrentRanges
        foreach ($Range in $CurrentRanges) {
            $ObservedBytes += $Range.Length
            if ($Range.Length -eq $ChunkBytes) {
                $DurableBytes += $Range.Length
                $CompletedChunks++
            }
        }
    }
}

$DownloadMarker = Join-Path $FastqRoot "DOWNLOAD_COMPLETE.tsv"
$AlignmentMarker = Join-Path $OutputRoot "ALIGNMENT_COMPLETE.tsv"
$DecisionMarker = Join-Path $OutputRoot "BENCHMARK_DECISION.tsv"
$Decision = ""
if (Test-Path -LiteralPath $DecisionMarker -PathType Leaf) {
    $DecisionTable = Import-Csv -LiteralPath $DecisionMarker -Delimiter "`t"
    $Decision = (
        $DecisionTable |
            Where-Object { $_.metric -eq "status" } |
            Select-Object -First 1
    ).value
}

$Processes = @(
    Get-CimInstance Win32_Process |
        Where-Object {
            $_.CommandLine -match [regex]::Escape($SampleId) -or
            $_.CommandLine -match "run_benchmark_supervisor"
        }
)
$ActiveDownload = @(
    $Processes | Where-Object { $_.Name -eq "curl.exe" }
).Count -gt 0
$WaitingSupervisors = @(
    $Processes | Where-Object {
        $_.Name -eq "powershell.exe" -and
        $_.CommandLine -match "run_benchmark_supervisor"
    }
).Count

$Phase = if (Test-Path -LiteralPath $DecisionMarker -PathType Leaf) {
    "benchmark_$($Decision.ToLowerInvariant())"
} elseif (Test-Path -LiteralPath $AlignmentMarker -PathType Leaf) {
    "alignment_complete_evaluation_pending"
} elseif (Test-Path -LiteralPath $DownloadMarker -PathType Leaf) {
    "download_verified_alignment_pending"
} elseif ($ActiveDownload) {
    "downloading"
} else {
    "download_inactive"
}

$LastDurableRange = $RangeFiles |
    Where-Object { $_.Length -eq $ChunkBytes } |
    Sort-Object LastWriteTime |
    Select-Object -Last 1
$DriveRoot = [System.IO.Path]::GetPathRoot($ResolvedWorkRoot)
$DriveName = $DriveRoot.TrimEnd("\").TrimEnd(":")
$Drive = Get-PSDrive -Name $DriveName
$OutputStoragePath = Join-Path $ResolvedWorkRoot "output"
$OutputStorageFreeGiB = $Drive.Free / 1GB
if (Test-Path -LiteralPath $OutputStoragePath -PathType Container) {
    $OutputStorageItem = Get-Item -LiteralPath $OutputStoragePath
    if (
        $OutputStorageItem.LinkType -eq "Junction" -and
        @($OutputStorageItem.Target).Count -gt 0
    ) {
        $OutputStoragePath = @($OutputStorageItem.Target)[0]
        $OutputStorageDriveRoot = [System.IO.Path]::GetPathRoot(
            $OutputStoragePath
        )
        $OutputStorageDriveName = $OutputStorageDriveRoot.
            TrimEnd("\").TrimEnd(":")
        $OutputStorageDrive = Get-PSDrive -Name $OutputStorageDriveName
        $OutputStorageFreeGiB = $OutputStorageDrive.Free / 1GB
    }
}
$CohortStatusPath = Join-Path $ResolvedWorkRoot (
    "cohort\LOCAL9_COHORT_STATUS.tsv"
)
$AnalysisCompletionPath = Join-Path $ResolvedWorkRoot (
    "cohort\LOCAL9_RAW_ANALYSIS_COMPLETE.tsv"
)
$AutopilotStatusPath = Join-Path $ResolvedWorkRoot (
    "cohort\LOCAL9_AUTOPILOT_STATUS.tsv"
)
$CohortRows = @()
if (Test-Path -LiteralPath $CohortStatusPath -PathType Leaf) {
    $CohortRows = @(
        Import-Csv -LiteralPath $CohortStatusPath -Delimiter "`t"
    )
}
$AutopilotState = ""
if (Test-Path -LiteralPath $AutopilotStatusPath -PathType Leaf) {
    $AutopilotRow = Import-Csv `
        -LiteralPath $AutopilotStatusPath `
        -Delimiter "`t" |
        Select-Object -First 1
    if ($null -ne $AutopilotRow) {
        $AutopilotState = $AutopilotRow.state
    }
}

[pscustomobject]@{
    sample_id = $SampleId
    run_accession = $Sample.run_accession
    phase = $Phase
    durable_mib = [math]::Round($DurableBytes / 1MB, 1)
    observed_mib = [math]::Round($ObservedBytes / 1MB, 1)
    expected_gib = [math]::Round($ExpectedBytes / 1GB, 3)
    durable_percent = [math]::Round(
        100 * $DurableBytes / $ExpectedBytes,
        2
    )
    completed_full_chunks = $CompletedChunks
    last_durable_chunk = if ($null -ne $LastDurableRange) {
        $LastDurableRange.Name
    } else {
        ""
    }
    last_durable_local = if ($null -ne $LastDurableRange) {
        $LastDurableRange.LastWriteTime.ToString("s")
    } else {
        ""
    }
    active_download = $ActiveDownload
    supervisor_processes = $WaitingSupervisors
    download_marker = Test-Path -LiteralPath $DownloadMarker -PathType Leaf
    alignment_marker = Test-Path -LiteralPath $AlignmentMarker -PathType Leaf
    benchmark_decision = $Decision
    cohort_status_marker = Test-Path `
        -LiteralPath $CohortStatusPath `
        -PathType Leaf
    cohort_passed_libraries = @(
        $CohortRows | Where-Object { $_.state -eq "PASS" }
    ).Count
    cohort_running_libraries = @(
        $CohortRows | Where-Object { $_.state -eq "RUNNING" }
    ).Count
    analysis_complete_marker = Test-Path `
        -LiteralPath $AnalysisCompletionPath `
        -PathType Leaf
    autopilot_state = $AutopilotState
    free_disk_gib = [math]::Round($Drive.Free / 1GB, 1)
    output_storage_path = $OutputStoragePath
    output_free_disk_gib = [math]::Round($OutputStorageFreeGiB, 1)
}
