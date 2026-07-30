[CmdletBinding()]
param(
    [ValidateSet(
        "preflight",
        "reference",
        "download",
        "index",
        "quantify",
        "quantify-verified",
        "analyze",
        "all"
    )]
    [string]$Phase = "preflight",
    [string]$WorkRoot = "E:\smn2_gse290979_local9",
    [ValidateRange(1, 2)]
    [int]$Threads = 2,
    [ValidateRange(1, 8)]
    [int]$DownloadParts = 8,
    [ValidateRange(4, 64)]
    [int]$DownloadChunkMiB = 16,
    [ValidateRange(1, 50)]
    [int]$CurlRetries = 20,
    [string]$WslDistro = "Ubuntu",
    [switch]$NoWaitForAlignment
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Rscript = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
$Local9Runner = Join-Path $PSScriptRoot "run_local9.ps1"
$SampleSheetPath = Join-Path $RepoRoot (
    "config\GSE290979_local9_sample_sheet.tsv"
)
$TranscriptomeName = "gencode.v47.transcripts.fa.gz"
$GencodeBaseUrl = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/" +
    "Gencode_human/release_47"
)

function Write-MetricMarker {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path,
        [Parameter(Mandatory = $true)]
        [System.Collections.IDictionary]$Metrics
    )
    $Lines = [System.Collections.Generic.List[string]]::new()
    $Lines.Add("metric`tvalue")
    foreach ($Entry in $Metrics.GetEnumerator()) {
        $Name = ([string]$Entry.Key).Replace("`t", " ").Replace("`n", " ")
        $Value = ([string]$Entry.Value).Replace("`t", " ").Replace("`n", " ")
        $Lines.Add("$Name`t$Value")
    }
    [System.IO.File]::WriteAllLines(
        $Path,
        $Lines,
        [System.Text.UTF8Encoding]::new($false)
    )
}

function Read-MarkerValue {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path,
        [Parameter(Mandatory = $true)]
        [string]$Metric
    )
    $Rows = Import-Csv -LiteralPath $Path -Delimiter "`t"
    $Row = @($Rows | Where-Object { $_.metric -eq $Metric })
    if ($Row.Count -ne 1) {
        throw "Metric '$Metric' is missing or duplicated in $Path"
    }
    return [string]$Row[0].value
}

function Test-PassMarker {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        return $false
    }
    return (Read-MarkerValue -Path $Path -Metric "status") -eq "PASS"
}

function Convert-ToWslPath {
    param([string]$Path)
    $Resolved = (Resolve-Path -LiteralPath $Path).Path
    if ($Resolved -notmatch "^([A-Za-z]):\\(.*)$") {
        throw "Cannot convert path to WSL form: $Resolved"
    }
    $Drive = $Matches[1].ToLowerInvariant()
    $Tail = $Matches[2].Replace("\", "/")
    return "/mnt/$Drive/$Tail"
}

function Invoke-WslSalmon {
    param([string[]]$Arguments)
    $PreviousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & wsl.exe -d $WslDistro -- salmon @Arguments
    $SalmonExitCode = $LASTEXITCODE
    $ErrorActionPreference = $PreviousErrorActionPreference
    if ($SalmonExitCode -ne 0) {
        throw (
            "Salmon failed with exit code ${SalmonExitCode}: " +
            ($Arguments -join " ")
        )
    }
}

function Invoke-RScript {
    param([string[]]$Arguments)
    & $Rscript @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw (
            "Rscript failed with exit code ${LASTEXITCODE}: " +
            ($Arguments -join " ")
        )
    }
}

function Get-NextAttemptDirectory {
    param([string]$BasePath)
    for ($Index = 1; $Index -le 999; $Index++) {
        $Candidate = "{0}.attempt-{1:D3}" -f $BasePath, $Index
        if (-not (Test-Path -LiteralPath $Candidate)) {
            return $Candidate
        }
    }
    throw "No free attempt directory remains for $BasePath"
}

function Get-SalmonQuantMutexName {
    param(
        [string]$WorkRoot,
        [string]$SampleId
    )
    $Hasher = [System.Security.Cryptography.SHA256]::Create()
    try {
        $Bytes = [System.Text.Encoding]::UTF8.GetBytes(
            "$WorkRoot|$SampleId"
        )
        $Hash = [System.BitConverter]::ToString(
            $Hasher.ComputeHash($Bytes)
        ).Replace("-", "")
    }
    finally {
        $Hasher.Dispose()
    }
    return "Local\SMN2_LOCAL9_SALMON_$($Hash.Substring(0, 32))"
}

function Test-CompletedSalmonQuantification {
    param(
        [string]$QuantMarker,
        [string]$SampleId
    )
    if (-not (Test-PassMarker -Path $QuantMarker)) {
        return $false
    }
    $ExistingPath = Read-MarkerValue `
        -Path $QuantMarker `
        -Metric "quant_path"
    if (
        -not (Test-Path -LiteralPath (
            Join-Path $ExistingPath "quant.sf"
        ) -PathType Leaf)
    ) {
        throw "Completed quant marker points to missing quant.sf"
    }
    Write-Host "Using completed Salmon quantification for $SampleId."
    return $true
}

function Wait-ForGenomicAlignment {
    if ($NoWaitForAlignment) {
        return
    }
    while ($true) {
        $Active = @(
            Get-CimInstance Win32_Process |
                Where-Object {
                    $_.Name -like "Rscript*" -and
                    $_.CommandLine -match (
                        "23_gse290979_local_raw_pilot[.]R.*" +
                        "--phase=(align|benchmark)"
                    )
                }
        )
        if ($Active.Count -eq 0) {
            return
        }
        Write-Host (
            "Waiting for active genomic alignment PID(s) {0} before " +
            "building/running Salmon..." -f
            (($Active.ProcessId | Sort-Object) -join ",")
        )
        Start-Sleep -Seconds 120
    }
}

if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) {
    throw "Rscript was not found at $Rscript"
}
if (-not (Test-Path -LiteralPath $Local9Runner -PathType Leaf)) {
    throw "Local9 downloader was not found at $Local9Runner"
}
if (-not (Test-Path -LiteralPath $SampleSheetPath -PathType Leaf)) {
    throw "Frozen sample sheet was not found at $SampleSheetPath"
}

$null = New-Item -ItemType Directory -Path $WorkRoot -Force
$ResolvedWorkRoot = (Resolve-Path -LiteralPath $WorkRoot).Path
$SalmonRoot = Join-Path $ResolvedWorkRoot "salmon"
$ReferenceRoot = Join-Path (
    $ResolvedWorkRoot
) "reference\gencode_v47_salmon"
$QuantRoot = Join-Path $SalmonRoot "quant"
$CountRoot = Join-Path $SalmonRoot "gene_counts"
$null = New-Item -ItemType Directory -Path $SalmonRoot -Force
$null = New-Item -ItemType Directory -Path $ReferenceRoot -Force
$null = New-Item -ItemType Directory -Path $QuantRoot -Force
$null = New-Item -ItemType Directory -Path $CountRoot -Force
$SampleSheet = @(
    Import-Csv -LiteralPath $SampleSheetPath -Delimiter "`t" |
        Sort-Object { [int]$_.array_index }
)
if (
    $SampleSheet.Count -ne 9 -or
    @($SampleSheet.sample_id | Select-Object -Unique).Count -ne 9
) {
    throw "The frozen local9 sheet must contain nine unique libraries"
}

$PreviousErrorActionPreference = $ErrorActionPreference
$ErrorActionPreference = "Continue"
$SalmonVersionOutput = & wsl.exe `
    -d $WslDistro `
    -- salmon --version 2>$null
$SalmonVersionExitCode = $LASTEXITCODE
$ErrorActionPreference = $PreviousErrorActionPreference
$SalmonVersion = ($SalmonVersionOutput | Out-String).Trim()
if (
    $SalmonVersionExitCode -ne 0 -or
    $SalmonVersion -notmatch "salmon"
) {
    throw "Salmon is unavailable in WSL distribution '$WslDistro'"
}
Write-Host (
    "Local9 Salmon phase=$Phase work_root=$ResolvedWorkRoot " +
    "threads=$Threads distro=$WslDistro version='$SalmonVersion'"
)

function Invoke-ReferencePhase {
    $Md5Path = Join-Path $ReferenceRoot "MD5SUMS"
    $TranscriptomePath = Join-Path $ReferenceRoot $TranscriptomeName
    $PartialPath = "$TranscriptomePath.part"
    $MarkerPath = Join-Path $ReferenceRoot "REFERENCE_COMPLETE.tsv"

    if (Test-PassMarker -Path $MarkerPath) {
        $Observed = (
            Get-FileHash -LiteralPath $TranscriptomePath -Algorithm MD5
        ).Hash.ToLowerInvariant()
        $Expected = Read-MarkerValue -Path $MarkerPath -Metric "md5"
        if ($Observed -ne $Expected) {
            throw "Existing transcriptome no longer matches its locked MD5"
        }
        Write-Host "Using checksum-verified GENCODE transcriptome."
        return
    }

    & curl.exe --fail --location --retry $CurlRetries `
        --output $Md5Path "$GencodeBaseUrl/MD5SUMS"
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to download official GENCODE MD5SUMS"
    }
    $Md5Text = Get-Content -LiteralPath $Md5Path -Raw
    $Match = [regex]::Match(
        $Md5Text,
        (
            "(?im)^([0-9a-f]{32})\s+\*?" +
            [regex]::Escape($TranscriptomeName) +
            "\s*$"
        )
    )
    if (-not $Match.Success) {
        throw "GENCODE MD5SUMS does not contain $TranscriptomeName"
    }
    $ExpectedMd5 = $Match.Groups[1].Value.ToLowerInvariant()

    if (Test-Path -LiteralPath $TranscriptomePath -PathType Leaf) {
        $ObservedMd5 = (
            Get-FileHash -LiteralPath $TranscriptomePath -Algorithm MD5
        ).Hash.ToLowerInvariant()
        if ($ObservedMd5 -ne $ExpectedMd5) {
            throw (
                "Existing transcriptome has MD5 $ObservedMd5; expected " +
                "$ExpectedMd5. It was preserved for review."
            )
        }
    } else {
        & curl.exe --fail --location --retry $CurlRetries `
            --continue-at - --output $PartialPath `
            "$GencodeBaseUrl/$TranscriptomeName"
        if ($LASTEXITCODE -ne 0) {
            throw "Failed to download $TranscriptomeName"
        }
        $ObservedMd5 = (
            Get-FileHash -LiteralPath $PartialPath -Algorithm MD5
        ).Hash.ToLowerInvariant()
        if ($ObservedMd5 -ne $ExpectedMd5) {
            throw (
                "Downloaded transcriptome has MD5 $ObservedMd5; expected " +
                "$ExpectedMd5. The partial file was preserved."
            )
        }
        Move-Item -LiteralPath $PartialPath -Destination $TranscriptomePath
    }

    Write-MetricMarker -Path $MarkerPath -Metrics ([ordered]@{
        status = "PASS"
        source = "$GencodeBaseUrl/$TranscriptomeName"
        transcriptome = $TranscriptomePath
        md5 = $ExpectedMd5
        bytes = (Get-Item -LiteralPath $TranscriptomePath).Length
        gencode_release = "47"
        verified_utc = (Get-Date).ToUniversalTime().ToString("o")
    })
    Write-Host "Checksum-locked $TranscriptomeName."
}

function Invoke-DownloadPhase {
    foreach ($Sample in $SampleSheet) {
        $SampleId = [string]$Sample.sample_id
        $Marker = Join-Path (
            $ResolvedWorkRoot
        ) "fastq\$SampleId\DOWNLOAD_COMPLETE.tsv"
        if (Test-PassMarker -Path $Marker) {
            Write-Host "Using checksum-verified FASTQs for $SampleId."
            continue
        }
        Write-Host "Downloading and verifying FASTQs for $SampleId..."
        & $Local9Runner `
            -Phase download `
            -SampleId $SampleId `
            -WorkRoot $ResolvedWorkRoot `
            -DownloadMethod ena `
            -DownloadParts $DownloadParts `
            -DownloadChunkMiB $DownloadChunkMiB `
            -CurlRetries $CurlRetries
        if ($LASTEXITCODE -ne 0) {
            throw "FASTQ download failed for $SampleId"
        }
        if (-not (Test-PassMarker -Path $Marker)) {
            throw "Verified download marker was not created for $SampleId"
        }
    }
}

function Invoke-IndexPhase {
    Wait-ForGenomicAlignment
    $ReferenceMarker = Join-Path $ReferenceRoot "REFERENCE_COMPLETE.tsv"
    if (-not (Test-PassMarker -Path $ReferenceMarker)) {
        throw "Run the reference phase before building the Salmon index"
    }
    $IndexMarker = Join-Path $ReferenceRoot "SALMON_INDEX_COMPLETE.tsv"
    if (Test-PassMarker -Path $IndexMarker) {
        $IndexPath = Read-MarkerValue -Path $IndexMarker -Metric "index_path"
        if (-not (Test-Path -LiteralPath $IndexPath -PathType Container)) {
            throw "Locked Salmon index directory is missing: $IndexPath"
        }
        Write-Host "Using completed Salmon index: $IndexPath"
        return
    }

    $TranscriptomePath = Join-Path $ReferenceRoot $TranscriptomeName
    $IndexPath = Get-NextAttemptDirectory -BasePath (
        Join-Path $ReferenceRoot "salmon_index_v47"
    )
    $null = New-Item -ItemType Directory -Path $IndexPath -Force
    Invoke-WslSalmon -Arguments @(
        "index",
        "-t", (Convert-ToWslPath -Path $TranscriptomePath),
        "-i", (Convert-ToWslPath -Path $IndexPath),
        "-p", [string]$Threads,
        "--gencode"
    )
    $VersionInfo = Join-Path $IndexPath "versionInfo.json"
    if (-not (Test-Path -LiteralPath $VersionInfo -PathType Leaf)) {
        throw "Salmon index completed without versionInfo.json"
    }
    Write-MetricMarker -Path $IndexMarker -Metrics ([ordered]@{
        status = "PASS"
        index_path = $IndexPath
        transcriptome = $TranscriptomePath
        salmon_version = $SalmonVersion
        index_kind = "gencode_transcriptome_only"
        decoy_aware = "FALSE"
        threads = $Threads
        completed_utc = (Get-Date).ToUniversalTime().ToString("o")
    })
}

function Invoke-QuantifyPhase {
    param([switch]$VerifiedOnly)

    Wait-ForGenomicAlignment
    $IndexMarker = Join-Path $ReferenceRoot "SALMON_INDEX_COMPLETE.tsv"
    if (-not (Test-PassMarker -Path $IndexMarker)) {
        throw "Run the index phase before Salmon quantification"
    }
    $IndexPath = Read-MarkerValue -Path $IndexMarker -Metric "index_path"
    $IndexWsl = Convert-ToWslPath -Path $IndexPath

    foreach ($Sample in $SampleSheet) {
        $SampleId = [string]$Sample.sample_id
        $RunAccession = [string]$Sample.run_accession
        $DownloadMarker = Join-Path (
            $ResolvedWorkRoot
        ) "fastq\$SampleId\DOWNLOAD_COMPLETE.tsv"
        if (-not (Test-PassMarker -Path $DownloadMarker)) {
            if ($VerifiedOnly) {
                Write-Host (
                    "Skipping $SampleId because its FASTQs are not yet " +
                    "checksum-verified."
                )
                continue
            }
            throw "Verified FASTQs are missing for $SampleId"
        }
        $R1 = Join-Path (
            $ResolvedWorkRoot
        ) "fastq\$SampleId\${RunAccession}_1.fastq.gz"
        $R2 = Join-Path (
            $ResolvedWorkRoot
        ) "fastq\$SampleId\${RunAccession}_2.fastq.gz"
        if (
            -not (Test-Path -LiteralPath $R1 -PathType Leaf) -or
            -not (Test-Path -LiteralPath $R2 -PathType Leaf)
        ) {
            throw "FASTQ pair is incomplete for $SampleId"
        }

        $SampleQuantRoot = Join-Path $QuantRoot $SampleId
        $null = New-Item -ItemType Directory -Path $SampleQuantRoot -Force
        $QuantMarker = Join-Path (
            $SampleQuantRoot
        ) "SALMON_QUANT_COMPLETE.tsv"
        if (
            Test-CompletedSalmonQuantification `
                -QuantMarker $QuantMarker `
                -SampleId $SampleId
        ) {
            continue
        }

        $MutexName = Get-SalmonQuantMutexName `
            -WorkRoot $ResolvedWorkRoot `
            -SampleId $SampleId
        $QuantMutex = [System.Threading.Mutex]::new(
            $false,
            $MutexName
        )
        $HasQuantMutex = $false
        try {
            try {
                if ($VerifiedOnly) {
                    $HasQuantMutex = $QuantMutex.WaitOne(0)
                }
                else {
                    Write-Host (
                        "Waiting for the per-sample quantification lock " +
                        "for $SampleId..."
                    )
                    $HasQuantMutex = $QuantMutex.WaitOne()
                }
            }
            catch [System.Threading.AbandonedMutexException] {
                $HasQuantMutex = $true
            }
            if (-not $HasQuantMutex) {
                Write-Host (
                    "Skipping $SampleId because another Salmon worker " +
                    "is quantifying it."
                )
                continue
            }
            if (
                Test-CompletedSalmonQuantification `
                    -QuantMarker $QuantMarker `
                    -SampleId $SampleId
            ) {
                continue
            }

            $AttemptPath = Get-NextAttemptDirectory -BasePath (
                Join-Path $SampleQuantRoot "salmon_quant"
            )
            $null = New-Item -ItemType Directory -Path $AttemptPath -Force
            Write-Host "Quantifying $SampleId into $AttemptPath..."
            Invoke-WslSalmon -Arguments @(
                "quant",
                "-i", $IndexWsl,
                "-l", "A",
                "-1", (Convert-ToWslPath -Path $R1),
                "-2", (Convert-ToWslPath -Path $R2),
                "--validateMappings",
                "--seqBias",
                "--gcBias",
                "--threads", [string]$Threads,
                "-o", (Convert-ToWslPath -Path $AttemptPath)
            )
            $QuantFile = Join-Path $AttemptPath "quant.sf"
            $MetaPath = Join-Path $AttemptPath "aux_info\meta_info.json"
            if (
                -not (Test-Path -LiteralPath $QuantFile -PathType Leaf) -or
                (Get-Item -LiteralPath $QuantFile).Length -eq 0 -or
                -not (Test-Path -LiteralPath $MetaPath -PathType Leaf)
            ) {
                throw "Salmon outputs are incomplete for $SampleId"
            }
            $Meta = Get-Content -LiteralPath $MetaPath -Raw |
                ConvertFrom-Json
            Write-MetricMarker -Path $QuantMarker -Metrics ([ordered]@{
                status = "PASS"
                sample_id = $SampleId
                run_accession = $RunAccession
                quant_path = $AttemptPath
                quant_file = $QuantFile
                salmon_version = $SalmonVersion
                library_type = [string]$Meta.library_types
                num_processed = [string]$Meta.num_processed
                num_mapped = [string]$Meta.num_mapped
                percent_mapped = [string]$Meta.percent_mapped
                validate_mappings = "TRUE"
                sequence_bias_correction = "TRUE"
                gc_bias_correction = "TRUE"
                whole_genome_alignment = "FALSE"
                completed_utc = (Get-Date).ToUniversalTime().ToString("o")
            })
        }
        finally {
            if ($HasQuantMutex) {
                $QuantMutex.ReleaseMutex()
            }
            $QuantMutex.Dispose()
        }
    }
}

function Invoke-AnalyzePhase {
    $ResultRoot = Join-Path (
        $RepoRoot
    ) "results\r\raw_confirmation\local9_salmon"
    $null = New-Item -ItemType Directory -Path $ResultRoot -Force

    Invoke-RScript -Arguments @(
        (Join-Path $RepoRoot "r\27_import_gse290979_local9_salmon.R"),
        "--work-root=$ResolvedWorkRoot",
        "--quant-root=$QuantRoot",
        "--count-root=$CountRoot",
        (
            "--transcriptome=" +
            (Join-Path $ReferenceRoot $TranscriptomeName)
        )
    )
    Invoke-RScript -Arguments @(
        (Join-Path $RepoRoot "r\24_analyze_gse290979_local9_counts.R"),
        "--work-root=$ResolvedWorkRoot",
        "--count-root=$CountRoot",
        "--decision-marker=SALMON_QUANT_COMPLETE.tsv",
        "--quantification-method=salmon_transcriptome_quasimapping",
        "--output-root=$ResultRoot"
    )
    Invoke-RScript -Arguments @(
        (Join-Path $RepoRoot "r\25_compare_gse290979_local9_to_processed.R"),
        "--raw-output-root=$ResultRoot",
        "--output-root=$ResultRoot"
    )
    Invoke-RScript -Arguments @(
        (Join-Path $RepoRoot "r\26_validate_gse290979_local9_outputs.R"),
        "--output-root=$ResultRoot",
        (
            "--completion-marker=" +
            "LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv"
        ),
        (
            "--expected-quantification-method=" +
            "salmon_transcriptome_quasimapping"
        ),
        "--expected-whole-genome-alignment=FALSE"
    )
    $Completion = Join-Path (
        $ResultRoot
    ) "LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv"
    if (-not (Test-Path -LiteralPath $Completion -PathType Leaf)) {
        throw "Salmon analysis completion marker was not created"
    }
    Write-MetricMarker `
        -Path (Join-Path $SalmonRoot "LOCAL9_SALMON_PIPELINE_COMPLETE.tsv") `
        -Metrics ([ordered]@{
            status = "COMPLETE"
            libraries = 9
            quantification_method = "salmon_transcriptome_quasimapping"
            whole_genome_alignment = "FALSE"
            result_root = $ResultRoot
            completion_marker = $Completion
            completed_utc = (Get-Date).ToUniversalTime().ToString("o")
        })
}

switch ($Phase) {
    "preflight" {
        Write-Host "Preflight passed for the expedited local9 Salmon workflow."
    }
    "reference" {
        Invoke-ReferencePhase
    }
    "download" {
        Invoke-DownloadPhase
    }
    "index" {
        Invoke-IndexPhase
    }
    "quantify" {
        Invoke-QuantifyPhase
    }
    "quantify-verified" {
        Invoke-ReferencePhase
        Invoke-IndexPhase
        Invoke-QuantifyPhase -VerifiedOnly
    }
    "analyze" {
        Invoke-AnalyzePhase
    }
    "all" {
        Invoke-ReferencePhase
        Invoke-DownloadPhase
        Invoke-IndexPhase
        Invoke-QuantifyPhase
        Invoke-AnalyzePhase
    }
}
