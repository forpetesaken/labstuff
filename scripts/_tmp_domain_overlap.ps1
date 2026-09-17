$base = "c:\Users\Nat\Downloads\AIDEN Lab\Code\ctcf_bedmotif_analysis"
$thresholds = @('4.0','5.5','6.0','6.5')

function Get-AnchorRegions($bedpePath) {
    $anchors = New-Object System.Collections.Generic.List[object]
    Get-Content $bedpePath | ForEach-Object {
        if ([string]::IsNullOrWhiteSpace($_) -or $_.StartsWith('#')) { continue }
        $f = $_ -split "`t"
        if ($f.Count -lt 6) { continue }
        $anchors.Add([PSCustomObject]@{chrom=$f[0]; start=[int]$f[1]; end=[int]$f[2]})
        $anchors.Add([PSCustomObject]@{chrom=$f[3]; start=[int]$f[4]; end=[int]$f[5]})
    }
    return $anchors
}

function Get-MotifsByChrom($bedPath) {
    $map = @{}
    Get-Content $bedPath | ForEach-Object {
        if ([string]::IsNullOrWhiteSpace($_) -or $_.StartsWith('#')) { continue }
        $f = $_ -split "`t"
        if ($f.Count -lt 3) { continue }
        $chr = $f[0]
        if (-not $map.ContainsKey($chr)) { $map[$chr] = New-Object System.Collections.Generic.List[object] }
        $map[$chr].Add([PSCustomObject]@{start=[int]$f[1]; end=[int]$f[2]})
    }
    return $map
}

$config = @(
    @{Species='human'; MotifDir=Join-Path $base 'vertebrates\human_motifcalls\hs_fimo_out'; Bedpe='vertebrates\human_motifcalls\5000_blocks_top50.bedpe'},
    @{Species='arctic_lamprey'; MotifDir=Join-Path $base 'vertebrates\alamprey_motifcalls\ArcticLamprey_LJ_fimo_out'; Bedpe='vertebrates\alamprey_motifcalls\sl_25000_blocks_top50.bedpe'},
    @{Species='drosophila'; MotifDir=Join-Path $base 'invertebrates\drosophila_motifcalls\DM_fimo_out'; Bedpe='invertebrates\drosophila_motifcalls\dm_5000_blocks_top50.bedpe'}
)

$rows = New-Object System.Collections.Generic.List[object]
foreach ($c in $config) {
    $bedpePath = Join-Path $base $c.Bedpe
    if (!(Test-Path $bedpePath)) { continue }

    $anchors = Get-AnchorRegions $bedpePath
    $total = $anchors.Count
    if ($total -eq 0) { continue }

    foreach ($thr in $thresholds) {
        $motifPath = Join-Path $c.MotifDir ("fimo_ctcf_{0}.bed" -f $thr)
        if (!(Test-Path $motifPath)) { continue }

        $motifsByChrom = Get-MotifsByChrom $motifPath
        $withCtcf = 0

        foreach ($a in $anchors) {
            $hit = $false
            if ($motifsByChrom.ContainsKey($a.chrom)) {
                foreach ($m in $motifsByChrom[$a.chrom]) {
                    if (-not ($m.end -lt $a.start -or $m.start -gt $a.end)) { $hit = $true; break }
                }
            }
            if ($hit) { $withCtcf++ }
        }

        $rows.Add([PSCustomObject]@{
            Species = $c.Species
            Domain_Set = Split-Path $c.Bedpe -Leaf
            Threshold = $thr
            Boundaries_Total = $total
            Boundaries_With_CTCF = $withCtcf
            Percent_Overlap = [math]::Round((100.0 * $withCtcf / $total), 2)
        })
    }
}

$rows | Sort-Object Species, @{Expression={[double]$_.Threshold}} | Export-Csv -Path $outPath -NoTypeInformation -Delimiter "`t"
$rows | Sort-Object Species, @{Expression={[double]$_.Threshold}} | Format-Table -AutoSize
