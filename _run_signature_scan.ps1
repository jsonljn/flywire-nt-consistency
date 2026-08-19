# Temporary runner: generate signature_scan.csv without Python.
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
if (-not $root) { $root = Get-Location }

$entropyPath = Join-Path $root "results\entropy_raw.csv"
$n10Path = Join-Path $root "results\entropy_raw_n10.csv"
if (Test-Path $n10Path) { $entropyPath = $n10Path }

$ntOrder = @("ACH","GABA","GLUT","DA","SER","OCT")
$histSeeds = @("R7","R8","R1-6")
$ornSeeds = @("ORN_V","ORN_VM3","ORN_VA2","ORN_DA3","ORN_DA4m","ORN_DA4l","ORN_DM2","ORN_DM3","ORN_DL4","ORN_DL3")
$dmSeeds = @("Dm12","Dm19","Dm1")
$patternSeeds = @{
  "histamine_blindspot" = $histSeeds
  "ORN_SER_confusion" = $ornSeeds
  "Dm_GLUT_confusion" = $dmSeeds
}
$jsFloor = 0.12
$looMargin = 1.35

function Parse-Counts([string]$raw) {
  $map = @{}
  if ([string]::IsNullOrWhiteSpace($raw)) { return $map }
  $t = $raw -replace "np\.int(?:64|32)\((\d+)\)", '$1'
  [regex]::Matches($t, "'([^']+)':\s*(\d+)") | ForEach-Object {
    $map[$_.Groups[1].Value] = [int]$_.Groups[2].Value
  }
  return $map
}

function To-Vec($counts) {
  $v = New-Object double[] 6
  $sum = 0.0
  for ($i=0; $i -lt 6; $i++) {
    $k = $ntOrder[$i]
    $c = 0.0
    if ($counts.ContainsKey($k)) { $c = [double]$counts[$k] }
    $v[$i] = $c
    $sum += $c
  }
  if ($sum -le 0) { return $v }
  for ($i=0; $i -lt 6; $i++) { $v[$i] = $v[$i] / $sum }
  return $v
}

function JS-Div([double[]]$p, [double[]]$q) {
  $eps = 1e-12
  $ps = 0.0; $qs = 0.0
  for ($i=0; $i -lt 6; $i++) { $ps += $p[$i]; $qs += $q[$i] }
  if ($ps -le 0 -or $qs -le 0) { return 1.0 }
  $kl = 0.0
  for ($i=0; $i -lt 6; $i++) {
    $pi = $p[$i] / $ps
    $qi = $q[$i] / $qs
    $m = 0.5 * ($pi + $qi)
    if ($pi -gt $eps) { $kl += 0.5 * $pi * [Math]::Log(($pi + $eps) / ($m + $eps)) / [Math]::Log(2) }
    if ($qi -gt $eps) { $kl += 0.5 * $qi * [Math]::Log(($qi + $eps) / ($m + $eps)) / [Math]::Log(2) }
  }
  return $kl
}

$rows = Import-Csv $entropyPath
$lookup = @{}
$records = @()
foreach ($r in $rows) {
  $counts = Parse-Counts $r.nt_distribution
  $vec = To-Vec $counts
  $lookup[$r.cell_type] = $vec
  $records += [pscustomobject]@{
    cell_type = $r.cell_type
    n_neurons = [int]$r.n_neurons
    entropy = [double]$r.entropy
    dominant_nt = $r.dominant_nt
    dominant_frac = [double]$r.dominant_frac
    vec = $vec
  }
}

$known = New-Object System.Collections.Generic.HashSet[string]
$threePath = Join-Path $root "results\three_confusion_patterns.csv"
if (Test-Path $threePath) {
  Import-Csv $threePath | ForEach-Object { [void]$known.Add($_.fafb_cell_type) }
}
$lit = New-Object System.Collections.Generic.HashSet[string]
$litPath = Join-Path $root "results\literature_validated_candidates.csv"
if (Test-Path $litPath) {
  Import-Csv $litPath | Where-Object { $_.agrees_with_literature -eq "True" } | ForEach-Object { [void]$lit.Add($_.fafb_cell_type) }
}

function Nearest-JS([double[]]$p, [string]$pattern, [string]$exclude) {
  $best = [double]::PositiveInfinity
  foreach ($seed in $patternSeeds[$pattern]) {
    if ($seed -eq $exclude) { continue }
    if (-not $lookup.ContainsKey($seed)) { continue }
    $d = JS-Div $p $lookup[$seed]
    if ($d -lt $best) { $best = $d }
  }
  return $best
}

$looMax = @{}
$thresholds = @{}
foreach ($pattern in $patternSeeds.Keys) {
  $dists = @()
  foreach ($seed in $patternSeeds[$pattern]) {
    if (-not $lookup.ContainsKey($seed)) { continue }
    $d = Nearest-JS $lookup[$seed] $pattern $seed
    if ([double]::IsInfinity($d) -eq $false) { $dists += $d }
  }
  $mx = if ($dists.Count -gt 0) { ($dists | Measure-Object -Maximum).Maximum } else { $jsFloor }
  $looMax[$pattern] = $mx
  $th = [Math]::Max($jsFloor, $mx * $looMargin)
  $thresholds[$pattern] = $th
  Write-Host ("{0}: LOO max={1:N4} threshold={2:N4}" -f $pattern, $mx, $th)
}

$out = New-Object System.Collections.Generic.List[object]
foreach ($rec in $records) {
  $name = $rec.cell_type
  $p = $rec.vec
  $jsH = Nearest-JS $p "histamine_blindspot" $name
  $jsO = Nearest-JS $p "ORN_SER_confusion" $name
  $jsD = Nearest-JS $p "Dm_GLUT_confusion" $name
  $bestPattern = "histamine_blindspot"
  $bestJs = $jsH
  if ($jsO -lt $bestJs) { $bestPattern = "ORN_SER_confusion"; $bestJs = $jsO }
  if ($jsD -lt $bestJs) { $bestPattern = "Dm_GLUT_confusion"; $bestJs = $jsD }
  $th = $thresholds[$bestPattern]
  $isSeed = $histSeeds.Contains($name) -or $ornSeeds.Contains($name) -or $dmSeeds.Contains($name)
  $inNb = $bestJs -le $th
  $already = $known.Contains($name)
  $novel = $inNb -and (-not $isSeed) -and (-not $already)
  $out.Add([pscustomobject]@{
    cell_type = $name
    n_neurons = $rec.n_neurons
    entropy = $rec.entropy
    dominant_nt = $rec.dominant_nt
    dominant_frac = $rec.dominant_frac
    p_ACH = $p[0]
    p_GABA = $p[1]
    p_GLUT = $p[2]
    p_DA = $p[3]
    p_SER = $p[4]
    p_OCT = $p[5]
    p_fast_transmitters = $p[0]+$p[1]+$p[2]
    js_histamine_blindspot = $jsH
    js_ORN_SER_confusion = $jsO
    js_Dm_GLUT_confusion = $jsD
    best_pattern = $bestPattern
    best_js = $bestJs
    pattern_threshold = $th
    in_neighborhood = $inNb
    is_seed = $isSeed
    already_name_matched = $already
    literature_confirmed = $lit.Contains($name)
    is_novel_candidate = $novel
  }) | Out-Null
}

$sorted = $out | Sort-Object @{Expression="is_novel_candidate";Descending=$true}, @{Expression="in_neighborhood";Descending=$true}, @{Expression="best_js";Descending=$false}
$csvPath = Join-Path $root "results\signature_scan.csv"
$sorted | Export-Csv -Path $csvPath -NoTypeInformation -Encoding UTF8
$novel = $sorted | Where-Object { $_.is_novel_candidate -eq $true }
$novel | Export-Csv -Path (Join-Path $root "results\signature_scan_novel.csv") -NoTypeInformation -Encoding UTF8

Write-Host ("Scored {0} types, novel candidates {1}" -f $sorted.Count, @($novel).Count)
Write-Host "Seeds:"
$sorted | Where-Object { $_.is_seed -eq $true } | Select-Object cell_type, best_pattern, best_js, in_neighborhood | Format-Table -AutoSize
Write-Host "Top novel:"
$novel | Select-Object -First 25 cell_type, n_neurons, entropy, dominant_nt, best_pattern, best_js | Format-Table -AutoSize
