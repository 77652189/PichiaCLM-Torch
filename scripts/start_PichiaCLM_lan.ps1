param([switch]$ConfigureFirewallOnly)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$OutputEncoding = [System.Text.UTF8Encoding]::new($false)

$projectName = "PichiaCLM"
$port = 8501
$ruleName = "PichiaCLM Streamlit 8501 LAN"
$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$appPath = "Model_PichiaCLM/interfaces/streamlit_app.py"

function Test-IsAdministrator {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Test-FirewallRuleReady {
    $rule = Get-NetFirewallRule -DisplayName $ruleName -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $rule -or $rule.Enabled -ne "True" -or $rule.Direction -ne "Inbound" -or $rule.Action -ne "Allow") {
        return $false
    }
    $portFilter = $rule | Get-NetFirewallPortFilter
    $addressFilter = $rule | Get-NetFirewallAddressFilter
    return ($portFilter.Protocol -eq "TCP" -and $portFilter.LocalPort -eq "$port" -and $addressFilter.RemoteAddress -eq "LocalSubnet")
}

function Ensure-FirewallRule {
    if (-not (Test-IsAdministrator)) {
        Write-Host "需要管理员权限配置 Windows 防火墙规则，正在请求 UAC..." -ForegroundColor Yellow
        $args = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", "`"$PSCommandPath`"", "-ConfigureFirewallOnly")
        Start-Process -FilePath "powershell.exe" -ArgumentList $args -Verb RunAs -Wait
        return
    }
    $rules = Get-NetFirewallRule -DisplayName $ruleName -ErrorAction SilentlyContinue
    if (-not $rules) {
        New-NetFirewallRule -DisplayName $ruleName -Direction Inbound -Action Allow -Protocol TCP -LocalPort $port -RemoteAddress LocalSubnet -Profile Any -Enabled True | Out-Null
    } else {
        $rules | Set-NetFirewallRule -Enabled True -Direction Inbound -Action Allow -Profile Any
        $rules | Get-NetFirewallPortFilter | Set-NetFirewallPortFilter -Protocol TCP -LocalPort $port
        $rules | Get-NetFirewallAddressFilter | Set-NetFirewallAddressFilter -RemoteAddress LocalSubnet
    }
}

function Get-ListeningProcessId {
    $connection = Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue |
        Select-Object -First 1
    if ($connection) {
        return $connection.OwningProcess
    }
    return $null
}

function Test-StreamlitHealth {
    try {
        $response = Invoke-WebRequest -Uri "http://127.0.0.1:$port/_stcore/health" -UseBasicParsing -TimeoutSec 5
        return ($response.StatusCode -eq 200 -and $response.Content.Trim() -eq "ok")
    } catch {
        return $false
    }
}

function Get-LanAddresses {
    return Get-NetIPAddress -AddressFamily IPv4 |
        Where-Object { $_.IPAddress -notlike "127.*" -and $_.IPAddress -notlike "169.254.*" -and $_.InterfaceAlias -notlike "*WSL*" } |
        Select-Object -ExpandProperty IPAddress
}

function Write-AccessUrls {
    param([string[]]$LanAddresses)

    Write-Host ""
    Write-Host "$projectName 局域网访问地址：" -ForegroundColor Green
    foreach ($address in $LanAddresses) {
        Write-Host "  http://$address`:$port" -ForegroundColor Green
    }
    Write-Host "  http://127.0.0.1:$port" -ForegroundColor Green
    Write-Host ""
}

if ($ConfigureFirewallOnly) {
    Ensure-FirewallRule
    exit
}

$lanAddresses = Get-LanAddresses
$existingPid = Get-ListeningProcessId
if ($existingPid) {
    if (Test-StreamlitHealth) {
        Write-Host ""
        Write-Host "$projectName 已经在运行，端口 $port 当前由 PID $existingPid 使用。" -ForegroundColor Green
        Write-AccessUrls -LanAddresses $lanAddresses
        Write-Host "这个启动器不需要再次启动服务；关闭本窗口不会停止已经运行的服务。" -ForegroundColor Yellow
        exit 0
    }

    $process = Get-CimInstance Win32_Process -Filter "ProcessId = $existingPid" -ErrorAction SilentlyContinue
    Write-Host ""
    Write-Host "端口 $port 已被占用，但健康检查没有通过，未启动新服务。" -ForegroundColor Red
    if ($process) {
        Write-Host "占用进程 PID: $existingPid" -ForegroundColor Yellow
        Write-Host "命令行: $($process.CommandLine)" -ForegroundColor Yellow
    }
    Write-Host "请先确认是否要停止该进程，或改用其他端口。" -ForegroundColor Yellow
    exit 1
}

if (-not (Test-FirewallRuleReady)) {
    Ensure-FirewallRule
}

Set-Location $projectRoot

$python = Join-Path $projectRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $python)) {
    $python = "python"
}

Write-Host ""
Write-Host "$projectName 将在局域网启动：" -ForegroundColor Green
Write-AccessUrls -LanAddresses $lanAddresses
Write-Host "如果要停止服务，关闭这个窗口或按 Ctrl+C。" -ForegroundColor Yellow
Write-Host ""

& $python -m streamlit run $appPath --server.address=0.0.0.0 --server.port=$port
