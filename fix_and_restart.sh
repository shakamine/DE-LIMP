#!/usr/bin/env bash
# =============================================================================
# fix_and_restart.sh — One-shot recovery script for stuck DE-LIMP installs
# =============================================================================
#
# Diagnoses + fixes the common stuck states from v3.10.x install hotfix train:
#   - Port 3838 occupied by a zombie R process
#   - WSL-side clone at ~/.delimp/DE-LIMP stuck on stale code
#   - .NET 8 SDK missing (only runtime installed)
#   - Browser tab labeled "Docker" instead of "WSL"
#
# Usage (inside Ubuntu / WSL):
#
#   curl -sSL https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/fix_and_restart.sh | bash
#
# Or if you've cloned the repo:
#
#   bash fix_and_restart.sh
#
# This script is read-only on user data and idempotent — safe to re-run.
# =============================================================================

set -uo pipefail

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log()  { echo -e "${GREEN}[fix]${NC} $1"; }
warn() { echo -e "${YELLOW}[fix]${NC} $1"; }
err()  { echo -e "${RED}[fix]${NC} $1" >&2; }

REPO_DIR="${HOME}/.delimp/DE-LIMP"
PORT="${DELIMP_PORT:-3838}"

echo
echo "================================================================"
echo "  DE-LIMP fix-and-restart script"
echo "================================================================"
echo

# -----------------------------------------------------------------------------
# 1. Sanity check — must be inside Linux/WSL
# -----------------------------------------------------------------------------
if [ "$(uname -s)" != "Linux" ]; then
    err "This script must run inside Linux (or WSL Ubuntu)."
    err "If you're on Windows PowerShell, open Ubuntu first:"
    err "  wsl -d Ubuntu"
    err "Then re-run."
    exit 1
fi

# -----------------------------------------------------------------------------
# 2. Kill stale DE-LIMP / Shiny processes
# -----------------------------------------------------------------------------
log "Looking for stale R / Shiny processes..."
killed_any=0
for pat in 'shiny::runApp' 'delimp_wsl_setup.sh run' "Rscript.*delimp"; do
    pids=$(pgrep -f "${pat}" 2>/dev/null || true)
    if [ -n "${pids}" ]; then
        log "  Killing PIDs matching '${pat}': ${pids}"
        # shellcheck disable=SC2086
        kill ${pids} 2>/dev/null || true
        sleep 1
        # shellcheck disable=SC2086
        kill -9 ${pids} 2>/dev/null || true
        killed_any=1
    fi
done

# Also clear anything still holding the port
port_pids=$(sudo -n lsof -ti ":${PORT}" 2>/dev/null || true)
if [ -n "${port_pids}" ]; then
    log "  Killing PIDs on port ${PORT}: ${port_pids}"
    # shellcheck disable=SC2086
    sudo -n kill -9 ${port_pids} 2>/dev/null || true
    killed_any=1
fi

if [ "${killed_any}" = "1" ]; then
    log "Stale processes cleaned up."
else
    log "No stale processes found."
fi

# -----------------------------------------------------------------------------
# 3. Check / update DE-LIMP repo
# -----------------------------------------------------------------------------
if [ ! -d "${REPO_DIR}/.git" ]; then
    err "DE-LIMP not installed at ${REPO_DIR}."
    err "Run the regular WSL setup first:"
    err "  curl -sSL https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/delimp_wsl_setup.sh -o ~/delimp_wsl_setup.sh"
    err "  chmod +x ~/delimp_wsl_setup.sh"
    err "  bash ~/delimp_wsl_setup.sh"
    exit 1
fi

log "Updating DE-LIMP repo at ${REPO_DIR}..."
cd "${REPO_DIR}"
if ! git pull --ff-only 2>/dev/null; then
    warn "Fast-forward failed; trying force-reset to origin/main..."
    git fetch --depth 1 origin main
    git reset --hard origin/main
fi
repo_version=$(cat VERSION 2>/dev/null | tr -d '[:space:]')
repo_commit=$(git rev-parse --short HEAD 2>/dev/null)
log "  Now at: DE-LIMP v${repo_version:-unknown} (${repo_commit:-unknown})"

# -----------------------------------------------------------------------------
# 4. Verify .NET 8 SDK present
# -----------------------------------------------------------------------------
log "Checking .NET 8 SDK..."
if ! command -v dotnet >/dev/null 2>&1; then
    err "dotnet not on PATH. Install .NET 8 SDK first:"
    err "  bash ~/delimp_wsl_setup.sh"
    exit 1
fi

if ! dotnet --list-sdks 2>/dev/null | grep -qE '^8\.'; then
    warn ".NET 8 SDK NOT installed (you may have only the runtime)."
    warn "DIA-NN 2.x's Thermo .raw reader needs the SDK."
    warn ""
    warn "Quick fix from Microsoft's official installer:"
    warn "  sudo /tmp/dotnet-install.sh --channel 8.0 --install-dir /usr/share/dotnet"
    warn "Then re-run this script."
    warn ""
    warn "Or install via apt (Ubuntu 24.04+):"
    warn "  sudo apt install -y dotnet-sdk-8.0"
    warn ""
    warn "Or via the project's setup script (which handles both paths):"
    warn "  bash ~/delimp_wsl_setup.sh"
else
    sdk_version=$(dotnet --list-sdks 2>/dev/null | grep -E '^8\.' | head -1 | awk '{print $1}')
    log "  ✓ .NET 8 SDK ${sdk_version} present"
fi

# -----------------------------------------------------------------------------
# 5. Verify DIA-NN binary
# -----------------------------------------------------------------------------
DIANN_BIN="${HOME}/.delimp/diann/diann-linux"
if [ -x "${DIANN_BIN}" ]; then
    log "  ✓ DIA-NN binary at ${DIANN_BIN}"
else
    warn "DIA-NN binary missing at ${DIANN_BIN}"
    warn "Re-run setup script: bash ~/delimp_wsl_setup.sh"
fi

# -----------------------------------------------------------------------------
# 6. Find a free port (default 3838, fall through to 3839, 3840, ...)
# -----------------------------------------------------------------------------
find_free_port() {
    local start_port="$1"
    for p in $(seq "${start_port}" $((start_port + 20))); do
        if ! sudo -n lsof -ti ":${p}" >/dev/null 2>&1; then
            # Also check via /proc/net/tcp in case lsof needs sudo
            if ! ss -ln 2>/dev/null | grep -qE ":${p}\s"; then
                echo "${p}"
                return 0
            fi
        fi
    done
    return 1
}

free_port=$(find_free_port "${PORT}")
if [ -z "${free_port}" ]; then
    err "Could not find a free port near ${PORT}."
    err "Try: sudo lsof -i :3838-3858 to see what's holding them."
    exit 1
fi

if [ "${free_port}" != "${PORT}" ]; then
    warn "Port ${PORT} is busy. Using port ${free_port} instead."
    warn "Open: http://localhost:${free_port}"
else
    log "Port ${PORT} is free. Starting on http://localhost:${PORT}"
fi

# -----------------------------------------------------------------------------
# 7. Launch DE-LIMP
# -----------------------------------------------------------------------------
echo
log "================================================================"
log "  Launching DE-LIMP v${repo_version:-unknown} on port ${free_port}"
log "  Browser: http://localhost:${free_port}"
log "  Press Ctrl+C in this terminal to stop the app."
log "================================================================"
echo

DELIMP_PORT="${free_port}" exec bash ~/delimp_wsl_setup.sh run
