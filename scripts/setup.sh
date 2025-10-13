#!/usr/bin/env bash
set -euo pipefail

DOCS=false
CLONE_DIR=""
REPO_URL="https://github.com/payoto/rsvs3d.git"

usage() {
  cat <<'USAGE'
Usage: setup.sh [options]

Options:
  --docs              Install documentation toolchain (Doxygen + Graphviz).
  --clone-dir <path>  Clone the repository into <path> if it is not already present.
  --repo-url <url>    Override the repository URL to clone (default: https://github.com/payoto/rsvs3d.git).
  -h, --help          Show this help and exit.
USAGE
}

log() {
  printf '[setup] %s\n' "$*"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --docs)
      DOCS=true
      shift
      ;;
    --clone-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for --clone-dir" >&2; exit 1; }
      CLONE_DIR=$2
      shift 2
      ;;
    --repo-url)
      [[ $# -ge 2 ]] || { echo "Missing value for --repo-url" >&2; exit 1; }
      REPO_URL=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if ! command -v git >/dev/null 2>&1; then
  echo "Error: git is required but not installed." >&2
  exit 1
fi

if [ -d ".git" ]; then
  REPO_ROOT=$(git rev-parse --show-toplevel)
  cd "$REPO_ROOT"
else
  TARGET_DIR=${CLONE_DIR:-rsvs3d}
  if [ -d "$TARGET_DIR/.git" ]; then
    cd "$TARGET_DIR"
  else
    log "Cloning $REPO_URL into $TARGET_DIR"
    git clone "$REPO_URL" "$TARGET_DIR"
    cd "$TARGET_DIR"
  fi
fi

log "Updating git submodules"
git submodule sync --recursive
ORIGIN_URL=$(git config --get remote.origin.url || echo "")
if [[ "$ORIGIN_URL" == https://* || -z "$ORIGIN_URL" ]]; then
  git config submodule.SRCC/modules/tetgen.url https://github.com/payoto/tetgen.git
  git config submodule.SRCC/modules/rsvs3d-externals.url https://github.com/payoto/rsvs3d-externals.git
fi
git submodule update --init --recursive --jobs 4

OS=$(uname -s)
case "$OS" in
  Linux)
    if command -v apt-get >/dev/null 2>&1; then
      log "Installing Linux build dependencies with apt-get"
      apt_get() {
        if [ "${EUID:-$(id -u)}" -eq 0 ]; then
          DEBIAN_FRONTEND=noninteractive apt-get "$@"
        elif command -v sudo >/dev/null 2>&1; then
          sudo DEBIAN_FRONTEND=noninteractive apt-get "$@"
        else
          echo "Error: run as root or install sudo to manage packages." >&2
          exit 1
        fi
      }
      apt_install() {
        apt_get install --no-install-recommends -y "$@"
      }
      apt_get update
      apt_install \
        build-essential g++-9 make cmake git \
        xorg-dev libglu1-mesa-dev xpra xserver-xorg-video-dummy freeglut3-dev \
        mesa-common-dev libgl1 libglfw3-dev libglew-dev libxi-dev libxmu-dev libxrandr-dev libxinerama-dev libxcursor-dev \
        pkg-config
      if $DOCS; then
        log "Installing Linux documentation toolchain"
        apt_install doxygen graphviz libclang-cpp9 libclang1-9
      fi
    else
      echo "Unsupported Linux distribution: apt-get not found." >&2
      exit 1
    fi
    ;;
  Darwin)
    if command -v brew >/dev/null 2>&1; then
      log "Installing macOS build dependencies with Homebrew"
      brew update
      brew install gcc make cmake git glfw glew pkg-config
      brew install --cask xquartz
      if $DOCS; then
        log "Installing macOS documentation toolchain"
        brew install doxygen graphviz
      fi
    else
      echo "Homebrew is required on macOS. Install it from https://brew.sh/ and re-run this script." >&2
      exit 1
    fi
    ;;
  MINGW*|MSYS*|CYGWIN*)
    if command -v powershell.exe >/dev/null 2>&1; then
      if $DOCS; then
        DOCS_FLAG="true"
      else
        DOCS_FLAG="false"
      fi
      log "Installing Windows build dependencies with Chocolatey"
      powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -Command "
\$docFlag = [bool]::Parse('$DOCS_FLAG');
if (-not (Get-Command choco -ErrorAction SilentlyContinue)) {
  Write-Error 'Chocolatey is required. Install it from https://chocolatey.org/install and re-run.';
  exit 1
}
choco install -y git make mingw cmake
if (\$docFlag) {
  choco install -y doxygen graphviz
}
" || {
        echo "Chocolatey-based installation failed." >&2
        exit 1
      }
    else
      echo "powershell.exe is required on Windows to install dependencies." >&2
      exit 1
    fi
    ;;
  *)
    echo "Unsupported operating system: $OS" >&2
    exit 1
    ;;
esac

log "Setup complete."
cat <<'NEXT'
Next steps:
  1. cd SRCC
  2. Run make targets such as:
       make all
       make testall
       make distribution
  3. To execute the test suite, use: ./RSVS3D --test --no-gui
If you enabled --docs, build API docs with: cd SRCC && doxygen docs/Doxyfile
NEXT
if [[ "$OS" == Linux ]] && grep -qi microsoft /proc/version 2>/dev/null; then
  cat <<'WSL'
WSL2 hint: remember to export DISPLAY and LIBGL_ALWAYS_INDIRECT as described in the README before launching the GUI.
WSL
fi
