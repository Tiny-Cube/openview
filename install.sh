#!/bin/sh
set -eu
APP_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
command -v uv >/dev/null 2>&1 || {
    printf '%s\n' 'uv is required. Install uv, then run this script again.' >&2
    exit 1
}
uv sync --project "$APP_ROOT" --frozen
printf '\n%s\n' 'Installed OpenView. Start with:' "  sh \"$APP_ROOT/run.sh\"" 'For X11/Wayland, use the system GPU driver with OpenGL 3.3 or newer.'
