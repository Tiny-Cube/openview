#!/bin/sh
set -eu
APP_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ ! -x "$APP_ROOT/.venv/bin/python" ]; then
    printf '%s\n' 'Run sh install.sh first.' >&2
    exit 1
fi
command -v uv >/dev/null 2>&1 || {
    printf '%s\n' 'uv is required to run OpenView.' >&2
    exit 1
}
# Prefer this checkout even if the environment contains an older installed build.
export PYTHONPATH="$APP_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
exec uv run --project "$APP_ROOT" --no-sync python -m openview "$@"
