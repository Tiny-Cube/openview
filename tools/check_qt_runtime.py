"""Report unresolved Linux Qt xcb libraries before starting a GUI process."""
from pathlib import Path
import subprocess
import sys

import PySide6


def main():
    if sys.platform != 'linux':
        raise SystemExit('This Qt xcb runtime check is for Linux.')
    plugin = Path(PySide6.__file__).parent / 'Qt/plugins/platforms/libqxcb.so'
    result = subprocess.run(['ldd', str(plugin)], capture_output=True, text=True, check=True)
    missing = [line.strip() for line in result.stdout.splitlines() if 'not found' in line]
    if missing:
        raise SystemExit('Missing Qt xcb libraries:\n' + '\n'.join(missing))
    print('Qt xcb library dependencies resolved; GUI/OpenGL checks can proceed.')


if __name__ == '__main__':
    main()
