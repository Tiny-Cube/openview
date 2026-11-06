"""Check release contents and runtime resources without executing original software."""
import hashlib
import json
from email.parser import BytesParser
from pathlib import Path
import re
import tarfile
import zipfile


def main():
    root = Path(__file__).resolve().parents[1]
    version = re.search(r'^version = "([^"]+)"', (root / 'pyproject.toml').read_text(), re.M)[1]
    wheel = root / 'dist' / f'openview-{version}-py3-none-any.whl'
    source = root / 'dist' / f'openview-{version}.tar.gz'
    forbidden_parts = {'.git', '.venv', '.uv-cache', 'node_modules', 'original', 'analysis', '__pycache__'}
    forbidden_suffixes = {'.exe', '.dll', '.so', '.key', '.pem', '.chk', '.fchk', '.zip'}

    def allowed(name):
        path = Path(name)
        assert not path.is_absolute() and '..' not in path.parts, name
        assert not forbidden_parts.intersection(path.parts), name
        assert path.suffix.lower() not in forbidden_suffixes, name
        assert '.so.' not in path.name, name

    expected = {p.relative_to(root / 'src').as_posix(): p for p in (root / 'src/openview').rglob('*')
                if p.is_file() and '__pycache__' not in p.parts}
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        for name in names:
            allowed(name)
            assert not archive.read(name).startswith((b'\x7fELF', b'MZ', b'PK\x03\x04')), name
        runtime = {name for name in names if name.startswith('openview/') and not name.endswith('/')}
        assert runtime == expected.keys(), (runtime ^ expected.keys())
        for name, path in expected.items():
            assert archive.read(name) == path.read_bytes(), name
        metadata = BytesParser().parsebytes(archive.read(f'openview-{version}.dist-info/METADATA'))
        assert metadata['Name'] == 'openview' and metadata['Version'] == version
        entry = archive.read(f'openview-{version}.dist-info/entry_points.txt').decode()
        assert 'openview = openview.app:main' in entry
        catalog = json.loads(archive.read('openview/data/element_fragments.json'))
        assert catalog == json.loads((root / 'src/openview/data/element_fragments.json').read_text())
    with tarfile.open(source) as archive:
        for member in archive.getmembers():
            allowed(member.name)
            assert member.name.startswith(f'openview-{version}/') or (member.name == f'openview-{version}' and member.isdir())
            assert not member.issym() and not member.islnk(), member.name
            if member.isfile():
                assert not archive.extractfile(member).read(4).startswith((b'\x7fELF', b'MZ', b'PK\x03\x04')), member.name
        for relative in ['README.md', 'LICENSE', 'uv.lock', 'run.sh', 'install.sh',
                         'src/openview/app.py', 'src/openview/data/element_fragments.json',
                         'docs/VALIDATION.zh-CN.md', 'docs/validation/interaction.json',
                         'docs/images/openview.png', 'tests/test_element_fragments.py']:
            assert archive.extractfile(f'openview-{version}/{relative}').read() == (root / relative).read_bytes()
        source_members = len(archive.getmembers())
    report = dict(version=version, runtime_files=len(expected), source_archive_members=source_members,
                  metadata_author=metadata['Author-email'], console_entry='openview = openview.app:main',
                  no_original_binaries_or_keys=True,
                  archives={p.name: dict(bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                            for p in [wheel, source]})
    target = root / 'validation/distribution.json'
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
