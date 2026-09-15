"""Restore optional historical artifacts from the exact conversation ZIP.

No network access, no execution of extracted code, no overwrite of changed files.
A digest match proves identity with the saved archive, not physical correctness.
"""
from pathlib import Path, PurePosixPath
import argparse, hashlib, stat, zipfile

EXPECTED = 'ef47e6bac6ac00f8df642ab90ad7f5b64b95df0237765d6a3f2e8afe633d66da'
PREFIX = 'curvature_predictors_and_state_audit/'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    payload = args.archive.read_bytes()
    if hashlib.sha256(payload).hexdigest() != EXPECTED:
        raise SystemExit('Archive SHA-256 mismatch; nothing extracted.')
    root = Path(__file__).resolve().parents[1]
    planned = []
    with zipfile.ZipFile(args.archive) as z:
        for entry in z.infolist():
            if entry.is_dir() or not entry.filename.startswith(PREFIX + 'legacy/'):
                continue
            relative = PurePosixPath(entry.filename[len(PREFIX):])
            if '..' in relative.parts or relative.is_absolute():
                raise SystemExit('Unsafe archive path; nothing extracted.')
            if stat.S_ISLNK(entry.external_attr >> 16):
                raise SystemExit('Symlink entry rejected; nothing extracted.')
            target = (root / relative).resolve()
            if not target.is_relative_to(root):
                raise SystemExit('Path escapes repository; nothing extracted.')
            data = z.read(entry)
            if target.exists():
                if not target.is_file() or target.read_bytes() != data:
                    raise SystemExit(f'Existing modified file: {relative}; nothing extracted.')
            else:
                planned.append((target, data))
    for target, data in planned:
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as handle:
            handle.write(data)
    print(f'Restored {len(planned)} historical files. No training or validation was run.')

if __name__ == '__main__':
    main()
