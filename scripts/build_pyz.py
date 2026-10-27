#!/usr/bin/env python3
"""Build dist/b2s.pyz: the whole package as one file that runs with plain `python3` (>= 3.8).

    python3 scripts/build_pyz.py            # -> dist/b2s.pyz
    python3 dist/b2s.pyz matrix -i busco_results/ > matrix.tsv

Only the standard library is needed. Plots (matplotlib) and --xlsx (openpyxl) are used when
the running Python has them, and skipped with a warning otherwise.
"""

import shutil
import sys
import tempfile
import zipapp
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    out = ROOT / "dist" / "b2s.pyz"
    out.parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(ROOT / "src" / "busco2summary", Path(tmp) / "busco2summary",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        (Path(tmp) / "__main__.py").write_text(
            "import sys\nfrom busco2summary.cli import main\nsys.exit(main())\n")
        zipapp.create_archive(tmp, out, interpreter="/usr/bin/env python3", compressed=True)
    print(f"{out}  {out.stat().st_size / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
