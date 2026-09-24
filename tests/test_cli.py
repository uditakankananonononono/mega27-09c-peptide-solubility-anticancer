import subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "pepx.cli", *args],
                          capture_output=True, text=True, cwd=ROOT,
                          env={"PYTHONPATH": "src", "PATH": "/usr/bin:/bin:/usr/local/bin"})


def test_describe_outputs_descriptors():
    r = run_cli("describe", "KLAKLAKKLAKLAK")
    assert r.returncode == 0
    assert '"gravy"' in r.stdout and '"pI"' in r.stdout


def test_score_table_shape():
    r = run_cli("score", "KLAKLAKKLAKLAK", "GGGGGGGGGG")
    assert r.returncode == 0
    lines = [l for l in r.stdout.strip().splitlines() if l]
    assert len(lines) == 3  # header + 2 rows
    assert lines[0].startswith("sequence")


def test_score_file(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">a\nKLAKLAK\n>b\nFEKEAKKIEIKRH\n")
    out = tmp_path / "out.tsv"
    r = run_cli("score-file", str(fa), "--out", str(out))
    assert r.returncode == 0
    body = out.read_text().strip().splitlines()
    assert len(body) == 3
