"""Create empty FASTQ files and a metadata CSV for the demo.

Run once: uv run python examples/demo_fastqs/setup.py
"""

from pathlib import Path

HERE = Path(__file__).parent

SAMPLES = [
    ("WT_rep1", "S1"),
    ("WT_rep2", "S2"),
    ("WT_rep3", "S3"),
    ("KO_rep1", "S4"),
    ("KO_rep2", "S5"),
    ("KO_rep3", "S6"),
]

fastq_dir = HERE / "fastqs"
fastq_dir.mkdir(exist_ok=True)

for sample, index in SAMPLES:
    for read in ("R1", "R2"):
        (fastq_dir / f"{sample}_{index}_L001_{read}_001.fastq.gz").touch()

metadata = HERE / "metadata.csv"
metadata.write_text(
    "sample,strandedness,condition\n"
    + "\n".join(
        f"{name},reverse,{'wildtype' if name.startswith('WT') else 'knockout'}"
        for name, _ in SAMPLES
    )
    + "\n"
)

print(f"Created {len(SAMPLES) * 2} FASTQ files in {fastq_dir}")
print(f"Created metadata at {metadata}")
