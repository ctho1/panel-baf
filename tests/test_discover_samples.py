import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from discover_samples import discover  # noqa: E402


class DiscoverSamplesTest(unittest.TestCase):
    def test_pairs_and_skips(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ["A_S1_R1_001.fastq.gz", "A_S1_R2_001.fastq.gz", "B_S2_R1_001.fastq.gz", "notes.txt",
                         "C_R1_001_trim.fq.gz", "C_R2_001_trim.fq.gz"]:
                (Path(d) / name).touch()
            pairs, skipped = discover(d)
            self.assertEqual([p[0] for p in pairs], ["A_S1", "C_trim"])
            self.assertEqual([p.name for p in skipped], ["B_S2_R1_001.fastq.gz"])


if __name__ == "__main__":
    unittest.main()
