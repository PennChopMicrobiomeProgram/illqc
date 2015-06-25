import unittest

from illqclib.main import Trimmomatic

class TrimmomaticTests(unittest.TestCase):
    def test_make_command(self):
        app = Trimmomatic({
            "trimmomatic_jar_fp": "trimmomatic-0.30.jar",
            "adapter_dir": "adapters",
            "adapter": "NexteraPE-PE",
            "leading": 3,
            "trailing": 3,
            "slidingwindow": (4, 15),
            "minlen": 36,
            })
        observed = app.make_command("a.fastq", "b.fastq", "mydir")
        expected = [
            'java', '-jar', 'trimmomatic-0.30.jar', 'PE', '-phred33',
            'a.fastq', 'b.fastq',
            'mydir/a.fastq', 'mydir/a_unpaired.fastq',
            'mydir/b.fastq', 'mydir/b_unpaired.fastq',
            'ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10',
            'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36',
            ]
        self.assertEqual(observed, expected)

