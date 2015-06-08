import unittest

from illqc.main import Trimmomatic

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
        observed = app.make_command("a.fastq.gz", "b.fastq.gz", "mydir")
        expected = [
            'java', '-jar', 'trimmomatic-0.30.jar', 'PE', '--phred33',
            'a.fastq.gz', 'b.fastq.gz',
            'mydir/a.fastq.gz', 'mydir/b.fastq.gz',
            'mydir/a_unpaired.fastq.gz', 'mydir/b_unpaired.fastq.gz',
            'ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10',
            'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36',
            ]
        self.assertEqual(observed, expected)

