import os
import shutil
import tempfile
import unittest
import StringIO

from illqclib.main import (
    Trimmomatic, get_config,
)

class ConfigTests(unittest.TestCase):
        def setUp(self):
            self.temp_home_dir = tempfile.mkdtemp()
            self._old_home_dir = os.environ['HOME']
            os.environ['HOME'] = self.temp_home_dir

        def tearDown(self):
            shutil.rmtree(self.temp_home_dir)
            os.environ['HOME'] = self._old_home_dir

        def test_default_config_locataion(self):
            """Config file in user home dir should be read and used"""
            with open(os.path.join(self.temp_home_dir, ".illqc.json"), "w") as f:
                f.write('{"adapter_dir": "SOMECRAZYVALUE"}')
            config = get_config(None)
            self.assertEqual(config["adapter_dir"], u"SOMECRAZYVALUE")


class TrimmomaticTests(unittest.TestCase):
    config_vals = {
        "trimmomatic_jar_fp": "trimmomatic-0.30.jar",
        "adapter_dir": "adapters",
        "adapter": "NexteraPE-PE",
        "leading": 3,
        "trailing": 3,
        "slidingwindow": (4, 15),
        "minlen": 36,
        "java_heapsize":"200M"
    }

    def test_make_command(self):
        app = Trimmomatic(self.config_vals)
        observed = app.make_command("a.fastq", "b.fastq", "mydir")
        expected = [
            'java', '-Xmx200M', '-jar', 'trimmomatic-0.30.jar', 'PE', '-phred33',
            'a.fastq', 'b.fastq',
            'mydir/a.fastq', 'mydir/a_unpaired.fastq',
            'mydir/b.fastq', 'mydir/b_unpaired.fastq',
            'ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10:8:true',
            'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36',
            ]
        self.assertEqual(observed, expected)

    def test_make_command_sliding_window_as_list(self):
        config_vals = self.config_vals.copy()
        config_vals["slidingwindow"] = [6, 32]
        app = Trimmomatic(config_vals)
        observed = app.make_command("a.fastq", "b.fastq", "mydir")
        expected = [
            'java', '-Xmx200M', '-jar', 'trimmomatic-0.30.jar', 'PE', '-phred33',
            'a.fastq', 'b.fastq',
            'mydir/a.fastq', 'mydir/a_unpaired.fastq',
            'mydir/b.fastq', 'mydir/b_unpaired.fastq',
            'ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10:8:true',
            'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:6:32', 'MINLEN:36',
            ]
        self.assertEqual(observed, expected)
