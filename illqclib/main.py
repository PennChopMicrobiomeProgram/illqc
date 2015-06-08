import argparse
import collections
import json
import os
import subprocess

from .version import __version__

default_config = {
    "trimmomatic_jar_fp": "trimmomatic-0.33.jar",
    "adapter_dir": "",
    "adapter": "NexteraPE-PE",
    "leading": 3,
    "trailing": 3,
    "slidingwindow": (4, 15),
    "minlen": 36,
    }


def remove_file_ext(fp):
    return os.path.splitext(fp)[0]


class Trimmomatic(object):
    def __init__(self, config):
        self.config = config

    def make_command(self, fwd_fp, rev_fp, out_dir):
        fwd_paired_fp = os.path.join(out_dir, os.path.basename(fwd_fp))
        rev_paired_fp = os.path.join(out_dir, os.path.basename(rev_fp))

        fwd_unpaired_fp = os.path.join(
            out_dir,
            "%s_unpaired.fastq.gz" % remove_file_ext(remove_file_ext(fwd_fp)))
        rev_unpaired_fp = os.path.join(
            out_dir,
            "%s_unpaired.fastq.gz" % remove_file_ext(remove_file_ext(rev_fp)))

        trimlog_fp = os.path.join(out_dir, "trimmomatic.log")
        
        adapter_fp = os.path.join(
            self.config["adapter_dir"],
            "%s.fa" % self.config["adapter"])

        return [
            "java", "-jar", self.config["trimmomatic_jar_fp"],
            "PE","-phred33",
            fwd_fp, rev_fp,
            fwd_paired_fp, rev_paired_fp,
            fwd_unpaired_fp, rev_unpaired_fp,
            "ILLUMINACLIP:%s:2:30:10" % adapter_fp,
            "LEADING:%d" % self.config["leading"],
            "TRAILING:%d" % self.config["trailing"],
            "SLIDINGWINDOW:%d:%d" % self.config["slidingwindow"],
            "MINLEN:%d" % self.config["minlen"],
            ]

    def run(self, fwd_fp, rev_fp, out_dir):
        args = self.make_command(fwd_fp, rev_fp, out_dir)
        stdout = subprocess.check_output(args, stderr=subprocess.STDOUT)
        return self.parse_trim_summary(stdout)

    @staticmethod
    def parse_trim_summary(output):
        for line in output.splitlines():
            if line.startswith("Input Read Pairs"):
                # Really hacky: split up all words and select the
                # values we want. See test for expected output.
                toks = line.split()
                keys = ("input", "both kept", "fwd only", "rev only", "dropped")
                vals = (toks[3], toks[6], toks[11], toks[16], toks[19])
                vals = tuple(int(x) for x in vals)
                return dict(zip(keys, vals))
        raise ValueError(
            "Summary line not found in trimming output: %s" % output)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--forward-reads", required=True,
        type=argparse.FileType("r"),
        help="Forward read filepath (Gzipped FASTQ format)")
    p.add_argument("--reverse-reads", required=True,
        type=argparse.FileType("r"),
        help="Forward read filepath (Gzipped FASTQ format)")
    p.add_argument("--output-dir", required=True,
        help="Output sequence data directory")
    p.add_argument("--summary-file", required=True,
        type=argparse.FileType("w"),
        help="Summary filepath")
    p.add_argument("--config-file",
        type=argparse.FileType("r"),
        help="Configuration file (JSON format)")
    args = p.parse_args(argv)

    config = default_config
    if args.config_file:
        user_config = json.load(args.config_file)
        config.update(user_config)
    
    fwd_fp = args.forward_reads.name
    rev_fp = args.reverse_reads.name
    args.forward_reads.close()
    args.reverse_reads.close()

    if os.path.exists(args.output_dir):
        p.error("Output directory already exists")
    os.mkdir(args.output_dir)

    app = Trimmomatic(config)
    summary_data = app.run(fwd_fp, rev_fp, args.output_dir)
    save_summary(args.summary_file, config, summary_data)


def save_summary(f, config, data):
    result = {
        "program": "illqc",
        "version": __version__,
        "config": config,
        "data": data,
        }
    json.dump(result, f)
