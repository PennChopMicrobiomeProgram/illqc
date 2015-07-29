import argparse
import collections
import json
import os
import subprocess
import re
import csv
import StringIO

from .version import __version__

default_config = {
    "trimmomatic_jar_fp": "trimmomatic-0.33.jar",
    "adapter_dir": "",
    "adapter": "NexteraPE-PE",
    "leading": 3,
    "trailing": 3,
    "slidingwindow": (4, 15),
    "minlen": 36,
    "fastqc_dir":""
    }


def remove_file_ext(fp):
    return os.path.splitext(fp)[0]

def build_paired_fp(out_dir, old_fp):
    return os.path.join(out_dir, os.path.basename(old_fp))

class Trimmomatic(object):
    def __init__(self, config):
        self.config = config

    def make_command(self, fwd_fp, rev_fp, out_dir):
        fwd_paired_fp = build_paired_fp(out_dir, fwd_fp)
        rev_paired_fp = build_paired_fp(out_dir, rev_fp)
        
        fwd_unpaired_fp = os.path.join(
            out_dir, "%s_unpaired.fastq" % os.path.basename(remove_file_ext(fwd_fp)))
        rev_unpaired_fp = os.path.join(
            out_dir, "%s_unpaired.fastq" % os.path.basename(remove_file_ext(rev_fp)))
        
        trimlog_fp = os.path.join(out_dir, "trimmomatic.log")
        
        adapter_fp = os.path.join(
            self.config["adapter_dir"],
            "%s.fa" % self.config["adapter"])
        #check if the adapter file is there
        with open(adapter_fp) as f:
            pass
        
        return [
            "java", "-jar", self.config["trimmomatic_jar_fp"],
            "PE","-phred33",
            fwd_fp, rev_fp,
            fwd_paired_fp, fwd_unpaired_fp,
            rev_paired_fp, rev_unpaired_fp,
            "ILLUMINACLIP:%s:2:30:10:8:true" % adapter_fp,
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

class Fastqc(object):
    def __init__(sefl,config):
        sefl.config = config

    def make_command(self, fwd_fp, rev_fp, out_dir):
        return [
            self.config["fastqc_dir"],
            fwd_fp, rev_fp,
            '-extract', '-q', '-o', out_dir]
    def run(self, fwd_fp, rev_fp, out_dir):
        args = self.make_command(fwd_fp, rev_fp, out_dir)
        subprocess.check_call(args, stderr=subprocess.STDOUT)
        return (self.parse_fastqc_quality(make_report_fp(out_dir, fwd_fp)),
                self.parse_fastqc_quality(make_report_fp(out_dir, rev_fp)))

    @staticmethod
    def parse_fastqc_quality(output):
        with open(output) as f_in:
            report = f_in.read()
        tableString  = re.search('\>\>Per base sequence quality.*?\n(.*?)\n\>\>END_MODULE', report, re.DOTALL).group(1)
        table = parseTable(tableString)
        return table2dict(table, 0, 1)

def parseTable(tableString):
    f_s = StringIO.StringIO(tableString)
    reader = csv.reader(f_s, delimiter='\t')
    next(reader) #skip header
    return [row for row in reader]
                
def table2dict(table, keyIdx, valIdx):
    return {row[keyIdx]:row[valIdx] for row in table}

def make_report_fp(out_dir, file_dir):
    return os.path.join(out_dir, remove_file_ext(os.path.basename(file_dir))+'_fastqc', 'fastqc_data.txt')

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--forward-reads", required=True,
        type=argparse.FileType("r"),
        help="Forward reads file (FASTQ format)")
    p.add_argument("--reverse-reads", required=True,
        type=argparse.FileType("r"),
        help="Reverse reads file (FASTQ format)")
    p.add_argument("--output-dir", required=True,
        help="Output sequence data directory")
    p.add_argument("--qc-output-dir", required=True,
        help="Fastqc results directory")
    p.add_argument("--summary-file", required=True,
        type=argparse.FileType("w"),
        help="Summary file")
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

    before_trim_dir = os.path.join(args.qc_output_dir, 'before_trim')
    after_trim_dir = os.path.join(args.qc_output_dir, 'after_trim')

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    if not os.path.exists(before_trim_dir):
        os.makedirs(before_trim_dir) # make folders recursively
    if not os.path.exists(after_trim_dir):
        os.makedirs(after_trim_dir)
    
    app = Trimmomatic(config)
    summary_data_trim = app.run(fwd_fp, rev_fp, args.output_dir)

    fastqc = Fastqc(config)
    summary_qc_bef_fwd, summary_qc_bef_rev = fastqc.run(fwd_fp, rev_fp, before_trim_dir)

    summary_qc_aft_fwd, summary_qc_aft_rev = fastqc.run(build_paired_fp(args.output_dir, fwd_fp),
                                                        build_paired_fp(args.output_dir, rev_fp),
                                                        after_trim_dir)
    
    save_summary(args.summary_file, config,
                 build_summary(summary_data_trim, summary_qc_bef_fwd, summary_qc_bef_rev, summary_qc_aft_fwd, summary_qc_aft_rev))

def build_summary(trim, qc_bef_fwd, qc_bef_rev, qc_aft_fwd, qc_aft_rev):
    return {"illqc":trim, "fastqc_bef_fwd":qc_bef_fwd, "fastqc_bef_rev":qc_bef_rev, "fastqc_aft_fwd":qc_aft_fwd, "fastqc_aft_rev":qc_aft_rev}

def save_summary(f, config, data):
    result = {
        "program": "illqc",
        "version": __version__,
        "config": config,
        "data": data,
        }
    json.dump(result, f)
