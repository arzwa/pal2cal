#!/usr/bin/env python3
# coding: utf-8
# author: Arthur Zwaenepoel
# Because everybody writes a pal2nal at some point! Some people can't seem to
# get enough of them pal2nalling.
# Seriously though, I'd like a robust script that can handle a directory of
# alignments and a sequence database for the CDSs. Perhaps also a name map for
# the CDSs
import os
import sys
import click
import logging
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

__version__ = "0.1.0"
__author__ = "Arthur Zwaenepoel"

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stderr)
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

def pal2cal(aln, cds):
    l = aln.get_alignment_length()
    n = len(aln)
    recs = []
    for i in range(n):
        cds_seq = cds[aln[i].id].seq
        seq = ""
        k = 0
        for j in range(l):
            if aln[i].seq[j] == "-":
                seq += "---"
            elif aln[i].seq[j] == "X":
                seq += "XXX"
                k += 3
            else:
                seq += str(cds_seq[k:k+3])
                k += 3
        recs.append(SeqRecord(Seq(seq, generic_dna), id=aln[i].id))
    return MultipleSeqAlignment(recs)

def tryread(aln, fmt):
    try:
        return AlignIO.read(aln, fmt)
    except ValueError:
        logging.error("Could not read MSA `{}` (format `{}`)".format(aln, fmt))
        logging.error("Make sure the alignment *format* is correctly specified")
        exit()

@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('aln', type=click.Path(exists=True))
@click.argument('cds', type=click.Path(exists=True))
@click.option('--alnfmt', '-f', default="fasta", help="alignment format")
@click.option('--outdir', '-o', default=None, help='output directory')
def main(aln, cds, alnfmt, outdir):
    """
    Protein alignment ⟶  codon alignment
    """
    alndir = os.path.isdir(aln)
    cdsdir = os.path.isdir(cds)
    if alndir and not outdir:
        logging.error("Alignment input is a dir, but no output dir specified")
        logging.error("Please set the `--outdir/-o` option")
        exit()
    if cdsdir:
        cdsdata = {}
        for f in os.listdir(cds):
            pth = os.path.join(cds, f)
            logging.info("Indexing {}".format(f))
            cdsdata.update(SeqIO.to_dict(SeqIO.parse(pth, "fasta")))
    else:
        cdsdata = SeqIO.to_dict(SeqIO.parse(cds, "fasta"))
    if alndir:
        if os.path.isdir(outdir):
            logging.warning("Output directory {} exists".format(outdir))
        else:
            os.mkdir(outdir)
        for f in os.listdir(aln):
            pth = os.path.join(aln, f)
            out = os.path.join(outdir, f)
            pal = tryread(pth, alnfmt)
            nal = pal2cal(pal, cdsdata)
            AlignIO.write(nal, out, alnfmt)
    else:
        pal = tryread(aln, alnfmt)
        nal = pal2cal(pal, cdsdata)
        AlignIO.write(nal, sys.stdout, alnfmt)

if __name__ == "__main__":
    main()
