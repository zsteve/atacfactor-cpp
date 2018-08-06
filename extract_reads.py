#!/usr/bin/env python

import pysam
import re
import sys
import argparse

parser = argparse.ArgumentParser(description="usage: %prog [options] -p [peak positions bed files] -a [alignment sam/bam files] -o [<optional> output file] -s [set if bed files contain summit positions] -w [<optional> read window width]")
parser.add_argument("-p", "--peaks", nargs='+', help="bed files containing peak positions", required=True)
parser.add_argument("-a", "--alignments", nargs='+', help="sam/bam alignment files", required=True)
parser.add_argument("-o", "--output", nargs='?', type=argparse.FileType('w'), help="<optional> output file", default=sys.stdout)
parser.add_argument("-s", "--summit", action="store_true", help="set if bed files contain summit positions")
parser.add_argument("-w", "--width", help="<optional> set a custom read window width (100)", default=100)

args = parser.parse_args()

out = args.output

out.write("id\tchr\tpos\tstrand\tinsert\n")



for s in args.alignments:
    samfile = pysam.AlignmentFile(s, 'rb')
    refs = samfile.references
    lens = samfile.lengths
    for p in args.peaks:
        f = open(p, 'r')
        for line in f.readlines():
            sp = line.split()
            
            left = int(sp[1])
            right = int(sp[2])
            chrm = sp[0]
            if (args.summit):
                pos = int(sp[3])
                idn = sp[4]
            else:
                pos = left+(right-left)/2
                idn = sp[3]

            for read in samfile.fetch(chrm, max(0, left-args.width), min(right+args.width, lens[refs.index(chrm)])):
                if read.is_reverse:
                    strand="-"
                    start = read.reference_end-pos
                else:
                    strand = "+"
                    start = read.reference_start + 1 - pos
                out.write("%s\t%s\t%d\t%s\t%d\n" %(idn.split('_')[-1],chrm,start,strand,read.template_length)) 
        f.close()

out.close()
