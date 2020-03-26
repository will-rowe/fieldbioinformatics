import vcf
import sys
from operator import attrgetter
from collections import defaultdict
from .vcftagprimersites import read_bed_file

def check_filter(v):
   total_reads = float(v.INFO['TotalReads'])
   qual = v.QUAL
   strandbias = float(v.INFO['StrandFisherTest'])

   if qual / total_reads <= 4:
       return False

   strand_fraction_by_strand = v.INFO['SupportFractionByStrand']
   if float(strand_fraction_by_strand[0]) < 0.5: 
       return False

   if float(strand_fraction_by_strand[1]) < 0.5:
       return False

   #if strandbias >= 100:
   #    return False

   if total_reads < 20:
       return False

   return True

def go(args):
   vcf_reader = vcf.Reader(filename=args.inputvcf)
   vcf_writer = vcf.Writer(open(args.outputvcf, 'w'), vcf_reader)

   for v in vcf_reader:
      if check_filter(v):
          vcf_writer.write_record(v)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('inputvcf')
    parser.add_argument('outputvcf')

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()


