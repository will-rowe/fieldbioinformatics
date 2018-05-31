import sys
from Bio import SeqIO
import tempfile
import os
import glob
import shutil
import pandas as pd

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

def run(parser, args):
	if args.guppy:
		d = '%s' % (args.directory)
	else:
		d = '%s/workspace/pass' % (args.directory)

	if args.prefix:
		prefix = args.prefix
	else:
		prefix = os.path.split(args.directory)[-1]

	all_fastq_outfn = "%s_all.fastq" % (prefix)
	all_fastq_outfh = open(all_fastq_outfn, "w")

	for root, dirs, files in os.walk(d):
		paths = os.path.split(root)
		barcode_directory = paths[-1]

		fastq = [root+'/'+f for f in files if f.endswith('.fastq')]
		if len(fastq):
			fastq_outfn = "%s_%s.fastq" % (prefix, barcode_directory)
			outfh = open(fastq_outfn, "w")
			print >>sys.stderr, "Processing %s files in %s" % (len(fastq), barcode_directory)

			dups = set()
			uniq = 0
			total = 0	
			for f in fastq:
				for rec in SeqIO.parse(open(f), "fastq"):
					if args.max_length and len(rec) > args.max_length:
						continue
					if args.min_length and len(rec) < args.min_length:
						continue

					total += 1
					if rec.id not in dups:
						SeqIO.write([rec], outfh, "fastq")
						SeqIO.write([rec], all_fastq_outfh, "fastq")

						dups.add(rec.id)
						uniq += 1

			outfh.close()

			print "%s\t%s\t%s" % (fastq_outfn, total, uniq)

        all_fastq_outfh.close()

	print >>sys.stderr, "Collecting summary files\n"

	dfs = []

	summary_outfn = "%s_sequencing_summary.txt" % (prefix)
	summaryfh = open(summary_outfn, "w")

	d = '%s/sequencing_summary*.txt' % (args.directory,)
	for summaryfn in glob.glob(d):
		df = pd.read_csv(summaryfn, sep="\t")
		dfs.append(df)

	pd.concat(dfs).to_csv(summaryfh, sep="\t", index=False)
	summaryfh.close()


