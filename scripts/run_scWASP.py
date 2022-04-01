#!/usr/bin/env python3

import os
import glob
import pysam
import argparse
import subprocess
import logging
import pandas as pd
from multiprocessing import Pool

def split_sc_bams(args, singlecell):
	with pysam.AlignmentFile(args.bam_file, 'rb') as bamfile:
		hdr = bamfile.header.to_dict()
		del hdr['CO']
		obam = {bc:pysam.AlignmentFile(os.path.join(args.bam_output_path, '{}_{}.possorted.bam'.format(args.sample_name, bc)), 'wb', header=hdr) for bc in singlecell.index}
		for read in bamfile.fetch(until_eof=True):
			if read.has_tag('CB:Z') and read.mapping_quality>=10:
				bc = read.get_tag('CB:Z').split('-')[0]
				if bc in obam:
					obam[bc].write(read)
		for ob in obam:
			obam[ob].close()
	return args.bam_output_path

def create_snps_files(args):
	chr_pos = pd.DataFrame(columns=['#CHROM','POS'])
	for chrom in map(str, range(1,23)):
		vcf = pd.read_table(os.path.join('{}'.format(args.variants_dir), 'chr{}.dose.vcf.gz'.format(chrom)), sep='\t', header=0, skiprows=18)
		vcf['SAMPLE'] = vcf[args.sample_name].str.split(':', expand=True)[0]
		vcf = vcf.loc[vcf['SAMPLE'].isin(['1|0', '0|1'])]
		chr_pos = pd.concat([chr_pos, vcf[['#CHROM','POS']]], ignore_index=True, axis=0)
	if not os.path.exists(os.path.join(args.snps_output_path, 'snps.positions.txt')):
		chr_pos.to_csv(os.path.join(args.snps_output_path, 'snps.positions.txt'), sep='\t', header=False, index=False)
	vcf = chr_pos = None
	return

def find_intersecting_snps(args, bc):
	fis_script = '/usr/local/lib/python3.6/site-packages/wasp_map/WASP/mapping/find_intersecting_snps.py'
	subprocess.call(['python3', fis_script, os.path.join(args.bam_output_path, '{}_{}.possorted.bam'.format(args.sample_name, bc)), '--snp_tab', os.path.join(args.variants_dir, 'snp_tab.h5'), '--snp_index', os.path.join(args.variants_dir, 'snp_index.h5'), '--haplotype', os.path.join(args.variants_dir, 'haplotype.h5'), '--samples', args.sample_name, '--is_paired_end', '--is_sorted', '--output_dir', args.fis_output_path], stderr=subprocess.DEVNULL)
	return

def remap_reads(args, bc):
	fq1 = os.path.join(args.fis_output_path, '{}_{}.possorted.remap.fq1.gz'.format(args.sample_name, bc))
	fq2 = os.path.join(args.fis_output_path, '{}_{}.possorted.remap.fq2.gz'.format(args.sample_name, bc))
	if os.path.exists(fq1) and os.path.exists(fq2):
		bwa = subprocess.Popen(['bwa', 'mem', '-M', '-t', '1', args.genome_reference_file, fq1, fq2], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
		samtools_fixmate = subprocess.Popen(['samtools', 'fixmate', '-u', '--no-PG', '-r', '-', '-'], stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
		subprocess.call(['samtools', 'sort', '--no-PG', '-m', '2G', '-o', os.path.join(args.remap_output_path, '{}_{}.possorted.remapped.bam'.format(args.sample_name, bc)), '-'], stdin=samtools_fixmate.stdout)
	return

def filter_reads(args, bc):
	frr_script = '/usr/local/lib/python3.6/site-packages/wasp_map/WASP/mapping/filter_remapped_reads.py'
	toremap_bam = os.path.join(args.fis_output_path, '{}_{}.possorted.to.remap.bam'.format(args.sample_name, bc))
	remapped_bam = os.path.join(args.remap_output_path, '{}_{}.possorted.remapped.bam'.format(args.sample_name, bc))
	keep_bam = os.path.join(args.frr_output_path, '{}_{}.keep.bam'.format(args.sample_name, bc))
	if os.path.exists(toremap_bam) and os.path.exists(remapped_bam):
		subprocess.call(['python3', frr_script, toremap_bam, remapped_bam, keep_bam], stderr=subprocess.DEVNULL)
	return

def remove_duplicates(args, bc):
	rmdup_script = '/usr/local/lib/python3.6/site-packages/wasp_map/WASP/mapping/rmdup_pe.py'
	keep_bam1 = os.path.join(args.fis_output_path, '{}_{}.possorted.keep.bam'.format(args.sample_name, bc))
	keep_bam2 = os.path.join(args.frr_output_path, '{}_{}.keep.bam'.format(args.sample_name, bc))
	merged_bam = os.path.join(args.rmdup_output_path, '{}_{}.keep.merged.bam'.format(args.sample_name, bc))
	rmdup_bam = os.path.join(args.rmdup_output_path, '{}_{}.keep.merged.rmdup.bam'.format(args.sample_name, bc))

	if os.path.exists(keep_bam1) and os.path.exists(keep_bam2):
		merge = subprocess.Popen(['samtools', 'merge', '--no-PG', '-c', '-p', '-', keep_bam1, keep_bam2], stdout=subprocess.PIPE)
		subprocess.call(['samtools', 'sort', '--no-PG', '-m', '2G', '-o', merged_bam, '-'], stdin=merge.stdout)
		subprocess.call(['python3', rmdup_script, merged_bam, rmdup_bam]), stderr=subprocess.DEVNULL)
	return

def final_merge_pileup(args):
	rmdup_bams = glob.glob(os.path.join(args.rmdup_output_path, '{}_*.keep.merged.rmdup.bam'.format(args.sample_name)))
	final_bam = os.path.join(args.final_output_path, '{}.merged.rmdup.final.bam'.format(args.sample_name))
	merge = subprocess.Popen(['samtools', 'merge', '-u', '--no-PG', '-c', '-p', '-'] + rmdup_bams, stdout=subprocess.PIPE)
	subprocess.call(['samtools', 'sort', '-m', '96G', '--no-PG', '-o', final_bam, '-'], stdin=merge.stdout)
	snp_positions = os.path.join(args.snps_output_path, 'snps.positions.txt')
	pileup = os.path.join(args.pileup_output_path, '{}.pileup'.format(args.sample_name))
	subprocess.call(['samtools', 'mpileup', '-f', args.genome_reference_file, '-l', snp_positions, '--no-output-ins', '--no-output-ins', '--no-output-del', '--no-output-del', '--no-output-ends', '-o', pileup, final_bam])
	return

def cleanup_temp_files(args):
	return

def main(args):
	args.bam_output_path = os.path.join(args.output_prefix, args.sample_name, 'sc_bams')
	args.snps_output_path = os.path.join(args.output_prefix, args.sample_name, 'snps')
	args.fis_output_path = os.path.join(args.output_prefix, args.sample_name, 'find_snps')
	args.remap_output_path = os.path.join(args.output_prefix, args.sample_name, 'remap')
	args.frr_output_path = os.path.join(args.output_prefix, args.sample_name, 'filt')
	args.rmdup_output_path = os.path.join(args.output_prefix, args.sample_name, 'rmdup')
	args.final_output_path = os.path.join(args.output_prefix, args.sample_name, 'final')
	args.pileup_output_path = os.path.join(args.output_prefix, args.sample_name, 'pileup')
	for d in [args.bam_output_path, args.snps_output_path, args.fis_output_path, args.remap_output_path, 
			args.frr_output_path, args.rmdup_output_path, args.final_output_path, args.pileup_output_path]:
		try:
			os.makedirs(d)
		except:
			pass

	singlecell = pd.read_table(args.singlecell_file, sep='\t', header=0, index_col=0)
	singlecell = singlecell.loc[singlecell['sample']==args.sample_name]
	split_sc_bams(args, singlecell)
	create_snps_files(args)

	with Pool(processes=args.max_processes) as pool:
		pool.starmap(find_intersecting_snps, [(args, bc) for bc in singlecell.index])
		pool.starmap(remap_reads, [(args, bc) for bc in singlecell.index])
		pool.starmap(filter_reads, [(args, bc) for bc in singlecell.index])
		pool.starmap(remove_duplicates, [(args, bc) for bc in singlecell.index])
	final_merge_pileup(args)
	return
 
def process_args():
	parser = argparse.ArgumentParser(description='Split bam file and run WASP')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-b', '--bam-file', required=True, type=str, help='Path to possorted_bam.bam')
	io_group.add_argument('-s', '--singlecell-file', required=True, type=str, help='Path to barcodes file')
	io_group.add_argument('-v', '--variants-dir', required=True, type=str, help='Path to variants directory')
	io_group.add_argument('-o', '--output-prefix', required=True, type=str, help='Path to output directory')
	io_group.add_argument('--sample-name', required=True, type=str, help='Sample name')
	io_group.add_argument('--max-processes', required=False, type=int, default=8, help='Max number of processes to run')
	io_group.add_argument('--genome-reference-file', required=False, type=str, default='/home/joshchiou/references/male.hg19.fa', help='Path to genome reference file')
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
	args = process_args()
	main(args)
