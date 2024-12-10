#!/usr/bin/env python3
"""
Convert structural variants in a VCF to CGH (CytoSure) format
"""

import argparse
import pandas as pd
import logging
import gzip
import numpy as np
import os
from collections import namedtuple, defaultdict
from io import StringIO
from lxml import etree
from cyvcf2 import VCF

from .constants import *

from .__version__ import __version__

logger = logging.getLogger(__name__)

Event = namedtuple('Event', ['chrom', 'start', 'end', 'type', 'info'])

def wisecondorx_events(args, CONTIG_LENGTHS):

	variants  = pd.read_csv(args.wisecondorx_aberrations, sep="\t", header = 0)
	skipped = 0
	
	for _, variant in variants.iterrows():
		chrom = str(variant["chr"])
		if (chrom) not in CONTIG_LENGTHS:
			continue

		start = int(variant["start"])
		end = int(variant["end"])
		sv_type =str(variant["type"])
		
		if sv_type=='loss':
			sv_type='DEL'
		elif sv_type=='gain':
			sv_type='DUP'

		if args.wcx_size:
			variant_size = (end-start)
			if variant_size <= args.wcx_size:
				#Too short
				skipped += 1
				logger.debug(f'skipped variant start={start},end={end},type= {sv_type} for being too short. Size of skipped variant={variant_size}bp')
				logger.info(f'skipped {skipped} variants for being too short') 
				continue		
				
		logger.debug('%s at %s:%s-%s (%s bp), Total skipped due to size: %s', sv_type, chrom, start+1, end, end - start, skipped)
		yield Event(chrom=chrom, start=start, end=end, type=sv_type, info={})


def strip_template(path):
	"""
	Read in the template CGH file and strip it of everything that we don't need.

	Return the lxml.etree object.
	"""
	tree = etree.parse(path)

	# Remove all aberrations
	parent = tree.xpath('/data/cgh/submission')[0]
	for aberration in parent.xpath('aberration'):
		parent.remove(aberration)

	# Remove all except the first probe (in the order in which they occur in
	# the file) on each chromosome. Chromosomes without probes are not
	# clickable in the CytoSure UI.
	parent = tree.xpath('/data/cgh/probes')[0]
	seen = set()
	for probe in parent:
		chrom = probe.attrib.get('chromosome')
		if not chrom or chrom in seen:
			parent.remove(probe)
		else:
			seen.add(chrom)

	# Remove all segments
	parent = tree.xpath('/data/cgh/segmentation')[0]
	for segment in parent:
		parent.remove(segment)

	return tree


def make_probe(parent, chromosome, start, end, height, text, original_coverage=None):
	probe = etree.SubElement(parent, 'probe')
	probe.attrib.update({
		'name': text,
		'chromosome': CHROM_RENAME.get(chromosome, chromosome),
		'start': str(start + 1),
		'stop': str(end),
		'normalized': '{:.3f}'.format(-height),
		'smoothed': '0.0',
		'smoothed_normalized': '0.0',
		'sequence': 'AACCGGTT',
	})

	if original_coverage is not None:
		probe.attrib['original_coverage'] = '{:.3f}'.format(original_coverage)

	red = 1000
	green = red * 2**height

	spot = etree.SubElement(probe, 'spot')
	spot.attrib.update({
		'index': '1',
		'row': '1',
		'column': '1',
		'red': str(red),
		'green': '{:.3f}'.format(green),
		'gSNR': '100.0',
		'rSNR': '100.0',
		'outlier': 'false',
	})
	return probe


def make_segment(parent, chromosome, start, end, height):
	segment = etree.SubElement(parent, 'segment')
	segment.attrib.update({
		'chrId': CHROM_RENAME.get(chromosome, chromosome),
		'numProbes': '100',
		'start': str(start + 1),
		'stop': str(end),
		'average': '{:.3f}'.format(-height),  # CytoSure inverts the sign
	})
	return segment


def make_aberration(parent, chromosome, start, end, comment=None, method='converted from VCF',
		confirmation=None, n_probes=0, copy_number=99):
	"""
	comment -- string
	method -- short string
	confirmation -- string
	"""
	# Set gain to false for dels
	is_gain = 'true'
	if (confirmation == "DEL"):
		is_gain = 'false'

	aberration = etree.SubElement(parent, 'aberration')
	aberration.attrib.update(dict(
		chr=CHROM_RENAME.get(chromosome, chromosome),
		start=str(start + 1),
		stop=str(end),
		maxStart=str(start + 1),
		maxStop=str(end),
		copyNumber=str(copy_number),
		initialClassification='Unclassified',
		finalClassification='Unclassified',
		inheritance='Not_tested',
		numProbes=str(n_probes),
		startProbe='',
		stopProbe='',
		maxStartProbe='',
		maxStopProbe='',
		gain=is_gain,
		method=method,
		# TODO fill in the following values with something sensible
		automationLevel='1.0',
		baseline='0.0',
		mosaicism='0.0',
		inheritanceCoverage='0.0',
		logRatio='-0.4444',  # mean log ratio
		p='0.003333',  # p-value
		sd='0.2222',  # standard deviation
	))
	if comment:
		e = etree.SubElement(aberration, 'comments')
		e.text = comment
	if confirmation:
		e = etree.SubElement(aberration, 'confirmation')
		e.text = confirmation
	return aberration


def spaced_probes(start, end, probe_spacing=PROBE_SPACING):
	"""
	Yield nicely spaced positions along the interval (start, end).
	- start and end are always included
	- at least three positions are included
	"""
	l = end - start
	n = l // probe_spacing
	spacing = l / max(n, 2)  # float division
	i = 0
	pos = start
	while pos <= end:
		yield pos
		i += 1
		pos = start + int(i * spacing)


def probe_point(center, height=2.5, width=5001, steps=15):
	"""
	Yield (pos, height) pairs that "draw" a triangular shape (pointing upwards)
	"""
	pos_step = (width - 1) // (steps - 1)
	height_step = height / ((steps - 1) // 2)
	for i in range(-(steps // 2), steps // 2 + 1):
		yield center + i * pos_step, height - height_step * abs(i) + 0.1


def format_comment(info):
	comment = ''
	for k, v in sorted(info.items()):
		if k in ('CSQ', 'SVTYPE'):
			continue
		comment += '\n{}: {}'.format(k, v)
	return comment


def merge_intervals(intervals):
	"""Merge overlapping intervals into a single one"""
	events = [(coord[0], 'START') for coord in intervals]
	events.extend((coord[1], 'STOP') for coord in intervals)
	events.sort()
	active = 0
	start = 0
	for pos, what in events:
		# Note adjacent 'touching' events are merged because 'START' < 'STOP'
		if what == 'START':
			if active == 0:
				start = pos
			active += 1
		else:
			active -= 1
			if active == 0:
				yield (start, pos)


def add_probes_between_events(probes, chr_intervals, CONTIG_LENGTHS):
	for chrom, intervals in chr_intervals.items():
		if chrom not in CONTIG_LENGTHS:
			continue
		intervals = merge_intervals(intervals)
		for start, end in complement_intervals(intervals, CONTIG_LENGTHS[chrom]):
			for pos in spaced_probes(start, end, probe_spacing=200000):
				# CytoSure does not display probes at height=0.0
				make_probe(probes, chrom, pos, pos + 60, 0.01, 'between events')


class CoverageRecord:
	__slots__ = ('chrom', 'start', 'end', 'coverage')

	def __init__(self, chrom, start, end, coverage):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.coverage = coverage

def parse_wisecondorx_coverages(args):
	probe_data=[]
	opener=open
	first=True
	df = pd.read_csv(args.wisecondorx_cov, sep="\t", header=0)
	
	for _, row in df.iterrows():

		chrom = row["chr"]
		start = int(row["start"])
		end = int(row["end"])
		coverage =float(row["ratio"])
		
		if pd.isna(coverage):
			continue	

		yield CoverageRecord(chrom, start, end, coverage)


def group_by_chromosome(records):
	"""
	Group records by their .chrom attribute.

	Yield pairs (chromosome, list_of_records) where list_of_records
	are the consecutive records sharing the same chromosome.
	"""
	prev_chrom = None
	chromosome_records = []
	for record in records:
		if record.chrom != prev_chrom:
			if chromosome_records:
				yield prev_chrom, chromosome_records
				chromosome_records = []
		chromosome_records.append(record)
		prev_chrom = record.chrom
	if chromosome_records:
		yield prev_chrom, chromosome_records


def bin_coverages(coverages, n):
	"""
	Reduce the number of coverage records by re-binning
	each *n* coverage values into a new single bin.

	The coverages are assumed to be from a single chromosome.
	"""
	chrom = coverages[0].chrom
	for i in range(0, len(coverages), n):
		records = coverages[i:i+n]
		cov = sum(r.coverage for r in records) / len(records)
		yield CoverageRecord(chrom, records[0].start, records[-1].end, cov)



def add_coverage_probes(probes, args, CONTIG_LENGTHS):
	"""
	probes -- <probes> element
	path -- path to tab-separated file with coverages
	"""
	coverages = [r for r in parse_wisecondorx_coverages(args) if r.chrom in CONTIG_LENGTHS]

	non_zero_len = len([r for r in coverages if r.coverage != 0])
	mean_coverage = sum(r.coverage for r in coverages) / non_zero_len
	logger.info('Mean coverage excluding 0 values is %.2f', mean_coverage)

	n = 0
	for chromosome, records in group_by_chromosome(coverages):
		coverage_factor = 1

		iterable_records = bin_coverages(records,args.bins)
				
		for record in iterable_records:	
			#print("#"+record.chrom,record.start)

			height = record.coverage
			adjusted_height=(((2**height)-1)* 100) / 10

			adjusted_height = min(MAX_HEIGHT, adjusted_height)
			adjusted_height = max(MIN_HEIGHT, adjusted_height)
			
			if adjusted_height == 0.0:
				adjusted_height = 0.01			
				
#			make_probe(probes, record.chrom, record.start, record.end, adjusted_height, 'coverage')
			make_probe(probes, record.chrom, record.start, record.end, adjusted_height, 'coverage', record.coverage)

			n += 1
	logger.info('Added %s coverage probes', n)

#retrieve the sample id, assuming single sample vcf
def retrieve_sample_id(input, input_path):
	sample = os.path.basename(input_path).split(".")[0]
    
	return sample

def main():
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = argparse.ArgumentParser("VCF2cytosure - convert SV vcf files to cytosure")

	group = parser.add_argument_group('Input')
	group.add_argument('--genome',required=False, default=37, help='Human genome version. Use 37 for GRCh37/hg19, 38 for GRCh38 template.')
	group.add_argument('--wisecondorx_aberrations', type=str, required=False,help='path to aberrations.bed file from WisecondorX')
	group.add_argument('--bins',type=int,default=20,help='the number of coverage bins per probes default=20')
	group.add_argument('--wisecondorx_cov', type=str, help='path to bins.bed file')
	group.add_argument('--out',help='output file (default = the prefix of the input vcf)')
	group.add_argument('--wcx_size',type=int,help='Variants smaller than this size will be filtered out')

	group.add_argument('-V','--version',action='version',version="%(prog)s "+__version__ ,
			   help='Print program version and exit.')
	# parser.add_argument('xml', help='CytoSure design file')
	args= parser.parse_args()

	logger.info('vcf2cytosure %s', __version__)

	if not args.wisecondorx_aberrations:
		print("Provide variant file. --wisecondorx_aberrations. See -help")
		quit()	

	if int(args.genome) == 38:
		CGH_TEMPLATE = CGH_TEMPLATE_38
		CONTIG_LENGTHS = CONTIG_LENGTHS_38
	else:
		CGH_TEMPLATE = CGH_TEMPLATE_37
		CONTIG_LENGTHS = CONTIG_LENGTHS_37

	if not args.out:
		args.out=".".join(args.vcf.split(".")[0:len(args.vcf.split("."))-1])+".cgh"
	parser = etree.XMLParser(remove_blank_text=True)

	sex_male = "false"
	promega_sex = 'Female'

	sample_id=retrieve_sample_id(None, args.wisecondorx_aberrations)

	tree = etree.parse(StringIO(CGH_TEMPLATE.format(sample_id,sample_id,sample_id,sample_id,sex_male,promega_sex,sex_male)), parser)

	segmentation = tree.xpath('/data/cgh/segmentation')[0]
	probes = tree.xpath('/data/cgh/probes')[0]
	submission = tree.xpath('/data/cgh/submission')[0]
		

	chr_intervals = defaultdict(list)
	n = 0

	event_generator = wisecondorx_events(args, CONTIG_LENGTHS)

	for event in event_generator:
		end = event.end
		height = ABERRATION_HEIGHTS[event.type]
		make_segment(segmentation, event.chrom, event.start, end, height)

		comment = format_comment(event.info)

		make_aberration(submission, event.chrom, event.start, end, confirmation=event.type,
			comment=comment)		

		chr_intervals[event.chrom].append((event.start, event.end))
		# show probes at slightly different height than segments
		for pos in spaced_probes(event.start, event.end - 1):
			make_probe(probes, event.chrom, pos, pos + 60, height, event.type)
		n += 1
	if args.wisecondorx_cov:
		add_coverage_probes(probes, args, CONTIG_LENGTHS)
	else:
		add_probes_between_events(probes, chr_intervals, CONTIG_LENGTHS)

	tree.write(args.out, pretty_print=True)
	logger.info('Wrote %d variants to CGH', n)


if __name__ == '__main__':
	main()
