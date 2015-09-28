#!/usr/bin/env python
# filename: scpcr_demultiplexing.py


#
# Copyright (c) 2015 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


from __future__ import print_function

import argparse
import glob
import logging
import os
import sqlite3
import string
import subprocess as sp
import sys
import tempfile
import time
import urllib

from pymongo import MongoClient

import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', dest='output', required=True,
					help="Output directory for demultiplexed FASTA files. Required.")
parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
					help="Directory for temporary files. Required.")
parser.add_argument('-l', '--log', dest='log',
					help="Location for the log file. \
					Default is <output>/<database>.log.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="Name of the MongoDB database containing un-demultiplexed sequences. Required.")
parser.add_argument('-c', '--collection', dest='collection', default=None,
					help="Name of the MongoDB collection to query. \
					If not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-u', '--user', dest='user', default=None,
					help="Username for the MongoDB server. Not used if not provided.")
parser.add_argument('-p', '--password', dest='password', default=None,
					help="Password for the MongoDB server. Not used if not provided.")
parser.add_argument('-I', '--index-file', dest='index_file', required=True,
					help="File containing the indexes, one per line. \
					File must have either 96 or 384 indexes.")
parser.add_argument('-M', '--plate-map', dest='plate_map', required=True,
					help="Plate map. Can either be provided as a list of well names \
					(one per line, in row order starting at A01) or as a list of wells \
					and well names (one per line, separated by white space) in any order. \
					If providing simply a list of well names and not every well is used, \
					leave blank lines for unused wells.")
parser.add_argument('--index-position', default='start', choices=['start', 'end'],
					help="Position of the indexes. Choices are 'start', if they're \
					at the start of the raw merged read (start of read 1) or 'end' if \
					they're at the end of the raw merged read (start of read 2).")
parser.add_argument('--index-length', default=0, type=int,
					help="Length of the index, in nucleotides. \
					Default is to parse the index length from the file of index sequences.")
parser.add_argument('--index-reverse-complement', default=False, action='store_true',
					help="Set if the indexes in the supplied index file are the reverse \
					complement of the indexes as they will appear in the sequences. \
					Default is False.")
parser.add_argument('--score-cutoff-heavy', default=200, type=int,
					help="V-gene alignment score cutoff for heavy chains. \
					Alignment score must be equal to or higher than cutoff to be considered for clustering. \
					Default is 200.")
parser.add_argument('--score-cutoff-light', default=100, type=int,
					help="V-gene alignment score cutoff for kappa/lambda chains. \
					Alignment score must be equal to or higher than cutoff to be considered for clustering. \
					Default is 100.")
parser.add_argument('--cdhit-threshold', default=0.96, type=float,
					help="Threshold for CD-HIT clustering. \
					Default is 0.96.")
parser.add_argument('--minimum-cluster-size', default=50, type=int,
					help="Minimum size of a CD-HIT cluster of sequences. \
					Centroids will not be determined for clusters below this cutoff. \
					Default is 50.")
parser.add_argument('--minimum-cluster-fraction', default=0.5, type=float,
					help="Minimum fraction for a CD-HIT cluster (relative to the total \
					number of sequences in a given well) for a centroid to be determined. \
					For example, if set to 0.7, the largest CD-HIT cluster for a well must comprise \
					70 percent of all sequences in that well. \
					Default is 0.5.")
parser.add_argument('--debug', dest='debug', action='store_true', default=False,
					help="If set, will run in debug mode.")
args = parser.parse_args()



###############
#   MongoDB
###############

if args.user and args.password:
	password = urllib.quote_plus(password)
	uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
	conn = MongoClient(uri)
else:
	conn = MongoClient(args.ip, 27017)
db = conn[args.db]


def get_collections():
	if args.collection:
		return [args.collection, ]
	collections = db.collection_names(include_system_collections=False)
	return sorted(collections)


def get_sequences(collection, chain):
	score_cutoff = args.score_cutoff_heavy if chain == 'heavy' else args.score_cutoff_light
	seqs = db[collection].find({'chain': chain, 'prod': 'yes', 'v_gene.score': {'$gte': score_cutoff}}, {'seq_id': 1, 'raw_query': 1, 'vdj_nt': 1})
	return [s for s in seqs]



######################
#   SQLite database
######################


def build_seq_db(seqs):
	sys.stdout.flush()
	db_path = os.path.join(args.temp_dir, 'seq_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	create_cmd = get_seq_db_creation_cmd()
	insert_cmd = get_seq_db_insert_cmd()
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute(create_cmd)
	c.executemany(insert_cmd, seqs)
	sys.stdout.flush()
	start = time.time()
	c.execute('CREATE INDEX seq_index ON seqs (seq_id)')
	return c


def get_seq_db_creation_cmd():
	return '''CREATE TABLE seqs (seq_id text, vdj_nt text)'''


def get_seq_db_insert_cmd():
	return 'INSERT INTO seqs VALUES (?,?)'


def remove_sqlite_db():
	db_path = os.path.join(args.temp_dir, 'seq_db')
	os.unlink(db_path)



##############
#   CD-HIT
##############


def cdhit_clustering(seqs, bin_id):
	print('clustering...')
	legit = False
	seq_db = build_seq_db(seqs)
	temp_dir = args.temp_dir
	infile = make_cdhit_input(seqs)
	outfile = os.path.join(temp_dir, 'clust')
	logfile = open(os.path.join(temp_dir, 'log'), 'a')
	threshold = 0.9
	do_cdhit(infile.name, outfile, logfile)
	clust_handle = open('{}.clstr'.format(outfile), 'r')
	seq, size, total_count = parse_clusters(clust_handle, seq_db)
	os.unlink(infile.name)
	os.unlink(os.path.join(args.temp_dir, 'log'))
	os.unlink(outfile)
	os.unlink(outfile + '.clstr')
	remove_sqlite_db()
	if size >= args.minimum_cluster_size and 1. * size / len(seqs) >= args.minimum_cluster_fraction:
		legit = True
	if legit:
		print('PASS')
		logging.info('{}: PASSED, {} total sequences, {} in largest cluster'.format(
			bin_id, total_count, size))
		return seq
	print('FAIL')
	logging.info('{}: FAILED, {} total sequences, {} in largest cluster'.format(
			bin_id, total_count, size))
	return None


def make_cdhit_input(seqs):
	fastas = ['>{}\n{}'.format(s[0], s[1]) for s in seqs]
	infile = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
	infile.write('\n'.join(fastas))
	infile.close()
	return infile


def do_cdhit(fasta, clust, log):
	sys.stdout.flush()
	start_time = time.time()
	cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(fasta, clust, args.cdhit_threshold)
	cluster = sp.Popen(cdhit_cmd, shell=True, stdout=log)
	cluster.communicate()


def parse_centroids(centroid_handle, sizes=None):
	counter = 0
	centroids = []
	for seq in SeqIO.parse(centroid_handle, 'fasta'):
		if sizes:
			size = sizes[counter]
			centroids.append('>{}_{}\n{}'.format(seq.id, size, str(seq.seq)))
		else:
			centroids.append('>{}\n{}'.format(seq.id, str(seq.seq)))
		counter += 1
	return centroids


def parse_clusters(cluster_handle, seq_db):
	sys.stdout.flush()
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	sys.stdout.flush()
	start = time.time()
	cluster_lengths = []
	for cluster in clusters:
		length = len(cluster) - 1
		cluster_id = get_cluster_id(cluster)
		cluster_seq = get_cluster_seq(cluster_id, seq_db)
		cluster_lengths.append((length, (cluster_id, cluster_seq)))
	cluster_lengths.sort(key=lambda x: x[0], reverse=True)
	total_seqs = sum([l[0] for l in cluster_lengths])
	biggest_cluster = cluster_lengths[0][1]
	print('{} sequences were clustered.\nThe largest cluster contained {} sequences'.format(total_seqs, cluster_lengths[0][0]))
	return biggest_cluster[1], cluster_lengths[0][0], total_seqs


def parse_cluster_sizes(cluster_handle):
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	lengths = []
	for cluster in clusters:
		lengths.append(len(cluster) - 1)
	return lengths


def get_cluster_id(cluster):
	ids = []
	for c in cluster[1:]:
		if c:
			ids.append(c.split()[2][1:-3])
	return ids[0]


def get_cluster_seq(seq_id, seq_db):
	seq_db.execute('''SELECT seqs.seq_id, seqs.vdj_nt
							FROM seqs
							WHERE seqs.seq_id LIKE "{}"'''.format(seq_id))
	return seq_db.fetchone()



###############
#   Indexes
###############


def parse_indexes():
	indexes = {}
	index_seqs = []
	# parse index sequences from file
	with open(args.index_file) as f:
		for line in f:
			index_seqs.append(line.strip())
	# define row and column ranges
	if len(index_seqs) == 384:
		rows = [c for c in string.ascii_uppercase[:16]]
		columns = [str(c) if len(str(c)) == 2 else '0{}'.format(c) for c in range(1, 25)]
	elif len(index_seqs) == 96:
		rows = [c for c in string.ascii_uppercase[:8]]
		columns = [str(c) if len(str(c)) == 2 else '0{}'.format(c) for c in range(1, 13)]
	else:
		print('Index file must contain either 96 or 384 index sequences.')
		sys.exit(1)
	# perform a couple of sanity checks on the index sequences
	index_lengths = list(set([len(i) for i in index_seqs]))
	if len(index_lengths) > 1 and not args.index_length:
		print('Indexes must all be the same length, or --index-length must be provided.')
		sys.exit(1)
	if not args.index_length:
		args.index_length = index_lengths[0]
	if min(index_lengths) < args.index_length:
		print('All indexes must be at least as long as --index-length.')
		sys.exit(1)
	# build a dictionary of well locations and index sequences
	wells = []
	for row in rows:
		for column in columns:
			wells.append('{}{}'.format(row, column))
	for well, index in zip(wells, index_seqs):
		indexes[index] = well
	return indexes


def parse_plate_map(wells):
	well_names = []
	well_map = {}
	with open(args.plate_map) as f:
		for line in f:
			well_names.append(line.strip().split())
	if len(well_names[0]) == 1:
		well_names = [(w, n[0]) for w, n in zip(wells, well_names) if n]
	for (w, n) in well_names:
		well_map[w] = n
	return well_map


def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(seq)
	rc = reversed([complement.get(b, b) for b in bases])
	return ''.join(rc)


def bin_by_index(sequences, indexes):
	bins = {w: [] for w in indexes.values()}
	pos = args.index_length if args.index_position == 'start' else -1 * args.index_length
	for seq in sequences:
		if args.index_position == 'start':
			index = seq['raw_query'][:pos]
		else:
			index = seq['raw_query'][pos:]
		if args.index_reverse_complement:
			index = reverse_complement(index)
		if index in indexes:
			well = indexes[index]
			bins[well].append((seq['seq_id'], seq['vdj_nt']))
	return bins




##########################
#   Logging and Output
##########################


def setup_logging():
	logfile = args.log if args.log else os.path.join(args.output, '{}.log'.format(args.db))
	if args.debug:
		logging.basicConfig(filename=logfile,
							filemode='w',
							format='[%(levelname)s] %(asctime)s %(message)s',
							level=logging.DEBUG)
	else:
		logging.basicConfig(filename=logfile,
							filemode='w',
							format='[%(levelname)s] %(asctime)s %(message)s',
							level=logging.INFO)
	logging.info('LOG LOCATION: {}'.format(logfile))
	log_options()


def log_options():
	logging.info('OUTPUT DIRECTORY: {}'.format(args.output))
	logging.info('DATABASE: {}'.format(args.db))
	logging.info('INDEX FILE: {}'.format(args.index_file))
	logging.info('PLATE-MAP FILE: {}'.format(args.plate_map))
	logging.info('HEAVY CHAIN SCORE THRESHOLD: {}'.format(args.score_cutoff_heavy))
	logging.info('LIGHT CHAIN SCORE THRESHOLD: {}'.format(args.score_cutoff_light))
	logging.info('USER-SUPPLIED INDEX LENGTH: {}'.format(
		args.index_length if args.index_length > 0 else 'None'))
	logging.info('INDEX POSITION: {}'.format(args.index_position))
	logging.info('INDEX REVERSE COMPLEMENT: {}'.format(args.index_reverse_complement))
	logging.info('CD-HIT THRESHOLD: {}'.format(args.cdhit_threshold))
	logging.info('MINIMUM CLUSTER SIZE: {}'.format(args.minimum_cluster_size))
	logging.info('MINIMUM CLUSTER FRACTION: {}'.format(args.minimum_cluster_fraction))


def log_output(bins, seqs):
	num_bins = len([b for b in bins.keys() if len(bins[b]) >= args.minimum_cluster_size])
	num_seqs = len(seqs)
	logging.info('RESULTS: Of {} bins with at least {} reads, {} passed filter'.format(
		num_bins, args.minimum_cluster_size, num_seqs))


def write_output(seqs, outname):
	outpath = os.path.join(args.output, outname)
	seq_string = '\n'.join(['>{}\n{}'.format(s[0], s[1][1]) for s in seqs])
	open(outpath, 'w').write(seq_string)




################
#   Printing
################


def print_plate_info(name, collection):
	name_string = '     Processing {} (collection {})     '.format(name, collection)
	print('\n\n')
	print('=' * len(name_string))
	print(name_string)
	print('=' * len(name_string))
	logging.info('')
	logging.info('COLLECTION: {}'.format(collection))
	logging.info('PLATE NAME: {}'.format(name))


def print_bin_info(b):
	bin_string = '  {}  '.format(b)
	print('\n')
	print(bin_string)
	print('-' * len(bin_string))




def main():
	indexes = parse_indexes()
	plate_map = parse_plate_map(sorted(indexes.values()))
	for collection in get_collections():
		if collection not in plate_map:
			print('\n\n{} was not found in the supplied plate map file.'.format(
				collection))
			continue
		plate_name = plate_map[collection]
		print_plate_info(plate_name, collection)
		plate_seqs = []
		for chain in ['heavy', 'kappa', 'lambda']:
			print('\n\nQuerying for {} chain sequences'.format(chain))
			logging.info('{} CHAIN'.format(chain.upper()))
			sequences = get_sequences(collection, chain)
			print('Retrieved {} sequences.\n'.format(len(sequences)))
			logging.info('QUERY RESULTS: {} {} chain sequences met the quality threshold'.format(
				len(sequences), chain.lower()))
			bins = bin_by_index(sequences, indexes)
			if max([len(b) for b in bins.values()]) < args.minimum_cluster_size:
				logging.info('The entire plate had fewer than {} sequences and was not processed'.format(args.minimum_cluster_size))
				continue
			for b in sorted(bins.keys()):
				if len(bins[b]) < 25:
					continue
				print_bin_info(b)
				centroid = cdhit_clustering(bins[b], b)
				if centroid:
					plate_seqs.append((b, centroid))
			log_output(bins, plate_seqs)
		write_output(plate_seqs, plate_name)
	print('\n')


if __name__ == '__main__':
	setup_logging()
	main()
