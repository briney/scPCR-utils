#!/usr/bin/python
# filename: excel_platemap_parser.py

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
import itertools
import os
import string

from openpyxl import load_workbook

from containers import Well, Sample, Plate


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='input', required=True,
					help="The input file, in Excel format. Required.")
parser.add_argument('-o', '--out', dest='output', required=True,
					help="The output directory, into which the map files will be deposited. \
					If the directory does not exist, it will be created. \
					Required.")
parser.add_argument('-e', '--experiment', dest='experiment', default=None,
					help="Name of the experiment. \
					If not supplied, the input filename (after removing the .xlsx extension) \
					will be used.")
parser.add_argument('--plate-prefix', default='plate',
					help="Prefix for plate filenames, if names aren't to be parsed from the platemap file. \
					Default is 'plate'.")
parser.add_argument('--plate-numbering-start', default=1, type=int,
					help="Number at which to start the plate numbering, if not parsing plate names. \
					Default is 1.")
parser.add_argument('--plate-delim-prefix', required=True,
					help='Text at the start of a well marking the boundary between two plates. \
					Often, the best choice is the well with the plate name (if part of the platemap) \
					or the sample row. Required.')
parser.add_argument('--plate-delim-well-number', required=True, type=int,
					help="The well number (using 1-based indexing) of the plate delimiter well. \
					Used in conjunction with --plate-delim-row-number. \
					Required.")
parser.add_argument('--plate-name-row-number', default=None, type=int,
					help="The number of rows between the plate delimiter row and the row containing \
					the plate name. Only required if parsing plate names. \
					Default is None, which results in plate names not being parsed.")
parser.add_argument('--plate-name-well-number', default=None, type=int,
					help="The well number (using 1-based indexing) of the plate name well. \
					Used in conjunction with --plate-name-row-number. \
					Default is None, which results in plate names not being parsed.")
parser.add_argument('--sample-row-number', required=True, type=int,
					help="The number of rows between the plate delimiter row and the row containing \
					the sample names. Required.")
parser.add_argument('--sample-well-number', required=True, type=int,
					help="The well number (using 1-based indexing) of the first sample name. \
					Used in conjunction with --sample-row-number. Required.")
parser.add_argument('--sample-well-offset', default=1, type=int,
					help="The number of wells each sample name occupies. \
					Default (1) works in most circumstances, unless wells containing \
					sample names are merged. If samples are in a merged 'well' that originally was \
					two wells, use 2. For a merged three-well, use 3.")
parser.add_argument('--max-sample-number', default=4, type=int,
					help="Maximum number of samples in a single plate. Default is 4. \
					Increase if there are more than four samples for any well in the platemap.")
args = parser.parse_args()


def get_experiment():
	if args.experiment:
		return args.experiment
	exp = os.path.basename(args.input).rstrip('.xlsx').rstrip('.xls')
	return exp


def get_plate_blocks(ws):
	num = args.plate_delim_well_number
	prefix = args.plate_delim_prefix
	plate_blocks = []
	curr_plate = []
	for row in ws.rows:
		if not row[num - 1].value:
			curr_plate.append(row)
			continue
		if row[num - 1].value.startswith(prefix):
			plate_blocks.append(curr_plate)
			curr_plate = [row, ]
		else:
			curr_plate.append(row)
	if curr_plate:
		plate_blocks.append(curr_plate)
	return plate_blocks


def parse_plates(raw_plates):
	# plates = []
	# for i, rp in enumerate(raw_plates):
	# 	plate_name = 'plate' + str(i + 1) if len(str(i + 1)) == 2 else 'plate0' + str(i + 1)
	# 	wells = []
	# 	for row in rp:
	# 		cells = [cell for cell in row]
	# 		cell_vals = [cell.value for cell in cells]
	# 		if 'Sample name(s)' in cell_vals:
	# 			start = cell_vals.index('Sample name(s)')
	# 			samples = [Sample(cell) for cell in cells[start + 2:start + 9:2]]
	# 		if not row[0].value:
	# 			continue
	# 		if str(row[0].value) not in string.ascii_uppercase[:8]:
	# 			continue
	# 		row_wells = [Well(cell) for cell in row[1:13]]
	# 		for c, rw in enumerate(row_wells):
	# 			r = str(row[0].value)
	# 			rw.set_well(r, c + 1)
	# 		wells += row_wells
	# 	plates.append(Plate(plate_name, wells, samples))
	# return plates

	plates = []
	for i, rp in enumerate(raw_plates):
		if args.plate_name_row_number is not None and args.plate_name_well_number is not None:
			row_num = args.plate_name_row_number
			well_num = args.plate_name_well_number - 1
			plate_name = rp[row_num][well_num].value
		else:
			num = str(i + args.plate_numbering_start)
			if len(num) < 2:
				num = '0' + num
			plate_name = '{}{}'.format(args.plate_prefix, num)
		print('Processing plate: {}'.format(plate_name))
		samples = parse_samples(rp)
		wells = parse_plate_grid(rp)
		plates.append(Plate(plate_name, wells, samples))
	return plates


def parse_samples(plate):
	samples = []
	row = plate[args.sample_row_number]
	start = args.sample_well_number - 1
	end = start + args.max_sample_number
	step = args.sample_well_offset
	for cell in row[start:end:step]:
		samples.append(Sample(cell))
	return samples


def parse_plate_grid(raw_plate):
	# trim empty rows from the end of the raw plate
	rp = [list(row) for row in raw_plate]
	rp_set = list(set([c.value for c in rp[-1]]))
	while len(rp_set) == 1 and rp_set[0] is None:
		rp = rp[:-1]
		rp_set = list(set([c.value for c in rp[-1]]))
	# find the column containing row labels
	labelc = 0
	for i in range(len(rp)):
		column = [r[i] for r in rp]
		cvals = [c.value for c in column if c.value is not None]
		if len([c for c in cvals if str(c) in string.ascii_uppercase[:8]]) == 8:
			labelc = i
	gridc = labelc + 1
	# find the row that contains the column labels
	labelr = 0
	for i, row in enumerate(raw_plate):
		if len([r for r in row if r.value in range(1, 13)]) == 12:
			labelr = i
	gridr = labelr + 1
	# parse the cells from the plate grid
	plate = []
	for row in raw_plate[gridr:gridr + 8]:
		rname = row[labelc].value
		for cname, cell in enumerate(row[gridc:gridc + 12]):
			well = Well(cell)
			well.set_well(str(rname), str(cname + 1))
			plate.append(well)
	return plate




def write_output(plates, experiment):
	for i, plate in enumerate(plates):
		# num = str(i + args.plate_numbering_start)
		# if len(num) < 2:
		# 	num = '0' + num
		# ohandle = open(os.path.join(args.output, '{}{}'.format(args.plate_prefix, num)), 'w')
		ohandle = open(os.path.join(args.output, plate.name), 'w')
		output = []
		for well in plate.wells:
			output.append('{}\t{}\t{}\t{}'.format(
				well.well, well.sample, well.value, experiment))
		ohandle.write('\n'.join(output))


def main():
	experiment = get_experiment()
	wb = load_workbook(args.input)
	ws = wb[wb.get_sheet_names()[0]]
	plate_blocks = get_plate_blocks(ws)
	plural = '' if len(plate_blocks) <= 2 else 's'
	print('\nFound {} plate{} in the input file'.format(len(plate_blocks) - 1, plural))
	print('Experiment name: {}\n'.format(experiment))
	plates = parse_plates(plate_blocks[1:])
	write_output(plates, experiment)
	print()



if __name__ == '__main__':
	main()

