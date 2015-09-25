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
parser.add_argument('--plate-prefix', default='plate',
					help="Prefix for plate filenames.\
					Default is 'plate'.")
parser.add_argument('--plate-start', default=1, type=int,
					help="Number at which to start the plate numbering. \
					Default is 1.")
args = parser.parse_args()


def get_plate_blocks(ws):
	return [list(g) for k, g in itertools.groupby(
		ws.rows, lambda x: x[0].value == 'BARCODE:') if not k]


def parse_plates(raw_plates):
	plates = []
	for i, rp in enumerate(raw_plates):
		plate_name = 'plate' + str(i + 1) if len(str(i + 1)) == 2 else 'plate0' + str(i + 1)
		wells = []
		for row in rp:
			cells = [cell for cell in row]
			cell_vals = [cell.value for cell in cells]
			if 'Sample name(s)' in cell_vals:
				start = cell_vals.index('Sample name(s)')
				samples = [Sample(cell) for cell in cells[start + 2:start + 9:2]]
			if not row[0].value:
				continue
			if str(row[0].value) not in string.ascii_uppercase[:8]:
				continue
			row_wells = [Well(cell) for cell in row[1:13]]
			for c, rw in enumerate(row_wells):
				r = str(row[0].value)
				rw.set_well(r, c + 1)
			wells += row_wells
		plates.append(Plate(plate_name, wells, samples))
	return plates


def write_output(plates):
	for i, plate in enumerate(plates):
		num = str(i + args.plate_start)
		if len(num) < 2:
			num = '0' + num
		ohandle = open(os.path.join(args.output, '{}{}'.format(args.plate_prefix, num)), 'w')
		output = []
		for well in plate.wells:
			output.append('{}\t{}\t{}'.format(well.well, well.sample, well.value))
		ohandle.write('\n'.join(output))


def main():
	wb = load_workbook(args.input)
	ws = wb[wb.get_sheet_names()[0]]
	plate_blocks = get_plate_blocks(ws)
	plates = parse_plates(plate_blocks[1:])
	write_output(plates)



if __name__ == '__main__':
	main()

