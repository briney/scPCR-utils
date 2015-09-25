#!/usr/bin/python
# filename: containers.py

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


from colormath.color_conversions import convert_color
from colormath.color_objects import LabColor, sRGBColor
from colormath.color_diff import delta_e_cmc


class BaseContainer(object):
	"""Base class for Sample and Well objects"""
	def __init__(self, cell):
		super(BaseContainer, self).__init__()
		self.cell = cell
		self.color = cell.fill.bgColor.rgb


class Sample(BaseContainer):
	"""Container for Samples"""
	def __init__(self, cell):
		super(Sample, self).__init__(cell)
		self.name = cell.value

	@property
	def name(self):
		return self._name

	@name.setter
	def name(self, value):
		if not value:
			self._name = value
		else:
			try:
				self._name = str(int(value))
			except ValueError:
				self._name = str(value)


class Well(BaseContainer):
	"""Container for Wells"""
	def __init__(self, cell):
		super(Well, self).__init__(cell)
		self.value = cell.value
		self._well = None
		self._sample = None

	@property
	def well(self):
		return self._well

	def set_well(self, row, col):
		col = str(col)
		if len(col) == 1:
			col = '0' + col
		self._well = str(row) + col

	@property
	def sample(self):
	    return self._sample

	@sample.setter
	def sample(self, sample):
		self._sample = sample


class Plate(object):

	def __init__(self, name, wells, samples):
		self.name = name
		self.wells = wells
		self.samples = samples
		self.sample_names = [s.name for s in samples if self.name]
		self.sample_colors = {s.name: s.color for s in samples if s.name}
		self.sample_name_by_color = {s.color: s.name for s in samples}
		self.well_names = [w.well for w in wells]
		self.value = [w.value for w in wells]
		self.color = [w.color for w in wells]
		self._set_well_samples()

	@property
	def values(self):
		return self._values

	@values.setter
	def values(self, values):
		self._values = self._build_values_dict(values)

	@property
	def colors(self):
		return self._colors

	@colors.setter
	def colors(self, colors):
		self._colors = self._build_values_dict(colors)

	@property
	def sample(self):
		return self._sample

	@sample.setter
	def sample(self, samples):
		self._sample = samples

	def _set_well_samples(self):
		for w in self.wells:
			try:
				w.sample = self.sample_name_by_color[w.color]
			except KeyError:
				w.sample = self._infer_sample(w)

	def _infer_sample(self, well):
		# return None if the well is empty
		if well.color == '00000000' and well.color not in [s.color for s in self.samples]:
			return None
		# if the whole plate contains a single sample, assign that sample
		samples = [s for s in self.samples if s.name]
		if len(samples) == 1:
			return samples[0].name
		# find the sample color that's closest to the well color
		sample_labs = [convert_color(
			sRGBColor(
				*self.hex_to_rgb(s.color)), LabColor) for s in samples]
		well_lab = convert_color(
			sRGBColor(
				*self.hex_to_rgb(well.color)), LabColor)
		diffs = []
		for slab in sample_labs:
			diffs.append(delta_e_cmc(well_lab, slab))
		best_match = samples[diffs.index(min(diffs))]
		return best_match.name

	def _build_values_dict(self, values):
		return {w: v for w, v in zip(self.well_names, values)}

	@staticmethod
	def hex_to_rgb(value):
		value = value[2:]
		lv = len(value)
		return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
