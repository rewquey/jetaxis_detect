#!/usr/bin/env python
# -*- encoding: utf-8

import netCDF4 as nc
import numpy as np

from jet_detect import jet_detect


# Script configuration
basepath = './in'
tvar, latvar, lonvar = 'time', 'latitude', 'longitude'
ufile, uvar = basepath + '/ei.ans.1979.avg_925_700.u.nc', 'u'
vfile, vvar = basepath + '/ei.ans.1979.avg_925_700.v.nc', 'v'
ofile = './out/1979_jet_axes.txt'

# Detection configuration
#jet_detect.jetint_thres = 0.55e-8	# K-threshold for instantaneous ERA-Interim, PV2-level, T84-resolution
jet_detect.jetint_thres = 0.124e-8	# K-threshold for instantaneous ERA-Interim, mean winds 925-700 hPa, T84-resolution
jet_detect.searchrad = 1.5		# Maximum distance of points along the jet axis in grid point indices
jet_detect.minlen = 2.0e6		# Minimum lenght of the jet axes in meters
jet_detect.grid_cyclic_ew = True	# Is grid periodic in x-direction?

MAX_POINTS = 10000
MAX_LINES = 500


def save_jetaxes_txt(ofile, time, ja, jaoff):
	f = open(ofile, 'w')

	f.write('# Detected jet axes positions\n')
	f.write('#\n')
	f.write('# Columns\n')
	f.write('# yidx\txidx\twind speed\n')
	f.write('\n')

	for tidx, time in zip(range(len(time)), time):
		startpoint = 0
		jetid = 1
		for endpoint in jaoff[tidx][1:]:
			endpoint = int(endpoint)
			if endpoint <= startpoint:
				break

			f.write('time=%d, jetid=%d\n' % (time, jetid))
			for ptinfo in ja[tidx,startpoint:endpoint]:
				f.write('%4.2f\t%4.2f\t%4.2f\n' % tuple(ptinfo))
			f.write('\n')

			jetid += 1
			startpoint = endpoint
	
	return


# Execute the detection only when script is run rather than imported
if __name__ == '__main__':
	# Load u-data and grid information
	fu = nc.Dataset(ufile)
	u = fu.variables[uvar][::].squeeze()
	time = fu.variables[tvar][::]
	lat = fu.variables[latvar][::]
	lon = fu.variables[lonvar][::]
	fu.close()

	ny, nx = len(lat), len(lon)
	
	# Load v-data
	fv = nc.Dataset(vfile)
	v = fv.variables[vvar][::].squeeze()
	fv.close()
	
	# Construct grid distances
	dx = np.ones((ny,nx)) * 111.111e3 * np.cos(np.pi/180.0 * lat)[:,np.newaxis]
	dy = np.ones((ny,nx)) * 111.111e3

	# Detect jets
	ja, jaoff = jet_detect.run_jet_detect(MAX_POINTS, MAX_LINES, u, v, dx, dy)

	# Save
	save_jetaxes_txt(ofile, time, ja, jaoff)


# the end
