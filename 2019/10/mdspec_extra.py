'''
Additional scripts to analyze MDS spectra

Author: Jake Nichols
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mdspec_py

import os.path
import random
import math
import csv
import pandas as pd
from scipy import optimize, signal
from lmfit import models
from collections import defaultdict

# first-cut analysis functions (now deprecated)
def cv_wavelength_plot(shotnumber=178350):
	"""
	Plot of spectrum at various times near the visible C V line.
	Based on example_lineshape_comparison.py
	"""

	############ DEFINE CASE ############


	shotnumber = str(shotnumber)
	#tracks = ['L1', 'L2','L3'] # tracks of interest
	#tracks = ['L4', 'L5','L6'] # tracks of interest
	tracks = ['T2'] # tracks of interest
	times = [1.5,2.5,3.5,4.5]  # times of interest [s]
	xlim_ROI = [491.0, 498.0] # wavelength of interest [nm]
	#xlim_ROI = [494.0, 495.0] # wavelength of interest [nm]
	#wvl_offset = -17.5  # wavelength offset [A]
	wvl_offset = 0.0  # wavelength offset [A]


	############ GET DATA #############


	data = mdspec_py.data.spectra(shotnumber)
	data.preprocess(wavelength_offset=wvl_offset, background_subtract=True, neutron_filt=True)

	# find array location for each time of interest
	times_loc = []
	for time in times:
		times_loc.append(abs(data.times - time).argmin())

	############ PLOTTING ############

	#First change matplotlib defaults for aesthetics
	colors = np.linspace(50,205,len(times_loc)).astype('i')
	matplotlib.rcParams['font.size']=16
	matplotlib.rcParams['lines.linewidth']=2


	legend_names = []
	for time_index in times_loc:
		legend_names.append('Time = {:.0f} ms'.format(1000*data.times[time_index]))



	plt.figure(facecolor='white', figsize=(18,8))



	for subplot_index,track in zip(range(len(tracks)),tracks):

		plt.subplot(1,len(tracks),subplot_index+1)

		for color_index, time_index in zip(range(len(times_loc)),times_loc):
			plt.plot(data.wavelengths[track]/10,
					 data.intensities[track][time_index,:],
						 color=plt.cm.viridis(colors[color_index]))

		plt.title(data.shotnumber + ': Track ' + track)
		plt.xlim(xlim_ROI)
		plt.ylabel('Intensity [arb]')
		plt.xlabel('Approx wavelength [nm]')
		plt.legend(legend_names, loc='upper right', frameon=False, fontsize=14)
		plt.gca().ticklabel_format(useOffset=False)

	plt.axvline(x=494.388, c='r')
	plt.axvline(x=494.456, c='r')

	wvl_str = 'Wavelength offset = {0} A'.format(wvl_offset)
	plt.text(0.05,0.95,wvl_str,transform=plt.gca().transAxes)

	plt.show(block=False)

def get_integrated_cv(spectra,track):
	"""
	Calculate signal vs time for C V doublet in a given spectra, integrated
	over a range of wavelengths near the lines.
	"""

	wvl_integ=[4944.4,4944.6] # wavelengths to integrate over [A]
	wvl_loc = []
	for wvl in wvl_integ:
		wvl_loc.append(abs(spectra.wavelengths[track] - wvl).argmin())

	wvl_loc_min=np.amin(wvl_loc)
	wvl_loc_max=np.amax(wvl_loc)
	wvl_min=np.amin([spectra.wavelengths[track][wvl_loc_min],spectra.wavelengths[track][wvl_loc_max]])
	wvl_max=np.amax([spectra.wavelengths[track][wvl_loc_min],spectra.wavelengths[track][wvl_loc_max]])

	wvl_range='{0:.2f} to {1:.2f} nm'.format(wvl_min/10,wvl_max/10)
	print 'integrating across wavelength range:', wvl_range

	cv_integ_raw=np.sum(spectra.intensities[track][:,wvl_loc_min:wvl_loc_max],axis=1)

	return cv_integ_raw, wvl_range

def plot_cv_time(track):
	"""
	Plot C V time traces for a set of shots, for the given track.
	"""

	shotlist=[178345,178346,178348,178349,178350]
	#shotlist=[178345,178340]
	wvl_offset=-17.1

	#First change matplotlib defaults for aesthetics
	matplotlib.rcParams['font.size']=16
	matplotlib.rcParams['lines.linewidth']=3
	plt.figure(facecolor='white', figsize=(18,8))

	legend_names = []

	# grab data for each shot, integrate over wavelengths, and add to plot
	for shotnumber in shotlist:
		data = mdspec_py.data.spectra(shotnumber)
		data.preprocess(wavelength_offset=wvl_offset, background_subtract=True, neutron_filt=True)
		cv, wvl_range=get_integrated_cv(data,track)
		plt.plot(data.times,cv)
		legend_names.append(str(shotnumber))



	plt.title('Track ' + track + ', Wavelength Range: '+ wvl_range)
	plt.xlim([0,6])
	plt.grid()
	plt.ylabel('Intensity [arb]')
	plt.xlabel('Time [s]')
	plt.legend(legend_names, loc='upper right', frameon=False, fontsize=14)

	wvl_str = 'Wavelength offset = {0} A'.format(wvl_offset)
	plt.text(0.02,0.95,wvl_str,transform=plt.gca().transAxes)

	plt.show()




# The following functions are based on peak fitting code written in	https://chrisostrouchov.com/post/peak_fit_xrd_python/
def generate_model(spec):
    composite_model = None
    params = None
    x = spec['x']
    y = spec['y']
    x_min = np.min(x)
    x_max = np.max(x)
    x_range = x_max - x_min
    y_max = np.max(y)
    for i, basis_func in enumerate(spec['model']):
        prefix = 'm{0}_'.format(i)
        model = getattr(models, basis_func['type'])(prefix=prefix)
        if basis_func['type'] in ['GaussianModel']: # for now only GaussianModel
            #model.set_param_hint('sigma', min=1e-6, max=x_range) # can't use b/c of bug in current version?
            model.set_param_hint('center', min=x_min, max=x_max)
            model.set_param_hint('amplitude', min=1e-6)

            # default guess is horrible!! do not use guess()
            default_params = {
                prefix+'center': x_min + x_range * random.random(),
                prefix+'height': y_max * random.random(),
                prefix+'sigma': x_range * random.random()
            }
        else:
            raise NotImplemented('model {0} not implemented yet'.format(basis_func["type"]))

        if 'help' in basis_func:  # allow override of settings in parameter
            for param, options in basis_func['help'].items():
                model.set_param_hint(param, **options)
        model.set_param_hint('{0}height'.format(prefix), min=1e-6, max=1.1*y_max, expr='{0}amplitude/{1}sigma/2.50663'.format(prefix,prefix))
        #model_params = model.make_params(**default_params, **basis_func.get('params', {}))
        if 'params' in basis_func:
		    model_params = model.make_params(**basis_func.get('params', {}))
        else:
		    model_params = model.make_params(**default_params)

        if params is None:
            params = model_params
        else:
            params = params + model_params
        if composite_model is None:
            composite_model = model
        else:
            composite_model = composite_model + model
    #params.pretty_print()
    return composite_model, params

def gaussian(x, A, mu, sigma):   # Gaussian
	return A / (sigma * math.sqrt(2 * math.pi)) * np.exp(-(x-mu)**2 / (2*sigma**2))

def test_fit(method='manual'):
	# sandbox function to create toy data, and fit it with various methods

	# method:
	#	definition: plot known gaussians to make sure data generation is working
	#	manual: fit to 2 gaussian model, manually defined and constrained
	#	auto: use generate_model, with constraints defined in spec
	method = method

	# Define toy data: 2 gaussians with noise
	g_0 = [250.0, 4.0, 5.0]
	g_1 = [20.0, -5.0, 1.0]
	n = 150
	x = np.linspace(-10, 10, n)
	y = gaussian(x, *g_0) + gaussian(x, *g_1) + np.random.randn(n)
	maxy=np.max(y)

	if method == 'definition':
		fig, ax = plt.subplots()
		ax.scatter(x, y, s=2)
		ax.plot(x, gaussian(x, *g_0))
		ax.plot(x, gaussian(x, *g_1))
		ax.plot(x, gaussian(x, *g_0)+gaussian(x, *g_1))

	elif method == 'manual':
		model_1 = models.GaussianModel(prefix='m1_')
		model_2 = models.GaussianModel(prefix='m2_')
		model = model_1 + model_2

		# have to add constraints to get decent result
		model_1.set_param_hint('m1_center', min=0,max=10,value=5.0)
		model_1.set_param_hint('m1_sigma', max=20, value=1.0)
		model_1.set_param_hint('m1_amplitude', min=0)
		model_1.set_param_hint('m1_height', expr='m1_amplitude/m1_sigma/2.50663', min=0, max=1.1*maxy)
		model_2.set_param_hint('m2_center', min=-10,max=0, value=-5.0)
		#model_2.set_param_hint('m2_sigma', max=20, value=1.0)
		model_2.set_param_hint('m2_amplitude', min=0)
		model_2.set_param_hint('m2_height', expr='m2_amplitude/m2_sigma/2.50663', min=0, max=1.1*maxy)


		params_1 = model_1.make_params()
		params_2 = model_2.make_params()
		params = params_1 + params_2
		print 'Before peak fitting:'
		params.pretty_print()

		output = model.fit(y, params, x=x)

		print 'After peak fitting:'
		output.params.pretty_print()
		fig = plt.figure(figsize=(8,12))
		ax1 = plt.subplot2grid((3,2),(0,0),colspan=2,rowspan=2)
		ax2 = plt.subplot2grid((3,2),(2,0),colspan=2)
		output.plot_fit(ax=ax1,data_kws={'markersize': 2})
		output.plot_residuals(ax=ax2,data_kws={'markersize': 2})

		components = output.eval_components(x=x)
		ax1.plot(x,components['m1_'])
		ax1.plot(x,components['m2_'])

		print 'Fit report:'
		print output.fit_report()

	elif method == 'auto':
		spec = {
			'x': x,
		    'y': y,
		    'model': [
		        {
					'type': 'GaussianModel',
					'params': {'center': 5.0, 'amplitude': 1, 'fwhm': 1.0, 'sigma':1.0},
					'help': {'center': {'min': 0, 'max': 10}, 'fwhm':{'min':0, 'max':20}, 'sigma':{'max':20}}
				},
		        {
					'type': 'GaussianModel',
					'params': {'center': -5.0, 'amplitude': 1, 'fwhm': 1.0, 'sigma':1.0},
					'help': {'center': {'min': -10, 'max':0}, 'fwhm':{'min':0, 'max':20}, 'sigma':{'max':20}}
				}
		    ]
		}

		model, params = generate_model(spec)
		print 'Before peak fitting:'
		params.pretty_print()

		output = model.fit(spec['y'], params, x=spec['x'])

		print 'After peak fitting:'
		output.params.pretty_print()
		fig = plt.figure(figsize=(8,12))
		ax1 = plt.subplot2grid((3,2),(0,0),colspan=2,rowspan=2)
		ax2 = plt.subplot2grid((3,2),(2,0),colspan=2)
		output.plot_fit(ax=ax1,data_kws={'markersize': 2})
		output.plot_residuals(ax=ax2,data_kws={'markersize': 2})

		components = output.eval_components(x=spec['x'])
		for i, model in enumerate(spec['model']):
			ax1.plot(spec['x'],components['m{0}_'.format(i)])

		print 'Fit report:'
		print output.fit_report()

	else:
		return

	plt.show()





# New, better functions to plot C V spectral data
def plot_data_compare(track='T3', wvl_offset=-17.1, dataset1=[178346,[2.01,2.99]], dataset2=[178346,[4.01,4.99]], xlim_ROI=None):
	# plot two sets of MDS spectral data on the same axes, then plot their difference
	# shots/times/tracks are supplied to function
	# dataset* takes form of [shot,[tmin,tmax]], and the spectrum is averaged over that time interval

	data1 = mdspec_py.data.spectra(dataset1[0])
	data1.preprocess(wavelength_offset=wvl_offset, background_subtract=True, neutron_filt=True)
	time_indices1 = [x for x in range(len(data1.times)) if (data1.times[x]>dataset1[1][0] and data1.times[x]<dataset1[1][1])]
	avg_intensities1 = np.mean(data1.intensities[track][time_indices1],axis=0)

	data2 = mdspec_py.data.spectra(dataset2[0])
	data2.preprocess(wavelength_offset=wvl_offset, background_subtract=True, neutron_filt=True)
	time_indices2 = [x for x in range(len(data2.times)) if (data2.times[x]>dataset2[1][0] and data2.times[x]<dataset2[1][1])]
	avg_intensities2 = np.mean(data2.intensities[track][time_indices2],axis=0)

	if xlim_ROI==None:
		min_wvl=np.min(data1.wavelengths[track])
		max_wvl=np.max(data1.wavelengths[track])
		xlim_ROI=[min_wvl-0.5,max_wvl+0.5]

	fig = plt.figure(facecolor='white', figsize=(18,8))
	ax1 = plt.subplot2grid((4,5),(0,0),colspan=5,rowspan=3)
	ax2 = plt.subplot2grid((4,5),(3,0),colspan=5,sharex=ax1)

	colors = ['r','b','g','c','m','y']

	ax1.plot(data1.wavelengths[track], avg_intensities1, label=dataset1,
							 color=colors[0],marker='o',markersize=4,linewidth=0)
	ax1.plot(data2.wavelengths[track], avg_intensities2, label=dataset2,
							 color=colors[1],marker='o',markersize=4,linewidth=0)
	ax1.set_xlim(xlim_ROI)
	ax1.set_ylabel('Intensity [arb]')
	ax1.set_title('Track ' + track)
	ax1.ticklabel_format(useOffset=False)
	ax1.set_yscale('log')
	ax1.set_ylim([1e1,4e4])
	ax1.legend(loc='upper right', frameon=False, fontsize=14)
	wvl_str = 'Wavelength offset = {0} A'.format(wvl_offset)
	plt.text(0.05,0.95,wvl_str,transform=ax1.transAxes)

	diff = avg_intensities2 - avg_intensities1
	ax2.plot(data1.wavelengths[track], diff, color='k',marker='o',markersize=4,linewidth=0)
	ax2.axhline(y=0, color='k')
	#ax2.set_yscale('symlog', linthreshy=1e2)
	#ax2.set_ylim([-3e3,3e3])
	ax2.set_ylabel('Diff')
	ax2.set_xlabel('Approx wavelength [A]')

	plt.show()

def fit_spectra_single(shot=178346, time=3.5, track='T3', wvl_offset=-17.1, xlim_ROI=None, plot=True, specname='get_spec_4890_4970_limited'):
	# fit peaks to a single MDS spectrum and return fit values
	# optionally plot data and fit peaks
	#

	# get data
	data = mdspec_py.data.spectra(shot)
	data.preprocess(wavelength_offset=wvl_offset, background_subtract=True, neutron_filt=True)

	if xlim_ROI==None:
		min_wvl=np.min(data.wavelengths[track])
		max_wvl=np.max(data.wavelengths[track])
		xlim_ROI=[min_wvl-0.5,max_wvl+0.5]

	# find array location for each time of interest
	time_index = abs(data.times - time).argmin()

	# specify fit parameters
	spec=globals()[specname](data,track,time_index)

	numpeaks = len(spec['model'])
	print 'number of peaks in model: {0}'.format(len(spec['model']))

	# prettify
	colors = np.linspace(0.1,0.9,numpeaks)
	matplotlib.rcParams['font.size']=16
	matplotlib.rcParams['lines.linewidth']=2
	shot_str = str(data.shotnumber)
	time_str = '{:.0f}'.format(data.times[time_index]*1000)

	# do fit
	model, params = generate_model(spec)
	#params.pretty_print()
	output = model.fit(spec['y'], params, x=spec['x'])
	components = output.eval_components(x=spec['x'])
	residuals = spec['y'] - output.eval()
	maxresidual = np.max(np.abs(residuals))
	print 'Fit report:'
	print output.fit_report(show_correl=False)

	# plot
	if plot:
		fig = plt.figure(facecolor='white', figsize=(18,8))
		ax1 = plt.subplot2grid((4,5),(0,0),colspan=5,rowspan=3)
		ax2 = plt.subplot2grid((4,5),(3,0),colspan=5,sharex=ax1)
		ax1.plot(data.wavelengths[track], data.intensities[track][time_index],
							 color='r',marker='o',markersize=4,linewidth=0)

		for i, model in enumerate(spec['model']):
			ax1.plot(spec['x'],components['m{0}_'.format(i)], color=plt.cm.viridis(colors[i]))
		ax1.plot(spec['x'],output.eval(),color='k')

		ax1.set_title(shot_str+'.'+time_str+': Track ' + track)
		ax1.set_xlim(xlim_ROI)
		ax1.set_ylabel('Intensity [arb]')
		ax1.set_xlabel('Approx wavelength [A]')

		ax1.ticklabel_format(useOffset=False)
		ax1.set_yscale('log')
		ax1.set_ylim([1e2,1e5])

		wvl_str = 'Wavelength offset = {0} A'.format(wvl_offset)
		plt.text(0.05,0.95,wvl_str,transform=ax1.transAxes)

		ax2.set_ylabel('Residuals')
		ax2.set_xlabel('Approx wavelength [A]')
		ax2.plot(spec['x'],residuals, color='r',marker='o',markersize=2,linewidth=0)
		ax2.axhline(y=0, color='k')
		ax2.set_ylim([-maxresidual,maxresidual])
		ax2.locator_params(axis='y', tight=True, nbins=4)

		plt.show(block=False)

	# return dictionary of best fit values for parameters
	return output.params.valuesdict()

def fit_spectra_multi(shot=178346, tlim_ROI=[0.0,5.5], track='T3', wvl_offset=-17.1, xlim_ROI=None, plot=True, save=False, savebase='fitparams',specname='get_spec_4890_4970_limited'):
	# fit peaks to a time series of MDS spectra
	# optionally plot peak evolution
	# saves fit values to <savebase>_<shot>_<track>.csv for future plotting if save=True

	# get data
	data = mdspec_py.data.spectra(shot)
	data.preprocess(wavelength_offset=wvl_offset, background_subtract=True, neutron_filt=True)

	# find array location for each time of interest
	time_index_min = abs(data.times - tlim_ROI[0]).argmin()
	time_index_max = abs(data.times - tlim_ROI[1]).argmin()
	time_indices = np.arange(time_index_min, time_index_max, 1)
	print 'number of time points in series: {0}'.format(len(time_indices))

	# specify fit structure from which to build time series
	basespec=globals()[specname](data, track, time_indices[0])

	numpeaks = len(basespec['model'])
	print 'number of peaks in model: {0}'.format(len(basespec['model']))

	# prettify
	colors = np.linspace(0.1,0.9,numpeaks)
	matplotlib.rcParams['font.size']=16
	matplotlib.rcParams['lines.linewidth']=2
	# These are the "Tableau 20" colors as RGB, and scaled to [0,1]
	tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
	for i in range(len(tableau20)):
		r, g, b = tableau20[i]
		tableau20[i] = (r / 255., g / 255., b / 255.)
	linestyles = ['-', '--', '-.', ':']
	shot_str = str(data.shotnumber)
	time_str = '{:.0f}-{:.0f}'.format(data.times[time_indices[0]]*1000,data.times[time_indices[-1]]*1000)
	print shot_str+'.'+time_str

	# initialize dictionary, where the default value is an empty list
	times_fitdict = defaultdict(list)

	# do fit for each time point and append to dictionary lists
	for i, time_index in enumerate(time_indices):
		print 'fitting {}.{:.0f}'.format(shot_str,data.times[time_index]*1000)
		spec = globals()[specname](data, track, time_index)
		model, params = generate_model(spec)
		output = model.fit(spec['y'], params, x=spec['x'])
		fitdict = output.params.valuesdict()
		times_fitdict['time'].append(data.times[time_index])
		times_fitdict['shot'].append(shot)
		times_fitdict['track'].append(track)
		times_fitdict['wvl_offset'].append(wvl_offset)
		times_fitdict['numpeaks'].append(numpeaks)
		for key, value in fitdict.iteritems():
			times_fitdict[key].append(value)

	#print times_fitdict

	if save:
		outfilename = '{0}_{1}_{2}.csv'.format(savebase,shot_str,track)
		keys = sorted(times_fitdict.keys())
		with open(outfilename, 'wb') as outfile:
			writer = csv.writer(outfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
			for key in keys:
				writer.writerow([key]+times_fitdict[key])

	if plot:
		fig = plt.figure(facecolor='white', figsize=(18,8))
		ax1 = plt.subplot(111)
		for i in range(0,len(basespec['model'])):
			ax1.plot(times_fitdict['time'],times_fitdict['m{0}_height'.format(i)], label='m{0}'.format(i),
				marker='o',markersize=4,linewidth=3,linestyle=linestyles[i%4],color=tableau20[i%20])
		ax1.legend(loc='upper right', frameon=False, fontsize=14)
		ax1.set_title(shot_str+': Track ' + track)
		ax1.set_xlim([tlim_ROI[0]-0.1,tlim_ROI[1]+0.5])
		ax1.set_ylabel('Peak Height [arb]')
		ax1.set_xlabel('Time [s]')

		plt.show()

	return times_fitdict

def plot_spectra_multi(shot=178346, track='T3', load=True, loadbase='fitparams', normalize=True, toplot=[6,9,10]):
	# plot time evolution of peak fitting parameters, typically by loading from file named <loadbase>_<shot>_<track>.csv.
	# optionally normalize peak params to the average from a given shot/time with normalize=True.
	# can limit which peaks are plotted with toplot.

	times_fitdict = defaultdict(list)

	fname = '{0}_{1}_{2}.csv'.format(loadbase,shot,track)
	if load and os.path.isfile(fname):
		with open(fname, 'rb') as infile:
			reader = csv.reader(infile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
			for row in reader:
				times_fitdict[row[0]] = row[1:]
	else:
		times_fitdict = fit_spectra_multi(shot, track, plot=False, save=True)

	tmin = np.min(times_fitdict['time'])
	tmax = np.max(times_fitdict['time'])

	numpeaks = int(times_fitdict['numpeaks'][0])


	if normalize:
		# normalize to shot 178345 t=2.5-4.5 s
		norm_name = 'fitparams_178345_{0}.csv'.format(track)
		norm_trange = [2.5,4.5]
		times_fitdict_norm = defaultdict(list)
		with open(norm_name, 'rb') as normfile:
			reader = csv.reader(normfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
			for row in reader:
				times_fitdict_norm[row[0]] = row[1:]
		norm_dict = defaultdict(float)
		# kludgy way to approximate numpy boolean masks for lists
		boolmask = [(t>=norm_trange[0]) and (t<=norm_trange[1]) for t in times_fitdict_norm['time']]
		for key in times_fitdict_norm.keys():
			maskedlist = [times_fitdict_norm[key][i] for i in xrange(len(times_fitdict_norm[key])) if boolmask[i]]
			# average masked values if they are numeric, but take only the first value if not (like track name)
			if type(maskedlist[0])==float or type(maskedlist[0])==int:
				avgvalue = sum(maskedlist)/len(maskedlist)
			else:
				avgvalue = maskedlist[0]
			norm_dict[key] = avgvalue

	# prettify
	colors = np.linspace(0.01,0.99,numpeaks)
	matplotlib.rcParams['font.size']=16
	matplotlib.rcParams['lines.linewidth']=2
	# These are the "Tableau 20" colors as RGB, and scaled to [0,1]
	tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
	for i in range(len(tableau20)):
		r, g, b = tableau20[i]
		tableau20[i] = (r / 255., g / 255., b / 255.)
	linestyles = ['-', '--', ':']
	shot_str = str(shot)
	time_str = '{:.0f}-{:.0f}'.format(tmin*1000,tmax*1000)

	fig1 = plt.figure(facecolor='white', figsize=(18,8))
	ax1 = fig1.add_subplot(111)
	for i in toplot:
		ax1.plot(times_fitdict['time'],times_fitdict['m{0}_height'.format(i)], label='m{0}'.format(i),
			marker='o',markersize=4,linewidth=3,linestyle=linestyles[i%3],color=plt.cm.jet(colors[i]))
	ax1.legend(loc='upper right', frameon=False, fontsize=14)
	ax1.set_title(shot_str+': Track ' + track)
	ax1.set_xlim([tmin-0.1,tmax+0.5])
	ax1.set_ylabel('Peak Height [arb]')
	ax1.set_xlabel('Time [s]')

	fig2 = plt.figure(facecolor='white', figsize=(18,8))
	ax2 = fig2.add_subplot(111)
	for i in toplot:
		ax2.plot(times_fitdict['time'],times_fitdict['m{0}_center'.format(i)], label='m{0}'.format(i),
			marker='o',markersize=4,linewidth=3,linestyle=linestyles[i%3],color=plt.cm.jet(colors[i]))
	ax2.legend(loc='upper right', frameon=False, fontsize=14)
	ax2.set_title(shot_str+': Track ' + track)
	ax2.set_xlim([tmin-0.1,tmax+0.5])
	ax2.set_ylabel('Peak Center [A]')
	ax2.set_xlabel('Time [s]')

	if normalize:
		fig3 = plt.figure(facecolor='white', figsize=(18,8))
		ax3 = fig3.add_subplot(111)
		for i in toplot:
			ax3.plot(times_fitdict['time'],[y/norm_dict['m{0}_height'.format(i)] for y in times_fitdict['m{0}_height'.format(i)]],
				label='m{0}'.format(i), marker='o',markersize=4,linewidth=3,linestyle=linestyles[i%3],color=plt.cm.jet(colors[i]))
		ax3.legend(loc='upper right', frameon=False, fontsize=14)
		ax3.set_title(shot_str+': Track ' + track)
		ax3.set_xlim([tmin-0.1,tmax+0.5])
		ax3.set_ylim([0.0, 2.0])
		ax3.set_ylabel('Peak Height / <Peak Height (178345.2500-4500)>')
		ax3.set_xlabel('Time [s]')

	plt.show()

def fit_spectra_mega():
	# big batch of spectra to fit overnight, save results to csv

	shotlist = [178345, 178346, 178348, 178349, 178350, 178340]
	tracklist = ['T2', 'T3', 'L3']

	for shot in shotlist:
		for track in tracklist:
			fit_spectra_multi(shot=shot, tlim_ROI=[0.0,5.5], track=track, wvl_offset=-17.1, plot=False, save=True, savebase='fitparams',
				specname='get_spec_4890_4970_limited')
	return

# Spectrum fitting specification functions
def get_spec_4910_4980(data, track, time_index):
	# spectrum template for unshifted MDS data 4910-4980 A
	# peak locations eyeballed from spectrum; Tz=100 eV
	# data is a mdspec_py.data.spectra object
	spec = {
			'x': data.wavelengths[track],
		    'y': data.intensities[track][time_index],
		    'model': [
		        {
					'type': 'GaussianModel',
					'params': {'center': 4912.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4911.0, 'max': 4913.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4914.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4913.0, 'max': 4915.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4916.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4915.0, 'max': 4917.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4918.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4917.0, 'max': 4919.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4922.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4921.0, 'max': 4923.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
		        {
					'type': 'GaussianModel',
					'params': {'center': 4924.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4923.0, 'max': 4925.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4928.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4927.0, 'max': 4929.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4934.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4933.0, 'max': 4935.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
		        {
					'type': 'GaussianModel',
					'params': {'center': 4935.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4934.0, 'max': 4936.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4939.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4938.0, 'max': 4940.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4941.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4940.0, 'max': 4942.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4942.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4941.0, 'max': 4943.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4947.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4946.0, 'max': 4948.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4950.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4949.0, 'max': 4951.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4954.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4953.0, 'max': 4955.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4957.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4956.0, 'max': 4958.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4960.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4959.0, 'max': 4961.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4962.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4961.0, 'max': 4963.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4968.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4967.0, 'max': 4969.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4973.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4972.0, 'max': 4974.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				},
				{
					'type': 'GaussianModel',
					'params': {'center': 4978.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'min': 4977.0, 'max': 4979.0}, 'fwhm':{'min':0.0}, 'sigma':{'vary':False}}
				}
		    ]
		}
	return spec

def get_spec_4890_4970_nist(data, track, time_index):
	# spectrum template for -17.1 A shifted MDS data 4890-4970 A
	# peak locations from NIST database
	# data is a mdspec_py.data.spectra object
	spec = {
			'x': data.wavelengths[track],
		    'y': data.intensities[track][time_index],
		    'model': [
		        { # O II
					'type': 'GaussianModel',
					'params': {'center': 4890.9, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O IV
					'type': 'GaussianModel',
					'params': {'center': 4902.3, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # Ar II
					'type': 'GaussianModel',
					'params': {'center': 4904.7, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4906.8, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # He I
					'type': 'GaussianModel',
					'params': {'center': 4910.7, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # B III
					'type': 'GaussianModel',
					'params': {'center': 4917.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # B III
					'type': 'GaussianModel',
					'params': {'center': 4918.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
		        { # He I
					'type': 'GaussianModel',
					'params': {'center': 4920.6, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # He I
					'type': 'GaussianModel',
					'params': {'center': 4921.9, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4924.5, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # Ar II
					'type': 'GaussianModel',
					'params': {'center': 4933.2, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # B II
					'type': 'GaussianModel',
					'params': {'center': 4940.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4941.1, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4943.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4955.7, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # Ar II
					'type': 'GaussianModel',
					'params': {'center': 4965.1, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				}
		    ]
		}
	return spec

def get_spec_4890_4970_limited(data, track, time_index):
	# spectrum template for -17.1 A shifted MDS data 4890-4970 A
	# selected reference peaks plus mystery peaks of interest
	# data is a mdspec_py.data.spectra object
	spec = {
			'x': data.wavelengths[track],
		    'y': data.intensities[track][time_index],
		    'model': [
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4906.8, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # B III
					'type': 'GaussianModel',
					'params': {'center': 4917.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # B III
					'type': 'GaussianModel',
					'params': {'center': 4918.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
		        { # He I
					'type': 'GaussianModel',
					'params': {'center': 4920.6, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # He I
					'type': 'GaussianModel',
					'params': {'center': 4921.9, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4924.5, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # mystery, C I?
					'type': 'GaussianModel',
					'params': {'center': 4930.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # B II
					'type': 'GaussianModel',
					'params': {'center': 4940.4, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4943.0, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # mystery, C V or B V?
					'type': 'GaussianModel',
					'params': {'center': 4944.5, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # mystery, C V?
					'type': 'GaussianModel',
					'params': {'center': 4950.9, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				},
				{ # O II
					'type': 'GaussianModel',
					'params': {'center': 4955.7, 'amplitude': 1, 'fwhm': 0.5, 'sigma':0.47},
					'help': {'center': {'vary':False}, 'fwhm':{'min':0.0}, 'sigma':{'min':0.2, 'max':0.7}}
				}
		    ]
		}
	return spec
