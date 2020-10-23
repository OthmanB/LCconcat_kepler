'''
---------------
This file regroups all the basic necessary function to compute a clean lightcurve and power spectrum
from the individual quarters of Kepler observations. 
Use do_LC() and do_TF() in order to do so. Other functions are just dependencies of those.
---------------
'''

from astropy.io import fits
from astropy.timeseries import LombScargle
import astropy.units.si as usi
import astropy.units.cds as ucds
import os, fnmatch
import numpy as np
import matplotlib.pyplot as plt
#from scipy.signal import savgol_filter # for smoothing the spectrum

from scipy.ndimage import gaussian_filter

def version():
	return 'v0.25'

def format_ID(ID, Ndigits=9):
	'''
		Small function that ensure that the ID number is on Ndigits digits
		If this is not the case, we fill with 0 before the existing ID
		With Kepler data, Ndigits=9 is fine. For TESS, Ndigits=10 shoud be OK (need check)
	'''
	NID=len(ID)
	delta=8-NID+1 # We want to format using 8 digits
	New_ID=ID
	for d in range(delta):
		New_ID='0' + New_ID
	return New_ID

def readfits_MAST_kepler(file_in, getraw=False, substract_t0=True):
	'''
	function that extract the time and flux from a fits file downloaded in the MAST database
	'''
	d=fits.open(file_in)
	# Note: Use d.info() and d[1].header to understand the data structure
	# d[0] Contains some stellar information such as the Teff etc...
	# d[1] is the lightcurve with a lot of information that can be revealer looking at the header
	# d[2] Contains some position information, such as the RA and DEC
	ID=d[0].header['OBJECT']
	ID=ID.split()[1] # Remove the 'KIC' and keep only the Number
	ID=format_ID(ID)
	t=d[1].data['TIME']   # This is nearly identical to the KASOC TIME. But there is a slight difference of the order of the time error (don't understand why)
	#tcorr=d[1].data['TIMECORR']  # TIME CORRECTION
	lc_sap=d[1].data['SAP_FLUX'] # This is the RAW FLUX IN KASOC LIGHTCURVES
	lc_pdcsap=d[1].data['PDCSAP_FLUX'] # This is the CORRECTED FLUX IN KASOC LIGHTCURVES
	if getraw == False:
		lc=np.asarray(lc_pdcsap, dtype=np.float)
	else:
		lc=np.asarray(lc_sap, dtype=np.float)
	if substract_t0 == False:
		time=np.asarray(t, dtype=np.float)
	else:
		time=np.asarray(t-min(t), dtype=np.float)
	return ID, time, lc
	
def rem_trend(time, flux, pol_order, ppm=True, debug=False):
	'''
		Perform a polynomial fit of the  flux over time
		Option allows you to decide whether you convert the result in ppm or not
	'''
	#t1,a1=keep_finite(time, flux) # remove bad points (nan and inf)
	fit_A=np.polyfit(time-np.min(time), flux, pol_order)
	fit_A_model=np.polyval(fit_A, time-np.min(time))	
	if ppm == True:
		flux_out=1e6*((flux/fit_A_model) -1.) # normalisation in 
	else:
		flux_out=flux-fit_A_model + np.mean(flux)
 	#r3=reform(result[2,where(finite(r2_0) ne 0)]) # FROM IDL CODE
 	#r3=1d6*(r3/fit_A_model) ; propagate the error on the flux (FROM IDL CODE)
	if debug == True:
		print('---')
		print('min(flux)=', np.min(flux))
		print('max(flux)=', np.max(flux))
		print('mean(flux)=', np.mean(flux))
		print('var(flux)=', np.var(flux))
		#
		print('min(flux_out)=', np.min(flux_out))
		print('max(flux_out)=', np.max(flux_out))
		print('mean(flux_out)=', np.mean(flux_out))
		print('var(flux_out)=', np.var(flux_out))
		print('fit_A= ', fit_A)
		print('---')
	return time, flux_out, fit_A_model

def eval_con_quarters(timeA, fluxA, timeB, fluxB, Dt_connect=5):
	pol_order=1
	posA=np.where(timeA >= (max(timeA)-Dt_connect)) 
	fit_A=np.polyfit(timeA[posA], fluxA[posA], pol_order) # we take only the last days
	fit_A_model=np.polyval(fit_A, timeA[posA]) #Just for compatibility with graphs
	#fit_B=np.polyfit(timeB[posA]-np.min(timeB), fluxB[posA], pol_order) # we take only the last days
	#fit_B_model=np.polyval(fit_B, timeB[posA]-np.min(timeB[posA])) #Just for compatibility with graphs
	delta=fit_A[1] + fit_A[0]*timeB[0] #  new initial position of the B curve
	#print(timeB[0] - max(timeA))
	#print(fit_A)
	new_fluxB=fluxB- (fluxB[0] - delta)
	#print(delta)
	#print(new_fluxB[0], fluxB[0])
	#plt.plot(timeA, fluxA, color='black')
	#plt.plot(timeA[posA], fit_A_model, color='blue')
	#plt.plot(timeB, new_fluxB, color='red')
	#plt.show()
	#exit()
	return new_fluxB, timeA[posA], fit_A_model#, timeB[posA], fit_B_model

def keep_finite(time, flux):
	'''
		This remove any bad point such as NaN and Inf
	'''
	a1=flux[np.isfinite(flux)] # remove bad points (nan and inf)
	t1=time[np.isfinite(flux)]
	return t1,a1

def get_files_list(dir_in, extension='.fits', prefix='kplr', ID='*'):
	'''
	This routine returns all files that have some extension and prefix and optionally a specific kic encoded in a 9-numerical byte format
	dir_in: The directory to look in for files. This must be a full path
	extension: The extension to look for. This is designed to provide files for prepare_lc_kepler::concatenate_kepler() so that it is fits by default
	prefix: the prefix that defines the names of the files. Kepler data from MAST and KASOC are preceded of 'kplr', which is the default prefix
	ID: The ID number of a specific star. Should not be set if the data structure is 1-star, 1 directory. Instead use the default
	'''
	listOfFiles = os.listdir(dir_in)
	if ID != '*': 
		if extension == '.fits' and prefix == 'kplr': # Ensure that the search incorporate the date sufix that is always in the quarter filenames 
			ID=ID+'-*'
		else:
			ID=ID +'*' # Case where the function is used in another context than Kepler quarters
	pattern = prefix + ID + extension
	matching=[]
	for entry in listOfFiles:
		if fnmatch.fnmatch(entry, pattern):
			matching.append(entry)
			#print(entry)
	return matching

def concatenate_kepler(dir_in, files_in_list, dir_out, remove_trend=True, pol_order=6, var_calibration=False, ignore_bigGaps=False, useraw=False, doplots=False, debug=False):
	'''
		Main program handling the concatenation of Kepler quarters. It is a direct adaptation of the IDL program (find_target_KEPLER_v2015.pro)
		One difference however: The IDL program had the ability to detect some large gaps and then create separate lightcurves if the gaps was bigger than a certain
		threshold. The default was an excessive threshold so that this was never happenning. Hence, the removal. 
		dir_in: Directory that contains the searched files
		files_in_list: list of fits files that have to be scanned and concatenated. Obviously, this list should be for a single KIC star and it is
					   up to the user to be careful to ensure that (no failsafe)
		dir_out: Directory in which the concatenated files (.npz and jpg if doplots = True) will be written
		remove_trends: If True, remove long trends by performing a polynomial fit of order defined by pol_order (default = 2)
		pol_order: The polynomial order used to perform the trend removal
		var_calibration: If True, the variance/standard deviation of each quarter is modified to match the lower one. This was introduced in an attempt to deal
						 with quarters that have very large variance as it happens sometimes. *** BE CAREFULL WITH THIS AS IT IS AN ALTERATION OF THE DATA ***
		ignore_bigGaps: If True, ignore gaps which are very large (e.g. if there is missing quarters) by sticking two segments that are far appart. This allows to mitigate
					    aliasing effects, but is based on the asumption that the phase information is not critical (e.g. for solar-like)
					    WARNING: ignore_biggaps is applied only if remove_trends = True
		useraw: If False, it will concatenate the corrected lightcurve (recommended). Otherwise, it is using the uncorrected lightcurve
		doplots: To get plots into files
	'''
	col=['plum', 'darkgray', 'gray', 'brown', 'cyan', 'purple', 'pink','red', 'orange', 'yellow', 'gold', 'magenta', 
	'green', 'darkgreen', 'deeppink' , 'aquamarine', 'olive', 'springgreen','greenyellow', 'skyblue', 'darkgoldenrod']
	for i in range(10):
		col=col+col
	Nfiles=len(files_in_list)
	# case of multiples file to concatenate. In that case, dir is included in the all_filenames table
	t0=np.zeros(Nfiles, dtype=np.float)
	tmax=np.zeros(Nfiles, dtype=np.float)
	amp_min=np.zeros(Nfiles, dtype=np.float)
	amp_max=np.zeros(Nfiles, dtype=np.float)
	timegroup_min=np.zeros(1, dtype=np.float)
	timegroup_max=np.zeros(1, dtype=np.float)
	print(' ------- Configuration ---------')
	print(' Data Directory       : ', dir_in)
	print(' Output Directory     : ', dir_out)
	print(' Removing trends      : ', remove_trend)
	if remove_trend == True:
		print('          Order of the polynom : ', pol_order)
	print(' Variance calibration                 : ', var_calibration)
	print(' Ignoring lare gaps (>1 day)          : ', ignore_bigGaps)
	if useraw == False:
		print(' Using corrected/uncorrected data     :  Corrected Data')
	else:
		print(' Using corrected/uncorrected data : Uncorrected Data')	
	print(' Doing plots                          : ', doplots)
	print(' Number of files to process           : ', Nfiles)
	print(' -------------------------------')
	#
	infos='File prepared using prepare_lc_kepler.py ' + version() + ' function concatenate_kepler(), see https://github.com/OthmanB/'
	config={'pol_order':pol_order, 'data_dir':dir_in, 'remove_trend':remove_trend, 'var_calibration':var_calibration, 'ignore_biggaps':ignore_bigGaps, 'useraw':useraw}
	short_config=str(int(useraw == True)) + str(int(remove_trend == True)) + str(int(var_calibration == True)) + str(int(ignore_bigGaps == True))
	#	
	print(' [1] Scanning files and setting up')
	plt.xlabel('time (days)', fontsize=10)
	plt.ylabel('Intensity', fontsize=10)
	for i in range(Nfiles):
		kic_number, t, lc=readfits_MAST_kepler(dir_in + files_in_list[i], getraw=useraw, substract_t0=False)
		t,lc=keep_finite(t,lc)
		t0[i]=np.min(t) # identifying the ordering of the files by extracting start / end time
		tmax[i]=np.max(t) # identifying the ordering of the files by extracting start / end time
		amp_min[i]=np.nanmin(lc) # identifying the minimum amplitude... usefull for plots only
		amp_max[i]=np.nanmax(lc) # identifying the maximum amplitude... usefull for plots only
	file_core=kic_number + '_' + str(int(useraw == True)) + str(int(remove_trend == True)) + str(int(var_calibration == True)) + str(int(ignore_bigGaps == True))
	if Nfiles == 1:
		ts=t # may not be required
		lcs=lc # may not be required
		r1=ts
		r2=lc
		if doplots == True:
			plt.plot(r1-r1[0], r2)
		print(' [2] Saving result')
		lightcurve=np.zeros((2,len(r1)), dtype=np.float)
		#file_out=dir_out + kic_number + '.npz'
		file_out=file_core +  '_LC.npz'
		np.savez_compressed(file_out, lightcurve=lightcurve, id_number=kic_number, config=config, short_config=short_config, infos=infos)
	else:
		files=files_in_list
		t=t0
		t0_new=t0
		tmax_new=tmax
		tmax2=tmax
		# Sorting files by ascending date, using the iterative removal method 
		sorted_files=[]
		i=0
		while len(files) !=1:
			posmin=np.where(t == min(t))[0][0]
			sorted_files.append(files[posmin])
			t0_new[i]=t[posmin]
			tmax_new[i]=tmax2[posmin]
			posnotmin=np.where(t != min(t))
			t2=t[posnotmin]
			tmax2=tmax[posnotmin]
			files=[files[i] for i in posnotmin[0]] 
			t=t2
			i=i+1
		sorted_files.append(files[0]) # We should not forget the last file...
		t0=t0_new
		tmax=tmax_new
		files=[]
		for f in sorted_files:
			files.append(dir_in + f)
		files_in_list=files
		# ---
		#--- Remark: These lines are a remainder of the part that was handling division into several LC if gaps were large in the IDL code ---
		# They are kept for ensuring compatibility (simplification could be done, but would require a full check of the code... again)
		files_packet=files_in_list		
		filegroup=files_in_list
		timegroup_min=np.asarray(t0, dtype=np.float)
		timegroup_max=np.asarray(tmax, dtype=np.float)
		# -------------------
		Nfilegroup=len(filegroup)
		delta=np.zeros(Nfilegroup)
		print(' [2] Reading elements and concatenating...')
		if Nfilegroup <= 0:
			print('ERROR: Nfilegroup is less than 1. Debug required')
			exit()
		if Nfilegroup == 1:
			kic, r1, r2=readfits_MAST_kepler(filegroup, getraw=useraw, substract_t0=False)
			r1,r2=keep_finite(r1,r2)
			if remove_trend == True:
				r1,r2, fit_A_model=rem_trend(r1, r2, pol_order, ppm=True, debug=False)
		if Nfilegroup >1:
			for i in range(Nfilegroup-1):
				print('------------------------')
				print('      Element' , i+1)
				if i == 0:
					kic, r1, r2=readfits_MAST_kepler(filegroup[i], getraw=useraw, substract_t0=False)
					r1,r2=keep_finite(r1,r2)
					kic, r1B, r2B=readfits_MAST_kepler(filegroup[i+1], getraw=useraw, substract_t0=False)
					r1B, r2B=keep_finite(r1B,r2B)
					t10=r1	# For plot purpose only
					a10=r2  # For plot purpose only
					if remove_trend == True:
						r1,r2,fit_A_model=rem_trend(r1, r2, pol_order, ppm=True, debug=False)
						if var_calibration == 1:
							var_tab=np.zeros(1) + np.nanstd(r2)
						if ignore_bigGaps == True:
							timegroup_min=np.zeros(1) + np.min(r1)
							timegroup_max=np.zeros(1) + np.max(r1)
					rfinal1=r1
					rfinal2=r2
				else:
					r1=r1B
					r2=r2B
					kic, r1B, r2B=readfits_MAST_kepler(filegroup[i+1], getraw=useraw, substract_t0=False)
					r1B,r2B=keep_finite(r1B,r2B)
					t10=r1	# For plot purpose only
					a10=r2  # For plot purpose only
				t1,a1=keep_finite(r1,r2)
				t2,a2=keep_finite(r1B,r2B)
				if remove_trend == False:
					Dt_connect=10 # Number of days taken to evaluate how to connect successive segments
					r2B, t1_model, fit_A_model=eval_con_quarters(t1, a1, t2, a2, Dt_connect=Dt_connect)
					rfinal1=np.concatenate((rfinal1, r1B))
					rfinal2=np.concatenate((rfinal2, r2B))
					t10=t1 # We save this only for graphical purpose
					a10=a1 # We save this only for graphical purpose
				else:
					resB_old=r2B
					r1B,r2B, fit_B_model=rem_trend(t2, a2, pol_order, ppm=True, debug=False) # COULD BE SOME ERROR HERE (r1B, r2B =...)
					#t2_model=t2
					if ignore_bigGaps == True: # In that case we stick lightcurves when we have gaps bigger than threshold_t
						Threshold_t=1 # in days
						if np.max(r1B) - np.max(r1) >= Threshold_t:
							r1BminOLD=np.min(r1B)
							r1BmaxOLD=np.max(r1B)
							r1B=r1B - min(r1B) + max(r1) + (r1[len(r1)-1] - r1[len(r1)-2])
							timegroup_min=np.concatenate( (timegroup_min, np.zeros(1) + min(r1B) )) # update all timetables
							timegroup_max=np.concatenate( (timegroup_max, np.zeros(1) + max(r1B) )) # update all timetables	
					rfinal1=np.concatenate((rfinal1, r1B))
					rfinal2=np.concatenate((rfinal2, r2B))	
					if var_calibration == 1:
						var_tab=np.concatenate( (var_tab, np.zeros(1) + np.nanstd(r2B)) )
				if doplots == True:
					if i == 0:
						plt.ylim=[min(amp_min), max(amp_max)]
						plt.xlim=[min(t0) , max(tmax)]
						plt.plot(t10, a10, color='black')
						if remove_trend == True:
							plt.plot(t1, fit_A_model, color='blue')  # an oplot
					plt.plot(t2, a2, color=col[i])
					if remove_trend == True:
						plt.plot(t2, fit_B_model, color='blue')
				if debug == True:
					# --- Debug ---
					print('min/max(r1)=', np.min(r1), np.max(r1))
					print('min/max(r1B)=', np.min(r1B), np.max(r1B))
					print('min/max(rfinal1)=', np.min(rfinal1), np.max(rfinal1))
					print('mean(r2)=', np.mean(r2),   'var(r2)=', np.var(r2))
					print('mean(r2B)=', np.mean(r2B),   'var(r2B)=', np.var(r2B))
					print('mean(rfinal2) =', np.mean(rfinal2),   'var(rfinal2)=', np.var(rfinal2))
					print('min(r2)=', np.min(r2),   'max(r2)=', np.max(r2))
					print('min(r2B)=', np.min(r2B),   'max(r2B)=', np.max(r2B))
					print('min(rfinal2) =', np.min(rfinal2),   'max(rfinal2)=', np.max(rfinal2))
			if var_calibration == 1:
				new_var=np.zeros( len(var_tab) ) 
				var_ref=min(var_tab) # the reference variance is the minimum variance
				for rescale in range(len(var_tab)):
					pos0=np.where(np.bitwise_and(r1 >= timegroup_min[rescale], r1 <= timegroup_max[rescale]))
					r2[pos0] = r2[pos0]*var_ref/var_tab[rescale]
					new_var[rescale]=np.nanstd(r2[pos0])	
		print('      Element' , Nfilegroup)
		print('------------------------')
		if doplots == True:
			plt.savefig(dir_out + file_core +  '_LC_fit.png', dpi=300) # Save the previous figure
			plt.xlabel('time (days)', fontsize=10)
			plt.ylabel('Intensity', fontsize=10)
			plt.figure()
			plt.plot(rfinal1-rfinal1[0], rfinal2, linewidth=0.5)
			plt.savefig(dir_out + file_core +  '_LC.png', dpi=300)
			
		print(' [3] Saving results...')
		lightcurve=np.zeros((2,len(rfinal1)), dtype=np.float)
		lightcurve[0,:]=rfinal1
		lightcurve[1,:]=rfinal2
		file_out=dir_out + file_core +  '_LC.npz'
		np.savez_compressed(file_out, lightcurve=lightcurve, id_number=kic_number, config=config, short_config=short_config, infos=infos)
		plt.close() # Closing the plot to avoid excessive memory usage
		print(' Root names of created lightcurve files: ', file_core, '*.*')
	return lightcurve

def compute_ls(time, flux, Nyquist_type='mean'):
	'''
		Core function performing the Lomb-Scargle
		Using time and flux, it compute the power over the proper frequency range
		and at the proper resolution (no oversampling). Then it normalizes it using the Shanon theorem
		so that units are in ppm^2/microHz. 
	'''
	passed=False
	dt=time[1]-time[0]
	if Nyquist_type == 'best':
		for i in range(1, len(time)):
			if dt > time[i]-time[i-1]: # sampling period
				dt=time[i] - time[i-1] # We try to find the highest resolution yet maximising independence of data
		passed=True

	if Nyquist_type == 'linfit':
		fit=np.polyfit(np.linspace(0, len(time)-1, len(time)), time-np.min(time), 1)  # We compute the mean time spacing through linear fiting. Useful for serie with a lot of gaps
		dt=fit[0]
		passed=True

	if (passed == False) or (Nyquist_type == 'mean'): # Mean is the default even if the wrong value of Nyquist_type was given
		dt=(max(time) - min(time))/len(time) # The Nyquist is evaluated using the mean delta time. This is the most common standard (from Press & Rybicki)
		passed=True
	#
	Nyquist=1e6 / dt / 86400. /2 # Nyquist.
	resol=1e6/(max(time) - min(time))/86400.
	#
	freq=np.linspace(0, Nyquist, np.ceil(Nyquist/resol))
	freq[0]=resol/100. # avoid 0 calculation at 0 by approximating it with resol/100.
	power = LombScargle(time*86400.*usi.second, flux).power(freq*1e-6*ucds.Hz) #.autopower(nyquist_factor=1, method='cython')
	Nflux=len(flux)
	#
	tot_MS=np.sum( (flux - np.mean(flux))**2) / np.float(Nflux) # Total power in the timeserie
	tot_lomb=np.sum(power[0:Nflux -1]) # Total power in the lomb-scargle
	#
	bw=freq[2]-freq[1]
	power=power*(tot_MS/tot_lomb)/bw # Normalised spectrum in ppm^2/microHz

	return freq, power,resol

def simple_filter(time, flux, sigma=4):
	''' 
		Filters the lightcurve by removing signals that are above some sigma threshold
	'''
	mean=np.mean(flux)
	stddev=sigma*np.std(flux)
	pos_ok=np.where(np.abs(flux) <= mean + sigma*stddev)
	print('--- Using simple filter ---- ')
	print(' min/max flux:', np.min(flux), ' / ', np.max(flux))
	print(' mean/var flux:', np.mean(flux), ' / ', np.std(flux))
	print(' Total number of data points: ', len(flux))
	print(' Number of rejected points  : ', len(flux) - len(pos_ok[0]))
	print(' Fraction of rejected points: ', (len(flux) - len(pos_ok[0]))/len(flux) * 100), ' %'
	print(' ---------------------------')
	t=time[pos_ok]
	f=flux[pos_ok]
	return t,f

def do_tf_ls(dir_file, filename, doplots=True, planets=False):
	''' 
		Performs the LombScargle periodogram and normalise it in ppm^2/microHz
		Note that lightcurves are assumed to be:
			- saved into an npz format
			- being in a 2D array
			- with the first index having the time in days
			- and the second index having the flux in ppm
	'''
	d=np.load(dir_file+filename)
	id_number=str(d['id_number'])
	numeral_config=str(d['short_config']) # used to properly define the output filename

	time=d['lightcurve'][0,:]
	flux=d['lightcurve'][1,:]
	
	Nyquist_type='mean'
	# ---- Filetering and computing the LS ----
	if planets == True:
		print('The proper function to filter planetary / Binary signal to reveal pulsations')
		print('Is not written yet. Please use planets == False for the moment')
	else:
		time, flux=simple_filter(time, flux, sigma=3.5)
		freq, power,resol=compute_ls(time, flux, Nyquist_type=Nyquist_type) 

	file_out=dir_file + id_number + '_' + numeral_config + '_TF'
	
	# ---- Making proper plots ----
	sfactor1_mu=1.
	sfactor2_mu=0.1
	sfactor1=np.int(sfactor1_mu/resol) # Smooth over 1 microHz
	sfactor2=np.int(sfactor2_mu/resol) # Smooth over 10 microHz
	y_smooth1=gaussian_filter(power, sfactor1, mode='mirror')
	y_smooth2=gaussian_filter(power, sfactor2, mode='mirror')
	
	plt.figure()
	plt.xlabel('Frequency ($\mu$' + 'Hz)', fontsize=10)
	plt.ylabel('Power ($ppm^2/ \mu$'+'Hz)', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	plt.plot(freq, power, color='black', label='spectrum')
	plt.plot(freq, y_smooth1, color='blue', label='sfactor1='+ str(sfactor1_mu) + '($\mu$Hz)')
	plt.plot(freq, y_smooth2, color='red', label='sfactor2='+ str(sfactor2_mu) + '($\mu$Hz)')
	plt.legend(loc='lower left')
	plt.savefig(file_out + '.png', dpi=300) # Save the previous figure
	
	# ---- Saving into npz format data and metadata ----
	infos='File prepared using prepare_lc_kepler.py ' + version() + '  function do_tf_ls(), see https://github.com/OthmanB/'
	config={'Nyquist_type':Nyquist_type, 'dir_file':dir_file, 'filename':filename}
	np.savez_compressed(file_out +'.npz', freq=freq, power=power, id_number=id_number, config=config, infos=infos)
	print(' Root names of created spectrum files: ', dir_file + id_number, '*.*')

	return freq, power

def do_LC(dir_in, dir_out, ID='*'):
	'''
		Main function to concatenate quarters together and compute the full lightcurve
		The function looks into dir_in for files matching the expected syntax, including the ID number (if provided)
		It then compute the lightcurve with and without having removed a trend and based on corrected data. Removing the trend is 
		inapropriate if one wants to extract the periods at low frequencies which justifies the creation of two lightcurves
		Output files (plots and npz files) are stored in the dir_out directory.
		WARNING on ID: If ID='*' make sure that the dir_in does contain only quarters for a single star. Otherwise, it will end up mixing stars!
	'''
	files_in_list=get_files_list(dir_in, extension='.fits', prefix='kplr', ID=ID)
	lc_0100=concatenate_kepler(dir_in, files_in_list, dir_out, remove_trend=True, pol_order=6, var_calibration=False, ignore_bigGaps=False, useraw=False, doplots=True)
	lc_0000=concatenate_kepler(dir_in, files_in_list, dir_out, remove_trend=False, pol_order=2, var_calibration=False, ignore_bigGaps=False, useraw=False, doplots=True)
	print('    Two lightcurve created with short configuration 0000 and 0100')
	print(' Note on filenames: The files syntax follows a strict format defined as:')
	print('           [ID Number]_[useraw][remove_trend][var_calibration][ignore_bigGaps]_*')
	return lc_0100, lc_0000

def do_TF(dir_lc, ID='*'):
	'''
		The Main function that generate a Power Spectrum using Lomb-Scargle
		dir_lc: directory in which the lightcurve is provided following the
				standard that I use (see do_LC()). If a file is not found,
			    an error will be returned. If multiple files are found, one tf
			    is made for each of them.
		ID: Optional filter using an Identifier in case different stars are in the same directory
	'''
	files_in_list=get_files_list(dir_lc, extension='_LC.npz', prefix='', ID=ID)
	print(files_in_list)
	for file in files_in_list:
		do_tf_ls(dir_lc, file, doplots=True, planets=False)


def do_test_LC():
	current_dir=os.getcwd()
	dir_in=current_dir + '/test/'
	dir_out=dir_in + 'out/'
	do_LC(dir_in, dir_out)

def do_test_tf_ls():
	current_dir=os.getcwd()
	dir_in=current_dir + '/test/out/'
	files_in_list=get_files_list(dir_in, extension='_LC.npz', prefix='', ID='*')
	do_tf_ls(dir_file, filename, doplots=True, planets=False)

	print('dir_in =', dir_in)
	print('files_in_list: ')
	print(files_in_list)
	do_tf_ls(dir_in, files_in_list[0])

def do_test_TF():
	current_dir=os.getcwd()
	dir_lc=current_dir + '/test/out/'
	do_TF(dir_lc, ID='*')

print('prepare_lc_kepler.py, version ', version(), ' loaded')

# ----- Testing programs ----
#do_test_LC() # testing the generation of a lightcurve with a provided example
#do_test_tf_ls()
#do_test_TF()
