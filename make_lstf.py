'''
---------------
Here it is high level functions that are based on do_LC() and do_TF() that provide power spectra.
Typically the functions bellow search for the location of some kic stars
extract their lightcurve and or spectrum and organize them in structured manner (or not)
---------------
'''
from scan_mastdir import scandir_mast_kepler_bylist, scandir_mast_kepler_byfile
from prepare_lc_kepler import do_LC, do_TF
from termcolor import colored

def version():
	return 'v0.1'

def make_LC(dir_base, cadence, dir_base_out, kic_list='', kic_file='', kic_file_col=0):#, force_recompute=False):
	''' 
		Using a provided list of kic number (kic_list), the code (1) search for the location 
		of the quarters within dir_base and then (2) concatenate the quarters to make the 
		final lightcurves as specified by do_LC() in prepare_lc_kepler.py
		dir_base: root directory with the data downloaded from MAST
		cadence: Cadence of the data. only 'lc' or 'sc' is allowed.
		dir_base_out: root directory where the output lighturves will be put
		kic_list: List of kic number for which the processing is requested. Cannot be used jointly with kic_file
		kic_file: List of kic number specified into a text file. Cannot be used jointly with kic_list
		kic_file_col: Column to be considered to extract the kic list when reading the kic_file text file
		##force_recompute: If True, recompute the lightcurves even if existing concatenated lighturves are found
		##				 If False, will skip the stars for which npz files associated to concatenated lightcurves are found (NOT IMPLEMENTED)
	'''

	print('[1] Locating the Kepler quarters...')
	if kic_list != '' and kic_file == '':
		paths, found_kic, missing_kic=scandir_mast_kepler_bylist(dir_base, kic_list, cadence=cadence)
	if kic_list == '' and kic_file != '':
		paths, found_kic, missing_kic=scandir_mast_kepler_byfile(dir_base, kic_file, kics_col=kic_file_col, cadence=cadence)
	if kic_list != '' and kic_file != '':
		print('Error: Both kic_list and kic_file are set. This is forbidden. Please use one or the other.')
		exit()
	if kic_list == '' and kic_file == '':
		print('Error: Please provide one of the following argument:')
		print('      [1]  kic_list : A list containing kic numbers')
		print('      [2]  kic_file : A text file containing kic numbers. Consider also the validity')
		print('           of kic_file_cols as this specifies the column that has to be read for the kic')
		exit()
	#print('[2] Verifying that LC does not already exist for the specified KIC stars...')
	print('[3] Computing the Lightcurves')
	l=len(paths)
	for i in range(l):
		print(' [' + str(i+1) + '/' + str(l) + '] Processing KIC ', found_kic[i], '...') 
		try:
			do_LC(paths[i] +'/', dir_base_out)
		except Exception as e:
			print(colored('Error: Could not process KIC ', found_kic[i]),"yellow")
			#exit()
	return found_kic

def make_ALL(dir_base_out, dir_base='../lightcurves/', cadence='sc', kic_list='', kic_file='', kic_file_col=0):
	''' 
		Using a provided list of kic number (kic_list), the code (1) search for the location 
		of the quarters within dir_base and then (2) concatenate the quarters to make the 
		final lightcurves as specified by do_LC() in prepare_lc_kepler.py and (3) finally
		compute the power spectrum using lombscargle
		dir_base: root directory with the data downloaded from MAST
		cadence: Cadence of the data. only 'lc' or 'sc' is allowed.
		dir_base_out: root directory where the output lighturves will be put
		kic_list: List of kic number for which the processing is requested. Cannot be used jointly with kic_file
		kic_file: List of kic number specified into a text file. Cannot be used jointly with kic_list
		kic_file_col: Column to be considered to extract the kic list when reading the kic_file text file
	'''
	valid_kic=make_LC(dir_base, cadence, dir_base_out + '/' + cadence + '/', kic_list=kic_list, kic_file=kic_file, kic_file_col=kic_file_col)

	print('[4] Computing the Power Spectrum...')
	l=len(valid_kic)
	for i in range(l):
		print(' [' + str(i+1) + '/' + str(l) + '] Processing KIC ', valid_kic[i], '...') 
		do_TF(dir_base_out + '/' + cadence + '/', ID=valid_kic[i]) # Filter by KIC as here there might be several KIC in a single directory
	return valid_kic

print('tf_from_kiclist.py, version ', version(), ' loaded')