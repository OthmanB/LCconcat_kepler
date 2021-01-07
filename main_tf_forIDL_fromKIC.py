'''
	A procedure that handle the generation of the power spectrum
	It creates a copy of the file into a sav format for compatibility
	with my IDL codes
'''
from make_lstf import *
from data_converter import npzTF2idlsav
from prepare_lc_kepler import get_files_list
from scan_mastdir import scandir_mast_kepler_byfile

rootdir_quarters='/Volumes/home/Kepler/lightcurves/'
#dir_out='/Volumes/home/2020/ML-Siddarth/spec_compare/'
#dir_out='/Volumes/home/2020/RGB-SG-Processed/data/'
dir_out='/Volumes/home/Kepler/products/lc/'
cadence='sc'

# ---- TURN ON IF YOU WANT TO COMPUTE AND CONVERT -----
#kic_file='/Volumes/home/2020/ML-Siddarth/mails-files/25-Sept-2020/details_of_mispredictions.txt'
#kic_file_col=0
#valid_kic=make_ALL(rootdir_quarters, cadence, dir_out, kic_list='', kic_file=kic_file, kic_file_col=kic_file_col)

# ---- TURN ON IF YOU WANT TO CONVERT ONLY ----
kic_file='/Volumes/home/2020/RGB-SG-Processed/todo_Siddarth_lc.txt'
kic_file_col=0
paths, found_kic, missing_kic=scandir_mast_kepler_byfile(rootdir_quarters, kic_file, kics_col=0, cadence='lc')
valid_kic=found_kic

#valid_kic=['005689820', '005689820', '007799349', '007799349', '008702606', '008702606', '008751420', '008751420', '009574283', '009574283', '012508433', '012508433']
for k in valid_kic:
	all_files_version=get_files_list(dir_out, extension='_TF.npz', prefix=k)
	for v in all_files_version: # Because we get serveral lightcurves and TF for each star (with different concatenation configuration), this loop is necessary
		print(' Converting ', v , ' ...')
		npzTF2idlsav(dir_out, v, dir_sav=dir_out)
