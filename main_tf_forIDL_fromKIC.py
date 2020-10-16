'''
	A procedure that handle the generation of the power spectrum
	It creates a copy of the file into a sav format for compatibility
	with my IDL codes
'''
from make_lctf import *
from data_converter import npzTF2idlsav
from prepare_lc_kepler import get_files_list

rootdir_quarters='/Volumes/home/Kepler/lightcurves/'
dir_out='/Volumes/home/2020/ML-Siddarth/spec_compare/'
cadence='lc'
kic_file='/Volumes/home/2020/ML-Siddarth/mails-files/25-Sept-2020/details_of_mispredictions.txt'
kic_file_col=0
valid_kic=make_ALL(rootdir_quarters, cadence, dir_out, kic_list='', kic_file=kic_file, kic_file_col=kic_file_col)

for k in valid_kic:
	all_files_version=get_files_list(dir_out, extension='_TF.npz', prefix=k)
	for v in all_files_version: # Because we get serveral lightcurves and TF for each star (with different concatenation configuration), this loop is necessary
		print(' Converting ', v , ' ...')
		npzTF2idlsav(dir_out, v, dir_sav=dir_out)
