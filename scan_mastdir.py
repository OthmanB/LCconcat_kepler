'''
---------------
Here it is an ensemble of functions that allows you to scan the main directory that 
contains all downloaded data from the MAST website and for Kepler. 
The two main functions are scandir_mast_kepler_bylist() and scandir_mast_kepler_byfile().
The first one allows you to retrieve the locations of the Kepler Quarters from a list containing kic numbers
The second one allows you to retrieve the locations of the Kepler Quarters from a file containing kic numbers
---------------
'''
import os
from prepare_lc_kepler import format_ID

def version():
	return 'v0.1'

def scandir_mast_kepler_bylist(dir_base, kic_list, cadence='lc'):
	'''
		This function looks into the structured mast kepler data (as they were downloaded on the MAST website)
		to all matches of the kic_list list and the selected cadence ('lc' or 'sc')
		It returns full paths for each kic_list values (list of same size as kic_list)
		dir_base: The root directory containing the data. It is expected that it has the 
		following substructure:
				lc ---  0001 --- 000111111
				             --- 000111112
				             --- ....
				             --- [KIC number] (9 digits)
				   ---  0002
				   			--- 000211111
							--- 000211112
							--- ....
							--- [KIC number] (9 digits)
				   ---- ....
				   ---- [4 first digits of kics (9 digits in total)]
				sc ---  0001  --- 000111111
							  --- 000111112
						 	  --- ....
						      --- [KIC number] (9 digits)
				   ---  0002 --- 000211111
							 --- 000211112
							 --- ....
							 --- [KIC number] (9 digits)
				   ---  ....
				   ---- [4 first digits of kics (9 digits in total)]
		kic_list: a simple list of kic number that you are looking for
		cadence: if 'lc', will look for long cadence. if 'sc' will look for short cadence
	'''
	if cadence != 'lc' and cadence != 'sc':
		print('Wrong cadence argument.')
		print('The cadence must be either "lc" for long cadence or "sc" for short cadene')
		exit()

	dir0=dir_base + cadence + '/Kepler/'
	missing_kic=[]
	found_kic=[]
	paths=[]
	for k in kic_list:
		kic=format_ID(k, Ndigits=9) # Ensure that we have a proper encoding into 9 digit of the kic number
		dir_search=dir0 + kic[0:4] # Find the subdir that regroups all kic using the first 4 digits
		subdirs=get_dirs_list(dir_search)
		matching=[]
		for s in subdirs: # Finding the subdir name matching the kic_list element k
			if s.split('/')[-1] == kic:
				matching.append(s)
		if len(matching) == 0:
			print('Warning: No directory found to be matching the kic!')
			print('         Check that data were properly downloaded')
			print('         Pursuing removing the missing kic from the retrieved list')
			missing_kic.append(kic)
		if len(matching) >1:
			print('Warning: multiple directory found matching the kic!')
			print('         Debug might be required')
			print('         Pursuing taking only the first element of the list', matching[0])
			for m in matching:
				print(m)
			found_kic.append(kic)
			paths.append(matching[0])
		if len(matching) == 1:
			found_kic.append(kic)
			paths.append(matching[0])
	if len(missing_kic) !=0 :
		print('List of missing KIC:')
		print(missing_kic)
	else:
		print('All requested paths to kic found')
	return  paths, found_kic, missing_kic

def scandir_mast_kepler_byfile(dir_base, kics_file, kics_col=0, cadence='lc'):
	'''
		Same as scandir_mast_kepler_bylist() but take a text file as input to get a list of kic numbers
		dir_base: The root directory containing the data. See scandir_mast_kepler_bylist() for more details
		kics_file: File with all the kic numbers. It may contain several columns. 
				   Header and comments are allowed and must be marked using a # as first caracter
		kics_col:  Defines the column that must contain the kic numbers. Default is 0.
		cadence:  if 'lc', will look for long cadence. if 'sc' will look for short cadence
	'''
	# Read the file to get a list of kic numbers
	kic_list=[]
	with open(kics_file, 'r') as fh:
		# Skip initial comments that starts with #
		for line in fh:
			if not is_comment(line):
				kic_list.append(line.split()[kics_col]) # If it is not a commented line, split it and take the value in the kics_col coloun
	#print(kic_list)
	# Run the program that scans dir_base to find where the data are
	paths, found_kic, missing_kic=scandir_mast_kepler_bylist(dir_base, kic_list, cadence=cadence)
	return paths, found_kic, missing_kic

def is_comment(line, comment_marker='#'):
	'''
		function to check if a line
	    starts with some character.
	    Here # for comment
	'''
	return line.startswith(comment_marker)

def test_scan_mast(byfile=False):
	'''
		Test function for functions scandir_mast_kepler_bylist() and scandir_mast_kepler_byfile()
	'''
	#home = os.path.expanduser("~")
	#dir_base=home+'/Kepler/lightcurves/'
	dir_base='/Volumes/home/Kepler/lightcurves/'
	cadence='lc' # 'sc'
	if byfile == False:
		kic_list=['2141436', '2715975','3426673']
		paths, found_kic, missing_kic=scandir_mast_kepler_bylist(dir_base, kic_list, cadence='lc')
	else:
		current_dir=os.getcwd()
		kics_file=current_dir + '/inputs/kic_lists.txt'
		paths, found_kic, missing_kic=scandir_mast_kepler_byfile(dir_base, kics_file, kics_col=0, cadence='lc')
	print(paths)

def get_dirs_list(rootdir):
	'''
		A tiny function that scan a directory in order to find all of 
		its sub directories
	'''
	#print('os.scandir(rootdir)=', os.scandir(rootdir))
	dirs = [f.path for f in os.scandir(rootdir) if f.is_dir()]
	return dirs

print('scan_mastdir.py, version ', version(), ' loaded')
#test_scan_mast(byfile=True)