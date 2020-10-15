import subprocess
import os, fnmatch

def version():
	return 'v0.1'

def npzTF2idlsav(dir_npz, npzfile, dir_sav=''):
	'''
		A small program that calls idl in order to save the npz file with the spectrum  
		into a sav file.
		Because the idlpy bridge is not always available (no pip install)
		here we pass by ascii temporary file
		Typically, the conversion is as follow:
			1. Read the npz file
			2. Write the npz content temporary files: Support only the spectrum content
			3. Call an IDL subroutine inside python that will save the file in sav format 
	'''
	if dir_sav == '':
		dir_sav=dir_npz

	current_dir=os.getcwd()
	ascii_file=npzTF2ascii(dir_npz, npzfile, dir_ascii=dir_sav)
	subprocess.call(["mv", ascii_file, "tmp_in.ascii"])
	subprocess.call(["idl", "npz2sav_script.pro"])  # WARNING
	file_list=get_files_list(current_dir, extension='_out.sav', prefix='')
	print('file_list=', file_list)
	if file_list[0] == '':
		print('Error: Could not find the temporary output sav file')
		#exit()
	else:
		# If everything went well, we can move the temporary sav file into 
		# the directory dedicated for this sav file
		savfile=npzfile.split('.')[0] + '.sav'
		subprocess.call(["mv", file_list[0], dir_sav+savfile])

def npzTF2ascii(dir_npz, npzfile, dir_ascii=''):
	
	if dir_ascii == '':
		dir_ascii=dir_npz

	core_file=npzfile.split('.')[0] # remove the extension

	d=np.load(dir_npz+npzfile)
	freq=d['freq']
	power=d['power']
	id_number=str(d['id_number'])
	infos=str(d['infos'])
	#freq=freq, power=power, id_number=id_number, config=config, infos=infos
	ascii_file=dir_ascii + core_file + '.ascii'
	print(ascii_file)
	try:
		f=open(ascii_file, 'w+')
		f.write('#' + str(infos) + '\n')
		f.write('id: ' + str(id_number) + '\n')
		for i in range(len(freq)):
			string='{:f}   {:f}'.format(freq[i], power[i])
			f.write(string + '\n')
		f.close()
		print('Content of ', npzfile, ' written into', ascii_file)
	except:
		print('Error: prepare_lc_kepler.py::npzTF2ascii(): Failed to write npz spectrum content into ascii file')
		print('The program will exit now')
		#exit()
	return ascii_file
