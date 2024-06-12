import numpy as np
from subprocess import call
import os
from scipy.io import readsav

def tf_lnp_idl(time, flux, ofac=1):
    '''
    Calls the IDL tf_lnp function to compute the Lomb-Scargle periodogram.
    Parameters:
        time (array-like): Array of time values.
        flux (array-like): Array of flux values.
        ofac (float): Oversampling factor.
    Returns:
        freq (array-like): Array of periods (in days).
        power (array-like): Lomb-Scargle power corresponding to each period.
    '''
    idl_prg='idl/'
    savfile=os.path.join(idl_prg, 'tf_lnp.sav')
    ftime=os.path.join(idl_prg,'time.txt')
    fflux=os.path.join(idl_prg,'flux.txt')
    np.savetxt(ftime, time)
    np.savetxt(fflux, flux)

    # Get the environment variables
    env = os.environ.copy()
    idl_path = env["IDL_DIR"]
    env["PATH"] = f"{idl_path}:{env['PATH']}"
    # Create the IDL command to execute the procedure
    idl_command = f"""
        .compile -v {idl_prg}/dotf.pro
        .run {idl_prg}/dotf.pro
        dotf, '{ftime}', '{fflux}','{ofac}','{savfile}'
        exit
        """
    # Write the command to a temporary IDL script file
    with open(f"{idl_prg}/tf_lnp.pro", "w") as f:
        f.write(idl_command)
    # Call IDL as a subprocess and execute the script
    call([f"{idl_path}/bin/idl", f"{idl_prg}/tf_lnp.pro"], env=env)
  
    data=readsav(savfile)
    freq=data['freq']
    power=data['spec_reg']

    os.remove(ftime)
    os.remove(fflux)
    os.remove(savfile)

    return freq, power
