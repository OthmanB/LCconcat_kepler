pro doTF, ftime, fflux, ofac, sav_name
  compile_opt idl2

  ; read a 1 column file with time
  time = reader(ftime)
  ; read a 1 column file with flux
  flux = reader(fflux)
  print, 'Computing the TF_lnp...'
  spec = TF_lnp(time, flux, ofac = ofac, sav_name = sav_name)
end

function reader, file
  compile_opt idl2
  openr, unit, file, /get_lun
  Nvalmax = 2000000
  values = dblarr(Nvalmax)

  ; Read the file line by line
  i = 0
  while not eof(unit) do begin
    line = ''
    readf, unit, line
    values[i] = double(line)
    i = i + 1
  endwhile
  close, unit
  values = values[0 : i - 1]
  return, values
end

; *******************************************************************
; fonction de creation d'un spectre ‡ partir d'une serie irreguliere
function TF_lnp, time, flux, ofac = ofac, sav_name = sav_name
  compile_opt idl2

  if n_elements(ofac) eq 0 then ofac = 1
  if n_elements(sav_name) eq 0 then sav_name = ''

  len = n_elements(flux)
  t0 = dblarr(len)
  f0 = dblarr(len)
  t0[*] = time[*] - min(time[*])
  f0[*] = flux[*]
  ptot = double(0)
  freq = double(0)
  res_lomb = lnp_test(t0 * 86400.0d+00, f0, $
    ofac = ofac, wk1 = freq, wk2 = ptot, /double)

  tot_MS = total((flux[0 : len - 1] - mean(flux[0 : len - 1])) ^ 2) / float(len)
  tot_lomb = total(ptot[0 : len / 2 - 1])

  freq = freq * 1.0d+06

  bw = freq[1] - freq[0]
  ptot = ptot * (tot_MS / tot_lomb) / bw

  result = dblarr(2, n_elements(freq))
  result[0, *] = freq
  result[1, *] = ptot

  if sav_name ne '' then begin
    ; freq=power_spec[0,*]
    spec_reg = ptot
    save, freq, spec_reg, filename = sav_name
  endif

  return, result
end