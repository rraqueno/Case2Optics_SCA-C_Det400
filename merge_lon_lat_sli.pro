pro merge_lon_lat_sli, lon_file_name, lat_file_name

   envi_open_file, lon_file_name, r_fid = lon_input_file
   envi_open_file, lat_file_name, r_fid = lat_input_file

   envi_file_query, lon_input_file, bnames=band_names, nb=n_bands, nl=n_lines, $
      ns=n_samples, spec_names = lon_spectra_names

   envi_file_query, lat_input_file, bnames=band_names, nb=n_bands, nl=n_lines, $
      ns=n_samples, spec_names = lat_spectra_names

   file_dims=[ -1L, 0, n_samples-1, 0, n_lines-1]
   
   lon_matrix = envi_get_data( dims=file_dims, fid=lon_input_file, pos=0)
   lat_matrix = envi_get_data( dims=file_dims, fid=lat_input_file, pos=0)

   lon_lat_matrix = !NULL
   lon_lat_spec_names = !NULL

   for i =0, n_lines - 1 do begin

     lon_lat_matrix = [ [lon_lat_matrix], [lon_matrix[*,i]],$
	[lat_matrix[*,i]]] 

     lon_lat_spec_names = [ lon_lat_spec_names, lon_spectra_names[i], $
	lat_spectra_names[i] ]


   endfor

file_type="ENVI Spectral Library"


envi_write_envi_file, lon_lat_matrix, file_type=envi_file_type(file_type), out_name=lon_file_name+'_'+lat_file_name, interleave=0, ns=n_samples, nl=2*n_lines, nb=1, spec_names = lon_lat_spec_names, bnames=[band_names]



end
