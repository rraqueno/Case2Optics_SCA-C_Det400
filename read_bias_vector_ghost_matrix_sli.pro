function read_bias_vector_ghost_matrix_sli, input_file_name

   envi_open_file, input_file_name, r_fid = input_file

   envi_file_query, input_file, bnames=band_names, nb=n_bands, nl=n_lines, $
      ns=n_samples, spec_names = spectra_names

   file_dims=[ -1L, 0, n_samples-1, 0, n_lines-1]
   
   ghost_matrix = envi_get_data( dims=file_dims, fid=input_file, pos=0)

   spectra_names = spectra_names + ';'

   bias = double((stregex(spectra_names,'Bias=([-.0-9]+)[;,]',/subexpr,/extract))[1,*])
   lon = double((stregex(spectra_names,'Lon=([-.0-9]+)[;,]',/subexpr,/extract))[1,*])
   lat = double((stregex(spectra_names,'Lat=([-.0-9]+)[;,]',/subexpr,/extract))[1,*])

   answer = {bias: bias, ghost_matrix:ghost_matrix, lon:lon, lat:lat}

   return, answer

end
