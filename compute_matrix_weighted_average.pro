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


function read_ghost_weights_sli, input_file_name

   envi_open_file, input_file_name, r_fid = input_file

   envi_file_query, input_file, bnames=band_names, nb=n_bands, nl=n_lines, $
      ns=n_samples, spec_names = spectra_names

   file_dims=[ -1L, 0, n_samples-1, 0, n_lines-1]
   
   ghost_weights = envi_get_data( dims=file_dims, fid=input_file, pos=0)

   return,  ghost_weights 

end

pro compute_matrix_weighted_average, vector_matrix_sli_file, weights_sli_file

    answer_filename = file_basename( vector_matrix_sli_file ) + '.csv'
    vector_matrix = read_bias_vector_ghost_matrix_sli( vector_matrix_sli_file )
    bias = reform(vector_matrix.bias)
    radiance_matrix = vector_matrix.ghost_matrix
    lon = reform(vector_matrix.lon)
    lat = reform(vector_matrix.lat)

    weights = read_ghost_weights_sli( weights_sli_file )

    weighted_average = weights##transpose(radiance_matrix)/total(weights)

    answer = transpose([[lon],[lat],[bias],[weighted_average]])

    write_csv, answer_filename, answer, header=['Lon','Lat','Bias','Weighted Average']


end
