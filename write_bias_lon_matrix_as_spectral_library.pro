pro write_bias_lon_matrix_as_spectral_library, input_image_filename, detector_positions_filename


  input_image_file_basename = file_basename( input_image_filename, '.img')

  envi_open_file, input_image_filename, r_fid=input_fid
  envi_file_query, input_fid, dims=dims, bnames=band_names, nb=n_bands, fname=fname
  case_name = file_basename(file_dirname(fname))

;
; The following magic number 3 represents the 3 extra bands we needed to
; include to handle the case where we have more than one LOS vector
; hitting a GOES pixel.  The bands are N, Sum, Sum^2
;
  weight_data = envi_get_data( dims=dims, fid=input_fid, pos=3)

  histogram_band = envi_get_data( dims=dims, fid=input_fid, pos=0)

  histogram_band_int = fix( histogram_band )

  histogram_data = histogram( histogram_band_int )

  n_values = n_elements( histogram_data )



;
; Compute total number of non-zero points in the array
;
  n_vector_contributors = transpose( histogram_data ) # indgen( n_values )


  ;gt_zero = where( weight_data gt 0.0, n_pixel_contributors )

  matrix = dblarr( n_bands-4, n_vector_contributors )

  ;matrix[ 0, * ] = weight_data[ gt_zero ]

;  header = ["SCA=2; Detector=400;"]

  header = !NULL
  weights = !NULL

    for j = 1, n_elements( histogram_data ) - 1 do begin

       select_points  = where( histogram_band_int eq j )

       for k = 0, j-1 do begin
           weights = [ weights, weight_data[ select_points  ] ]
       endfor 
    endfor

;matrix[0,*]=weights

  for i = 4, n_bands  - 1 do begin

   header = [header, 'TYPE:LON; '+band_names[i] ]

    radiance_data = envi_get_data( dims=dims, fid=input_fid, pos=i)

    contributor = !NULL

    for j = 1, n_elements( histogram_data ) - 1 do begin

       select_points  = where( histogram_band_int eq j )

       for k = 0, j-1 do begin
           contributor = [ contributor, radiance_data[ select_points  ] ]
       endfor 

    endfor

      matrix[i-4,*] = contributor

  endfor

;
; Write out the data as a CSV file
;
  ;write_csv, input_image_file_basename+'-Ghost_Weights_Radiances.csv', matrix, header=header 

;
; Write out the data as a text file
;
  openw,  lun,  input_image_file_basename+'-Ghost_Weights_Radiances.txt',/get_lun

printf, lun, matrix

free_lun,lun

help, matrix

transposed_matrix = transpose(matrix)

n_samples = (size(transposed_matrix))[1]
n_lines = (size(transposed_matrix))[2]

file_type="ENVI Spectral Library"

help,weights

;envi_write_envi_file, weights, file_type=envi_file_type(file_type), out_name='L8-weight_data.sli', interleave=0, ns=n_samples, nl=1, nb=1, spec_names = 'TYPE:LON;'+[case_name], bname='TYPE:LON;'+[case_name]


envi_write_envi_file, transposed_matrix, file_type=envi_file_type(file_type), out_name='GOES-LON_data.sli', interleave=0, ns=n_samples, nl=n_lines, nb=1, spec_names = header, bname='TYPE:LON; '+[case_name]

;
; Create a matrix that is higher resolution using the L8 positions and 
; bias values.
;
envi_read_cols, detector_positions_filename, l8_data


n_entries = n_elements(l8_data[0,*])

;
; We are adding one additional entry to hold the ghost weights
; in the first column.
; Also adding an unused 0 element in the beginning of band_index
; to match the additional entry of the ghost weights.
;
band_index = [0,reform(l8_data[3,*])]
l8_matrix = dblarr( n_entries, n_vector_contributors )

;l8_header = ["SCA=2; Detector=400;"]

;l8_matrix[0,*] = matrix[0,*]

l8_header = !NULL

for j = 0, n_entries-1 do begin

   i = band_index[j]
   l8_matrix[ j, * ] = matrix[ i, *] 
   l8_header = [l8_header, header[i] ]

endfor

;
; Write out the data as a CSV
;
  ;write_csv, input_image_file_basename+'-high_res_Ghost_Weights_Radiances.csv', l8_matrix, header=l8_header 
 
;
; Write out the data as a text file
;
  openw,  lun,  input_image_file_basename+'-high_res_Ghost_Weights_Radiances.txt',/get_lun

printf, lun, l8_matrix

free_lun,lun

help, l8_matrix

transposed_l8_matrix = transpose(l8_matrix)
n_samples = (size(transposed_l8_matrix))[1]
n_lines = (size(transposed_l8_matrix))[2]

envi_write_envi_file, transposed_l8_matrix, file_type=envi_file_type(file_type), out_name='L8-LON_data.sli', interleave=0, ns=n_samples, nl=n_lines, nb=1, spec_names = l8_header, bname=['TYPE:LON; '+case_name]

end
