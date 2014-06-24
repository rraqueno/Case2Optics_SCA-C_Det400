pro write_ghost_weights_radiances_as_matrix, input_image_filename


  input_image_file_basename = file_basename( input_image_filename, '.img')

  envi_open_file, input_image_filename, r_fid=input_fid
  envi_file_query, input_fid, dims=dims, bnames=band_names

  n_bands = n_elements(band_names)

  weight_data = envi_get_data( dims=dims, fid=input_fid, pos=0)

  gt_zero = where( weight_data gt 0.0, n_contributors )

  matrix = dblarr( n_bands, n_contributors )

  matrix[ 0, * ] = weight_data[ gt_zero ]

  header = ["SCA=2; Detector=400;"]

  for i = 1, n_bands - 1 do begin

   header = [header, strmid(band_names[i],48,48) ]

    radiance_data = envi_get_data( dims=dims, fid=input_fid, pos=i)

    matrix[i,*] = radiance_data[ gt_zero ]

  endfor

;
; Write out the data as a CSV file
;
  write_csv, input_image_file_basename+'-Ghost_Weights_Radiances.csv', matrix, header=header 

;
; Write out the data as a text file
;
  openw,  lun,  input_image_file_basename+'-Ghost_Weights_Radiances.txt',/get_lun

printf, lun, matrix

free_lun,lun

help,matrix
 
end
