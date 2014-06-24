function get_entries, x

  answer = x.entries.toarray()

  return, answer

end

function read_detector_maps, filename
;
; Read in all the data into generic array
;
data = read_csv(filename, header=columns )

;
; Based on the generic array, figure out how many detector entries per SCA
;
sca1_indices = where( data.field03 eq 1)
sca2_indices = where( data.field03 eq 2)
sca3_indices = where( data.field03 eq 3)

;
; Based on the generic array, figure out how many pointing location 
; angles per detector
;
sca1_detector_hist = histogram(data.field05[sca1_indices])
sca2_detector_hist = histogram(data.field05[sca2_indices])
sca3_detector_hist = histogram(data.field05[sca3_indices])

;
; Put all three histograms into an array
;
sca_detector_hist = transpose([[sca1_detector_hist],[sca2_detector_hist],[sca3_detector_hist]])

;
; Number of SCAs and number of detectors per SCA
;
n_sca = 3
n_detectors_per_sca = n_elements( sca_detector_hist[0,*] )

;
; Create the data structure for a detector entry
;
entry = {x_degrees:!values.F_NAN(), y_degrees:!values.F_NAN(), $
	level:!values.F_NAN() }

;
; Create the data structure for to hold detector entries
;
sca_detector_data = replicate({n_entries:0, entries:list()}, n_sca, n_detectors_per_sca) 

;
; Go through each SCA and allocate the number of entries based
; on the histogram.
;
for sca= 0, n_sca-1 do begin 
    for i = 0, n_elements( sca_detector_hist[sca,*] ) -1 do begin
	;
	; We extract the number of entries for each detector
	; 
        n_entries = sca_detector_hist[sca,i]

	;
	; We set the number of entries into the data structure and
	; allocate the space to hold the entry information.
	;
        sca_detector_data[sca,i].n_entries = n_entries
	detector_entries =  replicate(entry, n_entries )

	;
	; Figure out where all the entries are for a given detector
	; within an SCA
	;
        indices = where( data.field03 eq sca+1 and data.field05 eq i+1 )

	;
	; Transfer all those entries into a temporary data structure
	;
	for j = 0, n_entries-1 do begin
           detector_entries[j].level = data.field06[indices[j]]
           detector_entries[j].x_degrees = data.field08[indices[j]]
           detector_entries[j].y_degrees = data.field09[indices[j]]
        endfor

	;
	; Attach that data structure into the overall data structure
	; that will be returned as an answer.
	;
        sca_detector_data[sca,i].entries = list(detector_entries)
    endfor
endfor

;
; Return the data structure as our answer
;
return, sca_detector_data

end

;
;
;
pro fly_sca2_detector400_output_matrix_rotated, image_filename, ghost_map_filename, $
	flight_path_filename 


;
; Get basenames from filename
;
image_file_basename = file_basename( image_filename, '.img')

ghost_file_basename = file_basename( ghost_map_filename, '.csv')

flight_file_basename = file_basename( flight_path_filename,'.txt')

output_file_basename = image_file_basename+'-'+ghost_file_basename+'-'+flight_file_basename

  envi_open_file, image_filename, r_fid=input_fid

  envi_file_query, input_fid, dims=dims

  map_info = envi_get_map_info(fid=input_fid)

  n_samples = dims[2]+1
  n_lines = dims[4]+1


;
; Read in satellite path
;
  envi_read_cols,flight_path_filename, L8_info

  L8_lats = L8_info[0,*]
  L8_lons = L8_info[1,*]
  L8_bias = L8_info[2,*]

;
; Get the data from the GOES image
;
  image = envi_get_data( fid=input_fid, dims=dims, pos=0 )

;
; Get the projection information of the GOES image
;
  image_projection = envi_get_projection(fid=input_fid)

;
; Establish a geographic coordinate system to which the samples
; and lines will be converted
;
  point_projection = envi_proj_create(/geographic)

;
; Convert L8 positions to file coordinates
;
;
; Convert image map coordinates to geographic coordinates
;
      envi_convert_projection_coordinates, L8_lons, L8_lats, point_projection, ground_track_xs, ground_track_ys, image_projection

;
; Now we want to convert the ground_track map coordinates into file coordinates
;
; Convert the file sample,line positions to image map coordinates
;
      envi_convert_file_coordinates, input_fid, ground_track_samples, ground_track_lines, ground_track_xs, ground_track_ys

roi_id = envi_create_roi( name='L8_PATH39',ns=n_samples , nl=n_lines )

envi_define_roi, roi_id, /point, xpts =reform(ground_track_samples), ypts=reform(ground_track_lines)


;
; Figure out how many integer samples and lines are in the L8 coordinates
;
histogram = hist_2D( ground_track_samples, ground_track_lines )

;
; Figure out how many are non-zeroes in the sample/line combinations
;
non_zeroes = where( histogram ne 0 )

;
; Change the vector indices to 2D indices
; 
indices = array_indices( histogram, non_zeroes )

;
; Sort the positions
;
sorted_lines = sort(indices[1,*])

;
; Reorder sorted_samples
;
sorted_indices = indices[*,[sorted_lines]]

;
; The following lines will generate the L8 positions, the biases,
; and the GOES positions that they correspond to. This file will 
; be used to generate a high resolution inversion matrix to go with
; the L8 resolution biases.
;
n_ground_track_samples = n_elements( ground_track_samples )

openw,1,'L8_bias_full_resolution.txt'
for i = 0, n_ground_track_samples -1 do begin

index = where( sorted_indices[0,*] eq fix(ground_track_samples[i]) and sorted_indices[1,*] eq fix(ground_track_lines[i]))

;print, 'i=',i ,'index = ',index, sorted_indices[*,index]

printf,1, l8_lats[i], l8_lons[i], l8_bias[i], index

endfor 

close,1


;
; The following lines will take the L8 positions and sample it down 
; to a single L8 position/bias value to a given GOES pixel.
; We will try to go to extract the middle value in the population
; that corresponds to a GOES pixel.
; 
n_goes_positions = n_elements( sorted_indices[0,*] ) 
biases = fltarr( n_goes_positions )

openw,2,'L8_bias_goes_resolution.txt'
goes_lats = !NULL
goes_lons = !NULL
for i=0, n_goes_positions-1 do begin
  l8_positions = where( sorted_indices[0,i] eq fix( ground_track_samples ) $
	and sorted_indices[1,i] eq fix( ground_track_lines), n_l8_positions )

  index = l8_positions[ n_l8_positions/2] 

  printf,2, l8_lats[ index ], l8_lons[ index ], l8_bias[index], i

  goes_lats = [goes_lats, l8_lats[index]]
  goes_lons = [goes_lons, l8_lons[index]]

endfor
close,2
;
; Figure out the angle to rotate the image based on the 
; first and last entries of the satellite path
;
delta_points = float(sorted_indices[*,-1]-sorted_indices[*,0])

rotate_ghost_angle = atan( delta_points[0]/delta_points[1] ) /!DTOR

goes_roi_id = envi_create_roi( name='GOES_PATH39',color=3, ns=n_samples , nl=n_lines )

envi_define_roi, goes_roi_id, /point, xpts =reform(sorted_indices[0,*]), ypts=reform(sorted_indices[1,*])

window,0,xsize=n_samples,ysize=n_lines
        tvscl,/ord,image
window,1,xsize=n_samples,ysize=n_lines
	wset,1
 tvscl,/ord,image

positions = sorted_indices

n_positions = n_elements( positions[0,*] )

  image_array = fltarr(n_samples,n_lines,n_positions)
  image_array[*]=!VALUES.F_NAN()
  image_array_display = fltarr(n_samples,n_lines,n_positions)

  for i = 00, n_positions-1 do  begin

    image_array_display[*,*, i ] = image

  endfor

  sca_detector_data = read_detector_maps(ghost_map_filename)

  band=0

band_names=strarr(n_positions)

SCA=2-1
detector=400-1
;scan_pos = 444-1

;
; In order to write out the radiance matrix corresponding to the weights matrix,
; we need to allocate the number of positions that the ghost map will sample
; as well as a matrix for the ghost map weights.
;
; In this example, radiance[*,*,0] will contain the ghost map weights and
; radiance[*,*,1:positions] will contain the radiance values sampled by each ghost map
; weight.
;
  radiance_array = fltarr(205, 205, n_positions+4)
  radiance_array[*] = !VALUES.F_NAN()
  ghost_radiance_array = radiance_array
  ground_track_values= fltarr(n_positions)

;
; Populate the first band with all the weights of the image
; after extracting out the weights
;
entry=sca_detector_data[sca,detector].entries.toarray()
n = n_elements(entry) 

for i=0,n-1 do begin

    x=entry[i].x_degrees
    y=-entry[i].y_degrees
    sample = tan(x*!DTOR)*281.14+103
    line = tan(y*!DTOR)*281.14+103
    radiance_array[sample,line,3]=entry[i].level
if ~finite(radiance_array[sample,line,0]) then radiance_array[sample,line,0]=0.0
    radiance_array[sample,line,0]=radiance_array[sample,line,0]+1
if ~finite(radiance_array[sample,line,1]) then radiance_array[sample,line,1]=0.0
    radiance_array[sample,line,1]=radiance_array[sample,line,1]+entry[i].level
if ~finite(radiance_array[sample,line,2]) then radiance_array[sample,line,2]=0.0
    radiance_array[sample,line,2]=radiance_array[sample,line,2]+entry[i].level^2 
    
endfor

;
; 1. Rotate the original ghost map and
; 2. put it back to into the original radiance_array
;
;ghost_radiance_array[*,*,0] = rot(radiance_array[*,*,0], rotate_ghost_angle )
;radiance_array[*,*,0] = ghost_radiance_array[*,*,0]

;ghost_radiance_array[*,*,1] = rot(radiance_array[*,*,1], rotate_ghost_angle )
;radiance_array[*,*,1] = ghost_radiance_array[*,*,1]

;ghost_radiance_array[*,*,2] = rot(radiance_array[*,*,2], rotate_ghost_angle )
;radiance_array[*,*,2] = ghost_radiance_array[*,*,2]

;ghost_radiance_array[*,*,3] = rot(radiance_array[*,*,3], rotate_ghost_angle )
;radiance_array[*,*,3] = ghost_radiance_array[*,*,3]

;
; Figure out where all the values greater than zero exits in the weights band.
; We will use this to index the subsection of the GOES image as the ghost map
; is moved across the scene
;
;gt_zeroes = where(ghost_radiance_array[*,*,3] gt 0.0 )
gt_zeroes = where(radiance_array[*,*,3] gt 0.0 )

;
; Pull out the GOES subsection for the ground track locations
;
for band = 0,n_positions-1 do begin  
    scan_pos = sorted_indices[ 0, band ]
    line_pos = sorted_indices[ 1, band ]
      ground_track_sample = scan_pos
      ground_track_line = line_pos 
      ground_track_values[ band ] = image[ground_track_sample, ground_track_line]
 
      goes_subsection = image[ground_track_sample-102:ground_track_sample+102, ground_track_line-102:ground_track_line+102]

      rotated_goes_subsection = rot(goes_subsection,missing=0.0, -rotate_ghost_angle) 


      radiance_array[*,*,band+4] = rotated_goes_subsection

;
; Fill temp variable with all NaNs
;
      temp = ghost_radiance_array[*,*,band+4] 

      temp[gt_zeroes] = rotated_goes_subsection[gt_zeroes]

      ghost_radiance_array[*,*,band+4] = temp

;
; Convert the file sample,line positions to image map coordinates
;
      envi_convert_file_coordinates, input_fid, ground_track_sample, ground_track_line, ground_track_x, ground_track_y, /to_map

;
; Convert image map coordinates to geographic coordinates
;
      envi_convert_projection_coordinates, ground_track_x, ground_track_y, image_projection, ground_track_lon, ground_track_lat, point_projection
      
        print, 'POS =',band,' Sample =',ground_track_sample, ' Line =', ground_track_line, ' Lat = ', ground_track_lat, ' Lon = ', ground_track_lon

wset,1
        tvscl,/ord,image_array[*,*,band]
        xyouts,10,10,ghost_map_filename+' SCA:'+strtrim(string(sca+1),2)+' '+'Detector:'+strtrim(string(detector+1),2) + " "+ "Entries="+strtrim( string(n),2)+" "+"POS="+strtrim(band,2),/device
	wait,0.125

	band_names[band]='File:'+ghost_map_filename+'; SCA:'+strtrim(string(sca+1),2)+'; '+'Detector:'+strtrim(string(detector+1),2)+"; "+"Entries="+ strtrim(string(n),2)+"; BoresightLat="+strtrim(string(goes_lats[band]),2)+"; BoresightLon="+strtrim(string(goes_lons[band]),2)+"; BoresightLine="+strtrim(string(ground_track_line),2)+"; BoresightSample="+strtrim(string(ground_track_sample),2)

wset,0

        tvscl,/ord,image_array_display[*,*,band]
        xyouts,10,10,ghost_map_filename+' SCA:'+strtrim(string(sca+1),2)+' '+'Detector:'+strtrim(string(detector+1),2) + " "+ strtrim(string(n),2)+" entries",/device
	
endfor

  envi_write_envi_file, radiance_array,out_name=output_file_basename+'-unsampled'+'.img',bnames=["Vector Hits","Sum Weights","Sum^2 Weights","Ghost Weights",band_names], map_info=map_info

  envi_write_envi_file, ghost_radiance_array,out_name=output_file_basename+'-'+'sampled'+'.img',bnames=["Vector Hits","Sum Weights","Sum^2 Weights","Ghost Weights",band_names], map_info=map_info



;
; The following files represent the radiances corresponding to the tracks
; supplied by the user.  This can either be subsatellite or boresight 
; lat lon tracks and the corresponding pixel value is extracted.
;
  openw, lun, output_file_basename+'-ground_track_radiances.txt',/get_lun
  printf,lun,transpose(ground_track_values)
  free_lun, lun

;
; Write out the unsampled and sampled bands in interleaving form for 
; flicker visual inspection
;
flicker_array = !NULL
flicker_band_names = !NULL
for i = 0, n_positions+3 do begin

  flicker_array= [[[flicker_array]],[[radiance_array[*,*,i]]],[[ghost_radiance_array[*,*,i]]]]

endfor 

  envi_write_envi_file, flicker_array,out_name=output_file_basename+'-'+'flicker'+'.img'

end
