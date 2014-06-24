pro extract_modis_pixels_given_l8_detector_lat_lons, l8_filename, modis_filename
;
; Establish base filename for combined input data
;
l8_basename = file_basename(l8_filename,'.img')
modis_basename = file_basename(modis_filename,'.img')
combined_filename = l8_basename+'-' + modis_basename

;
; 1-2. Get the L8 image information
; 3-4. Get L8 image dimension
;
  envi_open_file, l8_filename, r_fid=l8_input_fid
  envi_file_query, l8_input_fid, dims=l8_dims

  l8_n_samples = l8_dims[2]+1
  l8_n_lines = l8_dims[4]+1

;
; 1-3. Get the data from the l8 detector image data
;
  l8_data = envi_get_data( fid=l8_input_fid, dims=l8_dims, pos=0 )
  l8_lats = envi_get_data( fid=l8_input_fid, dims=l8_dims, pos=1 )
  l8_lons = envi_get_data( fid=l8_input_fid, dims=l8_dims, pos=2 )

;
; Get the MODIS image information
;
  envi_open_file, modis_filename, r_fid=modis_input_fid
  envi_file_query, modis_input_fid, dims=modis_dims

  modis_n_samples = modis_dims[2]+1
  modis_n_lines = modis_dims[4]+1

;
; 1. Get the data from the l8 detector image data
; 2. Get the projection information of the MODIS image
  modis_image = envi_get_data( fid=modis_input_fid, dims=modis_dims, pos=0 )
  modis_projection = envi_get_projection(fid=modis_input_fid)

;
; Establish the point coordinates from the L8 image as geographic
;
  point_projection = envi_proj_create(/geographic)

;
; Convert l8 geographic coordinates to modis projection coordinates
;
      envi_convert_projection_coordinates, l8_lons, l8_lats, point_projection, modis_lons, modis_lats, modis_projection

;
; Convert modis projection coordinate to image file coordinates
;
      envi_convert_file_coordinates, modis_input_fid, modis_samples, modis_lines, modis_lons, modis_lats


;
; 1. Extract out modis data corresponding to the locations of the l8 data
; 2. Create bias map by subtracting modis radiance minus landsat radiance
;
modis_data = modis_image[ [modis_samples ],[modis_lines] ]
bias_map = modis_data - l8_data

;
; 1. Create ROI representing the L8 locations on the MODIS image.
; 2. Populate ROI points from L8/MODIS lat lons
; 3. Save out ROI for future reference
;
roi_id = envi_create_roi( name=file_basename(l8_filename,'.img'),ns=modis_n_samples , nl=modis_n_lines )

envi_define_roi, roi_id, /point, xpts =reform(modis_samples), ypts=reform(modis_lines)

envi_save_rois, combined_filename+'.roi', roi_id

;
; 1. Tack on the modis data to the L8 radiance, lat, lon image 
; 2. Write it out to an envi file
;
combined_data=[[[l8_data]],[[modis_data]],[[bias_map]],[[l8_lats]],[[l8_lons]]]

;
; Write out the lat lon bias data text file for reference
;
full_set = transpose([[[l8_lats]],[[l8_lons]],[[bias_map]]])
openw,lun,combined_filename+'-geo-bias'+'.txt',/get_lun
;printf,lun, full_set
free_lun,lun

;
; Figure out the lines that are not in the BLL file
;
envi_read_cols,l8_basename+'.bll',bad_lines_list

water_mask=fltarr(n_elements(full_set[0,*]))
water_mask[*] = -1
water_mask[bad_lines_list-1]=bad_lines_list-1
water_lines = where( water_mask eq -1 )

water_set = full_set[*,water_lines]

openw,lun,combined_filename+'-geo-bias-water'+'.txt',/get_lun
printf,lun, water_set
free_lun,lun

;print,n_elements(water_set)

path_midpoint = water_set[*,n_elements(water_set[0,*])/2-1]

openw,lun,combined_filename+'-path-midpoint'+'.txt',/get_lun
printf,lun, path_midpoint
free_lun,lun

;
; Write out the augmented ENVI file
;
envi_write_envi_file, combined_data, out_name=combined_filename+'.img',bnames=['L8 radiance','MODIS radiance','MODIS-L8 radiance','latitude','longitude']


end
