ghost_weights.sli GOES_ghost_matrix.sli L8_ghost_matrix.sli : MODIS.img L8.img L8.bll B11.csv
	echo "envi & STEP1 & STEP2 & STEP3 " | /cygdrive/c/Program\ Files/Exelis/IDL82/bin/bin.x86_64/idl.exe
	./post_process.bash 

MODIS.img L8.img L8.bll B11.csv : 
	./pre_process.bash

clean:
	rm -f ghost_weights.sli GOES_ghost_matrix.sli L8_ghost_matrix.sli MODIS.img MODIS.img.hdr L8.img L8.hdr L8.bll B11.csv
