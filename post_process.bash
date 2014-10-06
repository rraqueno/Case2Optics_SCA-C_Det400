pwd_dir=`pwd`
base=`basename $pwd_dir`
cp -uv ghost_weights.sli $base-ghost_weights.sli
cp -uv ghost_weights.hdr $base-ghost_weights.hdr
cp -uv GOES_ghost_matrix.sli $base-GOES_ghost_matrix.sli
cp -uv GOES_ghost_matrix.hdr $base-GOES_ghost_matrix.hdr
cp -uv L8_ghost_matrix.sli $base-L8_ghost_matrix.sli
cp -uv L8_ghost_matrix.hdr $base-L8_ghost_matrix.hdr
