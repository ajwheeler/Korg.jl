DIR=~/Dropbox/departure_coefficients/amarsi

for el in Al Ba C Ca Fe H K Li Mg Na O Si Ti
do
    julia pack_amarsi_coefficients.jl $DIR/$el $el.h5
done
