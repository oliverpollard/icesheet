#! /bin/bash

# The easiest way to create a shear stress file is to partition the area of interest into places where you want to be able to vary the shear stress.

# The way I have done this is to create a shapefile that partitions the area as desired. I then exported the shapefile into a CSV file.
# from QGIS:
# > Open shapefile
# > Layer dropdown menu > save as
# > in the options, change "Symbology Export" to "feature symbology", and "Geometry" as "AS_WKT".
# > Change the CRS to WGS 84, or else it will not export the points as latitude/longitude, if that is desired

binary_path="."

csv_file=polygons.csv

is_compiled=$(command -v ${binary_path}/convert_grid)


if [ -z "${is_compiled}" ];
then
	binary_path="."
	make convert_grid
	make create_ss_grid

fi

if [ ! -f ${csv_file} ]; then
    echo "File 'polygons.csv' does not exist, please go into the shear_stress folder"
    echo "and follow the instructions in the readme file to create it"
    exit 0
fi


awk 'BEGIN { FPAT = "([^,]+)|(\"[^\"]+\")"} {if(NR > 1) {print ">", NR-1; print $1} }' ${csv_file} | sed -e 's/\"MULTIPOLYGON (((//g'  | sed -e 's/)))\"//g' | sed -e 's/,/\n/g' | sed -e 's/(//g' | sed -e 's/)//g' |  sed -e 's/\"POLYGON //g'  |  sed -e 's/\"//g' > gmt_file.txt

# if you are running on a MAC, use this command instead of the one above. Thanks to Felicity Williams for pointing this out.
#awk 'BEGIN { FPAT = "([^,]+)|(\"[^\"]+\")"} {if(NR > 1) {print ">", NR-1; print $1} }' ${csv_file} | sed -e 's/\"MULTIPOLYGON (((//g'  | sed -e 's/)))\"//g' | sed -e 's/,/ \'$'\n/g' | sed -e 's/(//g' | sed -e 's/)//g' |  sed -e 's/\"POLYGON //g'  |  sed -e 's/\"//g' > gmt_file.txt

awk 'BEGIN { FPAT = "([^,]+)|(\"[^\"]+\")"} {if(NR > 1) {print NR-1, $2} }' ${csv_file}  > shear_stress_domain_values.txt



lat_spacing=5000
long_spacing=5000

${binary_path}/create_ss_grid ${lat_spacing} ${long_spacing}

${binary_path}/convert_grid shear_stress_domain_values.txt
