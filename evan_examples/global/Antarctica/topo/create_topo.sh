#! /bin/bash

# use RTopo (Schaffer et al 2016) doi:10.1594/PANGAEA.856844

# change the path to where you store RTopo (or use another topography grid if you want)
# RTopo is convenient because it already takes off the Greenland and Antarctica ice

# switch this to 'y' the first time you run, because Rtopo has an incorrect header. Note, it also does a filter on the data. The grdfilter step, which #
# runs a 5 km median filter on the topography takes like 2 hours!
# Probably could just run it without doing that, but it is likely good in order to prevent aliasing.
# note, if you don't want to do the filtering step, you have to change the grdproject command to use the variable ${topo} instead of ${filtered_topo}
first_run=n


# switch this to 'n' if you have already generated the reduced 
run_project=n


original_topo=/scratch/users/egowan-local/topo/RTopo/RTopo-2.0.1_30sec_bedrock_topography.nc

topo=bed_topo.nc

margin_file=../margins/20000.gmt

# Rtopo is not formatted correctly as a COARDS compliant netcdf file. This command fixes that. Only need to run this once.

if [ "${first_run}" == "y" ]
then
ncrename -O -d londim,x -d latdim,y -v lon,x -v lat,y ${original_topo} ${topo}

fi

# reading from projection_info.sh now

# For Lambert azimuthal projection

#center_longitude=-94
#center_latitude=60
#resolution=5 # grid resolution, in km!

# corner points of the grid (if we don't use this, gmt assumes a global grid, which will be huge!
# west corresponds to the bottom left corner, east corresponds to the top right corner
# probably easiest to pick off the cordinates off Google Earth, in a really zoomed out view
#west_latitude=25
#west_longitude=-135
#east_latitude=58
#east_longitude=3

#map_width=15c

source ../projection_info.sh

mapproject << END    ${R_options} ${J_options} -F -C  > corners.txt
${west_longitude} ${west_latitude}
${east_longitude} ${east_latitude}
END

spacing=${resolution}000

r1=$(awk '{if (NR==1) print $1}' corners.txt)
r2=$(awk '{if (NR==2) print $1}' corners.txt)
r3=$(awk '{if (NR==1) print $2}' corners.txt)
r4=$(awk '{if (NR==2) print $2}' corners.txt)

# round the numbers, should only need to do this for the top left corner, really

x_min=${r1}
y_min=${r3}
x_max_temp=$(printf '%.0f\n' $(echo "scale=2; ${r2} / ${spacing}" | bc ) )
x_max=$(echo "${x_max_temp} * ${spacing}" | bc)
y_max_temp=$(printf '%.0f\n' $(echo "scale=2; ${r4} / ${spacing}" | bc ) )
y_max=$(echo "${y_max_temp} * ${spacing}" | bc)



makecpt -Cglobe -T-10000/10000 > shades.cpt

plot=topo_plot.ps

region=Antarctica
#filtered_topo=filtered_topo.nc
#cheating
filtered_topo=/scratch/users/egowan-local/topo/RTopo/bed_topo.nc
area_grid=${region}.nc

# uncomment first time you run, but run only once because it takes a couple of hours!
if [ "${first_run}" == "y" ]
then
grdfilter ${topo} -G${filtered_topo} -Fm${resolution} -D4  -V
fi

# takes a lot less time
#-R${x_min}/${x_max}/${y_min}/${y_max}

grdproject ${filtered_topo}  ${R_options} ${J_options_project} -Garea_temp.nc -D${resolution}000= -Fe  -V -C 

grdsample area_temp.nc -G${area_grid}  -R${x_min}/${x_max}/${y_min}/${y_max}  -I${resolution}000=

#grdconvert area_temp.nc -G${area_grid} -R${west_longitude}/${west_latitude}/${east_longitude}/${east_latitude}r

#${west_longitude} ${west_latitude}
#${east_longitude} ${east_latitude}

#grdproject ${filtered_topo}  -R${west_longitude}/${west_latitude}/${east_longitude}/${east_latitude}r ${J_options} -G${area_grid} -D${resolution}000= -Fe  -V -C 

grdproject ${area_grid}  ${R_options} ${J_options} -G${region}_topo_geo.nc  -I -Fe  -V  -C


x_min=$(grdinfo -F ${area_grid} | grep x_min  | awk -F':' '{print int($3)}')
x_max=$(grdinfo -F ${area_grid} | grep x_max  | awk -F':' '{print int($3)}')
y_min=$(grdinfo -F ${area_grid} | grep y_min  | awk -F':' '{print int($3)}')
y_max=$(grdinfo -F ${area_grid} | grep y_max  | awk -F':' '{print int($3)}')

grdimage ${area_grid} -Y12  -R${x_min}/${x_max}/${y_min}/${y_max}  -JX${map_width}/0 -K -P -Cshades.cpt -V -nb > ${plot}

#grdimage ${area_grid} -Y12  -R${x_min}/${x_max}/${y_min}/${y_max}  -Jx1:${map_scale_factor} -K -P -Cshades.cpt -V -nb > ${plot}


psxy ${margin_file}  ${R_options} ${J_options} -K -O -P -V -Wthickest,white >> ${plot}
psxy ${margin_file}  ${R_options} ${J_options} -K -O -P -V -Wthin,blue >> ${plot}

pscoast -Bafg -O -K -J -R -P -Wthin -Di -A5000 >> ${plot}

psscale -X-2 -Y-3.5 -Dx9c/2c/9c/0.5ch -P -O -Bx4000f1000+l"Elevation (m)" --FONT_LABEL=14p -Cshades.cpt -V  >> $plot

# convert to gmt formatted binary file for use in ICESHEET

bin_file="${region}.bin"

grdconvert ${area_grid} ${bin_file}=bf 

echo ${bin_file} > elev_parameters.txt
echo ${x_min} >> elev_parameters.txt
echo ${x_max} >> elev_parameters.txt
echo ${y_min} >> elev_parameters.txt
echo ${y_max} >> elev_parameters.txt
echo ${resolution}000 >> elev_parameters.txt

# create NetCDF file with equivalent water load

makecpt -Crainbow -T0/4000 -I > shades.cpt

grdmath ${area_grid} 0 LT = ocean_mask.nc

# 1025 / 917 = 1.118 (ratio of density of water to density of ice, used in ICESHEET)
grdmath ${area_grid}  ocean_mask.nc MUL -1.118 MUL = ocean_equivalent_ice.nc

plot=ocean_equivalent.ps

grdimage ocean_equivalent_ice.nc ${shift_up}  -R${x_min}/${x_max}/${y_min}/${y_max}  -JX${map_width}/0 -K -P -Cshades.cpt -V -nb > ${plot}

#psxy ${margin_file}  ${R_options} ${J_options} -K -O -P -V -Wthickest,white >> ${plot}

pscoast -Bafg -O -K ${R_options} ${J_options} -P -Wthin -Di -A5000 >> ${plot}

psscale -X-1 -Y-3.5 -Dx9c/2c/9c/0.5ch -P -O -Bx1000f500+l"equivalent ice thickness (m)" -G0/4000 -Cshades.cpt --FONT_LABEL=14p -V  >> $plot
