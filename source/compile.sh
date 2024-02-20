gfortran -o global_parameters.o -O2 -c global_parameters.f90
gfortran -o grids.o -O2 -c grids.f90 $(nc-config --flibs)
gfortran -o read_icefile.o -O2 -c read_icefile.f90
gfortran -o flowline_location.o -O2 -c flowline_location.f90
gfortran -o find_flowline_fisher_adaptive_4.o -O2 -c find_flowline_fisher_adaptive_4.f90
gfortran -o icesheet -O2 icesheet.f90 global_parameters.o grids.o read_icefile.o find_flowline_fisher_adaptive_4.o flowline_location.o $(nc-config --flibs) 
