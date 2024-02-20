!	copyright Evan. J. Gowan, 2013, 2016
!
!	This file is part of ICESHEET 1.0
!
!	ICESHEET is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, version 3 of the License
!
!	ICESHEET is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with ICESHEET.  If not, see <http://www.gnu.org/licenses/>.

program icesheet
    ! last change: January 11, 2016

    ! Written by Evan J. Gowan (evangowan@gmail.com)


    ! This program had its beginnings in early 2011, and became a major part of my PHD project at the Australian National University
    ! working on the Laurentide Ice Sheet. I had to overcome a lot of hurdles, such as my inexperience with programming, and
    ! some major problems related to detecting when the individual flowlines overlapped. I finally got a fairly
    ! stable version of the program working in August 2013, and that was used for completing my PHD project.
    ! After the finishing the corrections of my thesis (May 2015), I went back to this program to get it in a nice
    ! enough shape to publish. I fixed several bugs in it, including a bug in the ice thickness at the margin
    ! if it is below sea level, and some problems related to places where the ice sheet is extremely flat.

    ! I hope that this program proves useful to make quick estimations of ice sheet geometry through a glacial cycle.

    ! if you use this program, cite:

    ! ICESHEET 1.0: A program to produce paleo-ice sheet models with minimal assumptions
    ! Evan J. Gowan, Paul Tregoning, Anthony Purcell, James Lea, Oscar J. Fransner, Riko Noormets, and J. A. Dowdeswell
    ! submitted to Geoscience Model Development

    use global_parameters
    use grids

    implicit none

    character(len=10) :: elevation_interval_char
    character(len=10) :: minimum_spacing_char
    character(len=255) :: margin_file
    character(len=255) :: output_directory
    character(len=10) :: output_suffix

    integer :: counter1, counter2, counter3, counter4, istat, boundary_counter
    double precision :: grid_x, grid_y, latitude, longitude, x, y

    call GET_COMMAND_ARGUMENT(1, elevation_interval_char)
    call GET_COMMAND_ARGUMENT(2, minimum_spacing_char)
    call GET_COMMAND_ARGUMENT(3, margin_file)
    call GET_COMMAND_ARGUMENT(4, ss_netcdf_file)
    call GET_COMMAND_ARGUMENT(5, elevation_netcdf_file)
    call GET_COMMAND_ARGUMENT(6, output_directory)
    call GET_COMMAND_ARGUMENT(7, output_suffix)

    read(elevation_interval_char,*) elevation_interval
    read(minimum_spacing_char,*) minimum_spacing

    write(6,*) "Margin file:       ", margin_file
    write(6,*) "Shear stress file: ", ss_netcdf_file
    write(6,*) "Elevation file:    ", elevation_netcdf_file
    
    write(6,*) "Reading elevation netcdf"
    call load_elevation()

    write(6,*) "Reading shear stress netcdf"
    call load_shear_stress()

    ! subroutine that reads in the ice file, converts the boundary to distance from the central point, and find the local slope
    write(6,*) "Reading ice margin file"
    call read_icefile(margin_file)

    write(*,*) "Input file read successful, beginning calculations"
    call find_flowline(output_directory, output_suffix)
end program icesheet
