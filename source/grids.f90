module grids
    implicit none

    ! name of the netcdf files to be read.
    character (len = 255) :: elevation_netcdf_file, ss_netcdf_file

    integer :: elev_grid_spacing, ss_grid_spacing
    double precision, allocatable, dimension(:) :: elev_x_coords, elev_y_coords, ss_x_coords, ss_y_coords
    double precision, allocatable, dimension(:,:) :: elev_grid, ss_grid

    integer, dimension(6) :: elev_local_grid_x_index, elev_local_grid_y_index
    integer, dimension(6) :: ss_local_grid_x_index, ss_local_grid_y_index
    double precision, dimension(16) :: elev_alpha_array, ss_alpha_array

    ! bicubic interpolation paramaters
    integer, parameter, dimension(16,16) :: alpha_parameters =	reshape((/1, 0, -3, 2, 0, 0, 0, 0, -3, 0, 9, -6, 2, 0, -6, &
        4, 0, 0, 3, -2, 0, 0, 0, 0, 0, 0, -9, 6, 0, 0, 6, -4, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, -9, 6, -2, 0, 6, -4, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 9, -6, 0, 0, -6, 4, 0, 1, -2, 1, 0, 0, 0, 0, 0, -3, 6, -3, 0, 2, -4, 2, 0, 0, -1, 1, 0, 0, 0, 0, &
        0, 0, 3, -3, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -6, 3, 0, -2, 4, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0,&
        0, 2, -2, 0, 0, 0, 0, 1, 0, -3, 2, -2, 0, 6, -4, 1, 0, -3, 2, 0, 0, 0, 0, 0, 0, 3, -2, 0, 0, -6, 4, 0, 0, 3, -2, 0,&
        0, 0, 0, 0, 0, 0, 0, -1, 0, 3, -2, 1, 0, -3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, 0, 3, -2, 0, 0, 0, 0, 0, 1,&
        -2, 1, 0, -2, 4, -2, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 2, -2, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,&
        2, -1, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1/), (/16,16/))

    contains
        subroutine load_elevation()
            call load_netcdf(elevation_netcdf_file, elev_x_coords, elev_y_coords, elev_grid)
            elev_grid_spacing = elev_x_coords(2) - elev_x_coords(1)
        end subroutine load_elevation

        subroutine load_shear_stress()
            call load_netcdf(ss_netcdf_file, ss_x_coords, ss_y_coords, ss_grid)
            ss_grid_spacing = ss_x_coords(2) - ss_x_coords(1)
        end subroutine load_shear_stress

        subroutine load_netcdf(netcdf_file, x_coords, y_coords, grid)
            use netcdf
            implicit none

            ! Name of the netcdf file to be read.
            character (len = *), intent(IN) :: netcdf_file
            integer :: ncid

            double precision, allocatable, dimension(:), intent(INOUT) :: x_coords, y_coords
            double precision, allocatable, dimension(:,:), intent(INOUT) :: grid

            ! Expected number and names of dims
            integer, parameter :: NDIMS = 2
            character (len = *), parameter :: Y_NAME = "y"
            character (len = *), parameter :: X_NAME = "x"
            character (len = *), parameter :: Z_NAME = "z"

            ! Actual number of dims, vars, attrs and unlim_dims
            integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in

            ! Dimension lengths and ids
            integer :: y_len, x_len
            integer :: y_dimid, x_dimid
            integer :: dimids(NDIMS)
            character (len = 128) :: x_dim_name, y_dim_name

            ! For the coordinate variables.
            integer :: y_varid, x_varid, z_varid

            ! Loop indices
            integer :: y_index, x_index, i, j

            ! Open netcdf
            call netcdf_check( nf90_open(netcdf_file, nf90_nowrite, ncid) )
            call netcdf_check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )

            ! Check we have the expected number of dims and that none are unlimited
            if (ndims_in /= 2 .or. unlimdimid_in /= -1) stop 2

            ! Get dim_ids for X and Y
            call netcdf_check( nf90_inq_dimid(ncid, Y_NAME, y_dimid) )
            call netcdf_check( nf90_inq_dimid(ncid, X_NAME, x_dimid) )

            ! Find length of dimensions
            call netcdf_check( nf90_inquire_dimension(ncid, y_dimid, y_dim_name, y_len) )
            call netcdf_check( nf90_inquire_dimension(ncid, x_dimid, x_dim_name, x_len) )

            ! Get the varids of the latitude and longitude coordinate variables.
            call netcdf_check( nf90_inq_varid(ncid, Y_NAME, y_varid) )
            call netcdf_check( nf90_inq_varid(ncid, X_NAME, x_varid) )

            ! Get the varids of the pressure and temperature netCDF variables.
            call netcdf_check( nf90_inq_varid(ncid, Z_NAME, z_varid) )
            allocate( y_coords(y_len) )
            allocate( x_coords(x_len) )
            allocate( grid(x_len, y_len) )

            ! Read the x and y data.
            call netcdf_check( nf90_get_var(ncid, y_varid, y_coords) )
            call netcdf_check( nf90_get_var(ncid, x_varid, x_coords) )

            ! read in topography data
            call netcdf_check( nf90_get_var(ncid, z_varid, grid) )

            ! Close the file.
            call netcdf_check( nf90_close(ncid) )
        end subroutine load_netcdf

        subroutine netcdf_check(status)
            use netcdf
            implicit none
            integer, intent ( in) :: status

            if(status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
        end subroutine netcdf_check

        subroutine grid_value(x_coord, y_coord, grid, grid_x, grid_y, grid_spacing, &
            local_grid_x_index, local_grid_y_index, alpha_array, value)
            implicit none

            double precision, intent(in) :: x_coord, y_coord
            double precision, dimension(:), intent(in) :: grid_x, grid_y
            double precision, dimension(:,:), intent(in) :: grid
            integer, intent(in) :: grid_spacing
            double precision, dimension(16), intent(inout) :: alpha_array
            double precision, intent(out) :: value

            integer :: leftmost_local_grid_x_index, bottom_local_grid_y_index
            integer, dimension(6), intent(inout) :: local_grid_x_index, local_grid_y_index
            double precision, dimension(6,6) :: local_grid
            double precision :: x_unit, y_unit
            integer :: counter, x_counter, y_counter

            leftmost_local_grid_x_index = (floor((x_coord-grid_x(1))/dble(grid_spacing)) + 1) - 2
            bottom_local_grid_y_index = (floor((y_coord-grid_y(1))/dble(grid_spacing)) + 1) - 2

            if( leftmost_local_grid_x_index /= local_grid_x_index(1) .or. &
                bottom_local_grid_y_index /= local_grid_y_index(1)) THEN

                local_grid_x_index(1) = leftmost_local_grid_x_index
                local_grid_y_index(1) = bottom_local_grid_y_index

                ! populate local grid index values
                do counter = 2, 6
                    local_grid_x_index(counter) = local_grid_x_index(counter-1) + 1
                    local_grid_y_index(counter) = local_grid_y_index(counter-1) + 1
                end do

                ! populate local grid values
                do x_counter = 1, 6
                    do y_counter = 1, 6
                        local_grid(x_counter, y_counter) = &
                            grid(local_grid_x_index(x_counter), local_grid_y_index(y_counter))
                    end do
                end do

                call bicubic_alpha(local_grid, grid_spacing, alpha_array)
            end if

            x_unit = (x_coord - grid_x(local_grid_x_index(3))) / dble(grid_spacing)
            y_unit = (y_coord - grid_y(local_grid_y_index(3))) / dble(grid_spacing)

            value = alpha_array(1) + alpha_array(2)*x_unit + alpha_array(3)*x_unit**2 + alpha_array(4)*x_unit**3&
                + alpha_array(5)*y_unit + alpha_array(6)*y_unit*x_unit +alpha_array(7)*x_unit**2*y_unit&
                + alpha_array(8)*x_unit**3*y_unit + alpha_array(9)*y_unit**2 &
                + alpha_array(10)*x_unit*y_unit**2 + alpha_array(11)*x_unit**2*y_unit**2 &
                + alpha_array(12)*x_unit**3*y_unit**2 + alpha_array(13)*y_unit**3 &
                + alpha_array(14)*x_unit*y_unit**3 + alpha_array(15)*x_unit**2*y_unit**3 &
                + alpha_array(16)*x_unit**3*y_unit**3
        end subroutine grid_value

        subroutine bicubic_alpha(in_array, bi_grid_spacing, bi_alpha_array)
            ! bicubic interpolation algorithm shamelessly lifted off of Wikipedia

            implicit none

            double precision, dimension(6,6), intent(in) :: in_array
            integer, intent(in) :: bi_grid_spacing
            double precision, dimension(16), intent(out) :: bi_alpha_array
            double precision, dimension(16) :: f_array
            double precision, dimension(6,6) :: x_der, y_der, xy_der
            integer :: x_counter, y_counter, counter

            do x_counter = 2, 5
                do y_counter = 2, 5

                    y_der(x_counter, y_counter) = ((in_array(x_counter+1,y_counter+1) + &
                        2.d0*in_array(x_counter,y_counter+1) + in_array(x_counter-1,y_counter+1)) - &
                        (in_array(x_counter+1,y_counter-1) + 2.d0*in_array(x_counter,y_counter-1) + &
                        in_array(x_counter-1,y_counter-1)))/(8.d0*bi_grid_spacing) 

                end do
            end do

            ! calculate the x and xy derivative

            do x_counter = 3, 4
                do y_counter = 3, 4

                    x_der(x_counter, y_counter) = ((in_array(x_counter+1,y_counter+1) + &
                        2.d0*in_array(x_counter+1,y_counter) + in_array(x_counter+1,y_counter-1)) - &
                        (in_array(x_counter-1,y_counter+1) + 2.d0*in_array(x_counter-1,y_counter) + &
                        in_array(x_counter-1,y_counter-1)))/(8.d0*bi_grid_spacing)

                    xy_der(x_counter, y_counter) = ((y_der(x_counter+1,y_counter+1) + &
                        2.d0*y_der(x_counter+1,y_counter) + y_der(x_counter+1,y_counter-1)) - &
                        (y_der(x_counter-1,y_counter+1) + 2.d0*y_der(x_counter-1,y_counter) + &
                        y_der(x_counter-1,y_counter-1)))/(8.d0*bi_grid_spacing)



                end do
            end do


            ! next determine the f values. The will be multiplied with the alpha_parameters to find the alpha values

            f_array(1) = in_array(3,3) ! f[0,0]
            f_array(2) = in_array(4,3) ! f[1,0]
            f_array(3) = in_array(3,4) ! f[0,1]
            f_array(4) = in_array(4,4) ! f[1,1]
            f_array(5) = x_der(3,3) ! fx[0,0]
            f_array(6) = x_der(4,3) ! fx[1,0]
            f_array(7) = x_der(3,4) ! fx[0,1]
            f_array(8) = x_der(4,4) ! fx[1,1]
            f_array(9) = y_der(3,3) ! fy[0,0]
            f_array(10) = y_der(4,3) ! fy[1,0]
            f_array(11) = y_der(3,4) ! fy[0,1]
            f_array(12) = y_der(4,4) ! fy[1,1]
            f_array(13) = xy_der(3,3) ! fxy[0,0]
            f_array(14) = xy_der(4,3) ! fxy[1,0]
            f_array(15) = xy_der(3,4) ! fxy[0,1]
            f_array(16) = xy_der(4,4) ! fxy[1,1]

            ! find alpha


            do counter = 1, 16

                bi_alpha_array(counter) = dot_product(alpha_parameters(counter,:), f_array)

            end do
        end subroutine bicubic_alpha

        double precision function elevation(x_coord, y_coord)
            implicit none
            double precision, intent(in) :: x_coord, y_coord
            double precision :: elevation_value

            call grid_value(x_coord, y_coord, elev_grid, elev_x_coords, elev_y_coords, elev_grid_spacing, &
                elev_local_grid_x_index, elev_local_grid_y_index, elev_alpha_array, elevation_value)

            elevation = elevation_value
        end function elevation

        double precision function shear_stress(x_coord, y_coord)
            implicit none
            double precision, intent(in) :: x_coord, y_coord
            double precision :: shear_stress_value

            call grid_value(x_coord, y_coord, ss_grid, ss_x_coords, ss_y_coords, ss_grid_spacing, &
                ss_local_grid_x_index, ss_local_grid_y_index, ss_alpha_array, shear_stress_value)

            shear_stress = shear_stress_value
        end function shear_stress

        double precision function y_gradient(x, y, direction)
            ! this function takes in a particular x and y location, and accounting for the flowline direction, calculates the gradient with respect to flowline y
            ! rotation is the direction that "x_flowline" points

            use global_parameters

            implicit none

            double precision, intent(in) :: x, y, direction
            double precision :: Q11, Q12, Q21, Q22, length, x_temp, y_temp

            length = sqrt(2.d0 * dx_l**2)

            x_temp = x + length*cos(5.d0*pi/4.d0+direction)
            y_temp = y + length*sin(5.d0*pi/4.d0+direction)
            Q11 = elevation(x_temp, y_temp)

            x_temp = x + length*cos(3.d0*pi/4.d0+direction)
            y_temp = y + length*sin(3.d0*pi/4.d0+direction)
            Q12 = elevation(x_temp, y_temp)

            x_temp = x + length*cos(-1.d0*pi/4.d0+direction)
            y_temp = y + length*sin(-1.d0*pi/4.d0+direction)
            Q21 = elevation(x_temp, y_temp)

            x_temp = x + length*cos(1.d0*pi/4.d0+direction)
            y_temp = y + length*sin(1.d0*pi/4.d0+direction)
            Q22 = elevation(x_temp, y_temp)

            y_gradient = (-Q11  - Q21 + Q12 + Q22) / (4.d0 * dx_l)
        end function y_gradient

        double precision function Hf_gradient(x, y, direction)	
            ! this function takes in a particular x and y location, and accounting for the flowline direction, calculates the gradient with respect to flowline y
            ! rotation is the direction that "x_flowline" points

            use global_parameters

            implicit none

            double precision, intent(in) :: x, y, direction
            double precision :: Q11, Q12, Q21, Q22, length, x_temp, y_temp

            length = sqrt(2.d0 * dx_l**2)

            x_temp = x + length*cos(5.d0*pi/4.d0+direction)
            y_temp = y + length*sin(5.d0*pi/4.d0+direction)
            Q11 = shear_stress(x_temp, y_temp)

            x_temp = x + length*cos(3.d0*pi/4.d0+direction)
            y_temp = y + length*sin(3.d0*pi/4.d0+direction)
            Q12 = shear_stress(x_temp, y_temp)

            x_temp = x + length*cos(-1.d0*pi/4.d0+direction)
            y_temp = y + length*sin(-1.d0*pi/4.d0+direction)
            Q21 = shear_stress(x_temp, y_temp)

            x_temp = x + length*cos(1.d0*pi/4.d0+direction)
            y_temp = y + length*sin(1.d0*pi/4.d0+direction)
            Q22 = shear_stress(x_temp, y_temp)

            Hf_gradient = (-Q11  - Q21 + Q12 + Q22) / ((4.d0 * dx_l) * rho_ice * g) 
        end function Hf_gradient
end module grids
