! September 21st, 2022
! Cary Rock
! This program program takes in a list of velocities (sampled randomly from an
! appropriate distribution) and averages the results for a given set of input
! parameters, dM and theta. This, too, is to aid in the analysis of the STEREO
! '22 results (indirectly) at the ILL.

! U = V +/- Dm +/- uB
! W = Wabs + v * Wsc
! Dm = m_n' - m_n in units of eV
! eps = 

! https://github.com/haniibrahim/f90getopt - ! GH repo for Fortran (>='03) 
! C-style commandline arguments
! 

program stereo_sim
    use exact_banfor_module 
    use uuid_module         ! Generates a UUID for each file
    use material_list       ! Methods for reading in an inventory list for
                            ! simulating on ("material.list")
    use get_parameters      ! Checks the command line for parameters, or reads
                            ! them from "parameter.list"
    use NNp_file_writing

! === Begin variable declarations ==============================================
    implicit none
    integer, parameter                  :: prec = 18
    character(2)                        :: ff58
    character(20)                       :: f58
    character(128)                      :: INFILE_1, INFILE_2, INFILE_3
    integer                             :: j, k, numSteps  
        ! numSteps is the # of desired time steps to take
    !integer                             :: N            
        ! Number of neutron collisions
    real(8)                             :: Dm, vel, theta0, Vopt, Wsc, Wabs
    real(8)                             :: thick, A, B, tStep, x!, ang
    !complex(8)                          :: i
    real(8),        dimension(2, 2)     :: rho  ! rho is the density matrix
    complex(8),     dimension(2)        :: psi  ! initial n-n' state, e.g.,(1,0)
    !character(80),  dimension(3)        :: arg_in
    
    character(prec), dimension(2, 2)    :: rhonn
    character(prec)                     :: uni, dist, thetaChar, massChar, xChar
    character(prec)                     :: DmChar, theta0char, thickness
    character(36)                       :: uuid
    character(64)                       :: filename

    integer                             :: num_lines, num_vels
    type(materiallist), dimension(:), allocatable   :: inventory
    logical                             :: only_endpoint, no_scattering

    real(8), dimension(:, :), allocatable   :: rho_results
    real(8), dimension(4)                   :: rho_ave, rho_var
    real(8), dimension(:), allocatable  :: vel_list

    type(mesh_file_data)                :: mfp
    character(256)                      :: directoryName, file
! === End variable declarations ================================================

! === Begin variable assignments ===============================================
    numSteps    = int(1e3)
    x           = 0.D0

    ! This is just for nice formatting
    xChar           = "Global X"
    rhonn(1, 1)     = "pnn"
    rhonn(2, 1)     = "pmn"
    rhonn(1, 2)     = "pnm"
    rhonn(2, 2)     = "pmm"
    thetaChar       = "theta"
    uni             = "unitar"
    massChar        = "mass"
    dist            = "distance"
    uuid            = generate_UUID()
    DmChar          = trim((DmChar))
    theta0char      = trim((theta0char))

    rho(1, 1)   = 1.0
    rho(2, 1)   = 0.0
    rho(1, 2)   = 0.0
    rho(2, 2)   = 0.0

    write(ff58, "(I2)") prec
    f58 = "(F" // trim(adjustl(ff58)) // ".11)"

    ! For the moment, hard-code location of "parameter.list" and "material.list"
    ! That can be addressed later
    INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"
    INFILE_3 = "velocity.list"

    ! Check the command line for parameters
    ! TODO: get_params needs num_lines to allocate stuff, but num_lines needs 
    ! the actual INFILE_2 file to get the number of lines. Circular.
    !num_lines = get_lines(INFILE_2)
    call get_params(INFILE_1, INFILE_2, inventory, Dm, &
                    &theta0, psi, vel, filename, num_lines, &
                    &only_endpoint, no_scattering)

    write( DmChar,     trim(adjustl(f58)) ) Dm
    write( theta0char, trim(adjustl(f58)) ) theta0

    ! Shuffle this into get_params for easier running
    !call get_materials(inventory, filename, vel, num_lines, INFILE_2)
    write( thickness,  f58 ) sum(inventory%d)

    ! Get the number of velocities to average over - this is definitely going
    ! to need to be refactored
    num_vels = get_vel_lines(INFILE_3)
    allocate(rho_results(num_vels, 4))
    allocate(vel_list(num_vels))
    call get_velocities(INFILE_3, num_vels, vel_list)

! === End variable assignments =================================================
    
! === Begin initial file I/O ===================================================
    mfp = create_msd_struct(psi, filename, thickness, DmChar, theta0char, &
        &uuid, massChar, thetaChar, Dm, theta0, xChar, dist, rhonn, &
        &x, rho, no_scattering, prec)
    ! Write each file to a directory of the appropriate mass
    call mesh_file_prepare(mfp, directoryname, file)

! === End initial file I/O =====================================================

! === Begin main loop ==========================================================
    ! Since there may be multiple materials, and there may be multiple runs
    ! over multiple parameters for those multiple materials, each material needs
    ! a full treatment - maybe not a million points each, but a full treatment
    ! nonetheless. 
    ! Loop over the number of materials in `inventory`, then loop from 1 to
    ! numSteps in each material.

    do 400 k = 1, num_vels
        vel = vel_list(k)

        call iterate_probability()
        rho_results(k, 1) = rho(1, 1)
        rho_results(k, 2) = rho(1, 2)
        rho_results(k, 3) = rho(2, 1)
        rho_results(k, 4) = rho(2, 2)
        
400 end do
        
    do k = 1, 4
        rho_ave(k) = sum(rho_results(:,k)) / size(rho_results(:,k))
        rho_var(k) = compute_variance(rho_results(:,k))
    end do

    ! write to averaged file here
    call mesh_file_write(rho_ave, rho_var)
! === End main loop ============================================================

    close(unit = 1)

contains
    real(8) function compute_variance(array) result(var)
    ! Use the Two-Pass formulation because dead simple
        implicit none
        real(8), dimension(:), intent(in)   :: array
        real(8)                             :: average
        integer                             :: i

        average = 0.D0
        var     = 0.D0

        average = sum(array) / size(array)

        do i = 1, size(array)
            var = var + (array(i) - average)**2 / (size(array) - 1)
        end do
        !var = var / (size(array) - 1)
    end function compute_variance

    subroutine iterate_probability()
        implicit none
        integer     :: n
        do 100 j = 1, num_lines
            Vopt    = inventory(j)%V    * 1.D-9 ! Given in neV, want eV
            Wsc     = inventory(j)%Wsc  * 1.D-9 ! ^; W_sc
            Wabs    = inventory(j)%Wabs * 1.D-9 ! ^; W_abs
            thick   = inventory(j)%d    / 100   ! Was given in cm, want in m

            tStep   = thick / ( vel * numSteps )

        
            !ang = 2*(0.5 * atan(0.5 * abs(Dm) * tan(2. * theta0) / (Vopt - Dm)))**2
    ! Generate (2theta^2, 1-2theta^2)
            !psi(1)      = ang
            !psi(2)      = 1 - ang
    ! Generate (1-2theta^2, theta^2)
            !psi(1)      = 1 - ang
            !psi(2)      = ang

            if (only_endpoint) then
                    call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, thick, A, B,&
                    &thick,   psi, rho)
            else
                do 200 n = 1, numSteps
                    call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, thick, A, B,&
                        &n*tStep, psi, rho)
                    write(1, 57) x, n*tStep*vel, rho(1, 1), rho(2, 1), rho(1, 2), &
                        &rho(2, 2)
                    x = x + tStep * vel
200             end do
        end if

    ! Update psi with the diagonal elements of rho
        psi(1) = sqrt(rho(1, 1))
        psi(2) = sqrt(rho(2, 2))
100     end do

57  format(6ES16.6E3)
    end subroutine iterate_probability

end program stereo_sim
