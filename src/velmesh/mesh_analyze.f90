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
!    use omp_lib
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
    character(128)                      :: INFILE_4, INFILE_5
    integer                             :: i, j, k, l, numSteps  
        ! numSteps is the # of desired time steps to take
    !integer                             :: N            
        ! Number of neutron collisions
    real(8)                             :: Dm, vel, theta0, Vopt, Wsc, Wabs
    real(8)                             :: thick, A, B, tStep, x, lambda!, ang
    !complex(8)                          :: i
    real(8),        dimension(2, 2)     :: rho ! rho is the dens. matrix
    complex(8),     dimension(2)        :: psi  ! initial n-n' state, e.g.,(1,0)
    !character(80),  dimension(3)        :: arg_in
    
    character(prec), dimension(2, 2)    :: rhonn
    character(prec)                     :: uni, dist, thetaChar, massChar, xChar
    character(prec)                     :: DmChar, theta0char, thickness
    character(36)                       :: uuid
    character(64)                       :: filename

    integer                             :: num_lines, num_vels, nMasses, nAngles
    type(materiallist), dimension(:), allocatable   :: inventory
    logical                             :: only_endpoint, no_scattering
    logical                             :: no_absorption

    real(8), dimension(:, :), allocatable   :: rho_results
    real(8), dimension(4)                   :: rho_ave, rho_var
    real(8), dimension(:), allocatable  :: vel_list, Masses, Angles

    type(mesh_file_data)                :: mfp
    character(256)                      :: directoryName
    character(256)                      :: file1, file2, file3, file4
    character(256)                      :: file5, file6, file7, file8
    character(256)                      :: errorlog

    real(8), dimension(:, :), allocatable    :: mesh1
    real(8), dimension(:, :), allocatable    :: mesh2
    real(8), dimension(:, :), allocatable    :: mesh3
    real(8), dimension(:, :), allocatable    :: mesh4
    real(8), dimension(:, :), allocatable    :: mesh5
    real(8), dimension(:, :), allocatable    :: mesh6
    real(8), dimension(:, :), allocatable    :: mesh7
    real(8), dimension(:, :), allocatable    :: mesh8
! === End variable declarations ================================================

! === Begin variable assignments ===============================================
    numSteps    = int(1e3)
    x           = 0.D0
    lambda      = 0.D0

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
    DmChar          = trim(DmChar)
    theta0char      = trim(theta0char)

    !rho(1, 1)   = 1.0
    !rho(2, 1)   = 0.0
    !rho(1, 2)   = 0.0
    !rho(2, 2)   = 0.0

    write(ff58, "(I2)") prec
    f58 = "(F" // trim(adjustl(ff58)) // ".11)"

    ! For the moment, hard-code location of "parameter.list" and "material.list"
    ! That can be addressed later
    INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"
    INFILE_3 = "velocity.list"
    INFILE_4 = "deltaMs"
    INFILE_5 = "thetas"

    errorlog = "error-" // uuid // ".log"

    ! Check the command line for parameters
    call get_vel_params(INFILE_2, INFILE_4, INFILE_5, inventory, psi, &
                    &filename, num_lines, only_endpoint, no_scattering, &
                    &no_absorption, Masses, Angles, nMasses, nAngles)

    rho(1, 1) = conjg(psi(1)) * psi(1)
    rho(2, 2) = conjg(psi(2)) * psi(2)
    rho(1, 2) = psi(1) * conjg(psi(2))
    rho(2, 1) = psi(2) * conjg(psi(1))

    write( DmChar,     trim(adjustl(f58)) ) Dm
    write( theta0char, trim(adjustl(f58)) ) theta0

    ! Shuffle this into get_params for easier running
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
        &x, rho, no_scattering, no_absorption, prec)
    ! Write each file to a directory of the appropriate mass
    
    call mesh_file_prepare_ave(mfp, directoryname, file1, rhonn(1, 1), 1)
    call mesh_file_prepare_var(mfp, directoryname, file5, rhonn(1, 1), 5)
    call mesh_file_prepare_ave(mfp, directoryname, file2, rhonn(1, 2), 2)
    call mesh_file_prepare_var(mfp, directoryname, file6, rhonn(1, 2), 6)
    call mesh_file_prepare_ave(mfp, directoryname, file3, rhonn(2, 1), 3)
    call mesh_file_prepare_var(mfp, directoryname, file7, rhonn(2, 1), 7)
    call mesh_file_prepare_ave(mfp, directoryname, file4, rhonn(2, 2), 4)
    call mesh_file_prepare_var(mfp, directoryname, file8, rhonn(2, 2), 8)

! === End initial file I/O =====================================================

! === Begin main loop ==========================================================
    ! Since there may be multiple materials, and there may be multiple runs
    ! over multiple parameters for those multiple materials, each material needs
    ! a full treatment - maybe not a million points each, but a full treatment
    ! nonetheless. 
    ! Loop over the number of materials in `inventory`, then loop from 1 to
    ! numSteps in each material.

!    ! write to averaged file here
!    call mesh_file_write(rho_ave, rho_var)

    ! Shared: Mesh array (critical), vel_list
    ! Private: values of rho_results
    allocate(mesh1(nMasses, nAngles))
    allocate(mesh2(nMasses, nAngles))
    allocate(mesh3(nMasses, nAngles))
    allocate(mesh4(nMasses, nAngles))
    allocate(mesh5(nMasses, nAngles))
    allocate(mesh6(nMasses, nAngles))
    allocate(mesh7(nMasses, nAngles))
    allocate(mesh8(nMasses, nAngles))

    !$OMP PARALLEL PRIVATE(Dm, theta0, Vopt, Wsc, Wabs, thick, rho_results, rho_ave, rho_var, i, j, k, l, vel)
    !$OMP DO
    num_vels = 1

    do i = 1, nMasses
        do j = 1, nAngles
            do k = 1, 1!num_vels
                vel = vel_list(k)
                
                Dm = Masses(i)
                theta0 = Angles(j)

                call iterate_probability(inventory, vel, errorlog, numSteps, &
                    &no_scattering, no_absorption, only_endpoint, Dm, theta0, &
                    &lambda, A, B, psi, rho, num_lines)

                rho_results(k, 1) = rho(1, 1)
                rho_results(k, 2) = rho(1, 2)
                rho_results(k, 3) = rho(2, 1)
                rho_results(k, 4) = rho(2, 2)

                exit
            end do
        
            do k = 1, 4
                rho_ave(k) = sum(rho_results(:,k)) / size(vel_list)
                rho_var(k) = compute_variance(rho_results(:,k))
            end do
            !$OMP CRITICAL
            mesh1(i, j) = rho_ave(1)
            mesh2(i, j) = rho_ave(2)
            mesh3(i, j) = rho_ave(3)
            mesh4(i, j) = rho_ave(4)
            mesh5(i, j) = rho_var(1)
            mesh6(i, j) = rho_var(2)
            mesh7(i, j) = rho_var(3)
            mesh8(i, j) = rho_var(4)
            !$OMP END CRITICAL
!            goto 001
        end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL
! === End main loop ============================================================
    ! Write out the goods, and also the corresponding masses/angles for each
001 continue
51  format(A1, A15)

    open(unit = 1, file = trim(adjustl(directoryName)) &
                        &// "masses.txt", status = "unknown")
    open(unit = 2, file = trim(adjustl(directoryName)) &
                        &// "angles.txt", status = "unknown")
    write(unit = 1, fmt = 51) "#", "Mass Delta"
    write(unit = 2, fmt = 51) "#", "Angle"

    do i = 1, nMasses
        write(unit = 1, fmt = 50) Masses(i)
    end do
    close(unit = 1)

    do i = 1, nAngles
        write(unit = 2, fmt = 50) Angles(i)
    end do
    close(unit = 2)

    do i = 1, nMasses
        write(unit = 10, fmt = 52) mesh1(i, :)
        write(unit = 20, fmt = 52) mesh2(i, :)
        write(unit = 30, fmt = 52) mesh3(i, :)
        write(unit = 40, fmt = 52) mesh4(i, :)
        write(unit = 50, fmt = 52) mesh5(i, :)
        write(unit = 60, fmt = 52) mesh6(i, :)
        write(unit = 70, fmt = 52) mesh7(i, :)
        write(unit = 80, fmt = 52) mesh8(i, :)
    end do

    close(unit = 10)
    close(unit = 20)
    close(unit = 30)
    close(unit = 40)
    close(unit = 50)
    close(unit = 60)
    close(unit = 70)
    close(unit = 80)
50  format(ES17.8E3)
52  format(1608ES17.8E3)
contains
    real(8) function compute_variance(array) result(var)
    ! Use the Two-Pass formulation because dead simple
        implicit none
        real(8), dimension(:), intent(in)   :: array
        real(8)                             :: average
        integer                             :: m

        average = 0.D0
        var     = 0.D0

        average = sum(array) / size(array)

        do m = 1, size(array)
            var = var + (array(m) - average)**2 / (size(array) - 1)
        end do
        !var = var / (size(array) - 1)
    end function compute_variance

    subroutine iterate_probability(inventory, vel, errorlog, numSteps, &
        &no_scattering, no_absorption, only_endpoint, Dm, theta0, lambda, A, &
        &B, psi, rho, num_lines)
        implicit none
        type(materiallist), dimension(:), intent(in)    :: inventory
        real(8)         :: vel, Dm, theta0, lambda, A, B
        integer         :: numSteps, num_lines, n!, u
        character(256)  :: errorlog
        logical         :: no_scattering, no_absorption, only_endpoint
        real(8), dimension(2, 2)    :: rho
        complex(8), dimension(2)    :: psi

        do l = 1, num_lines
            Vopt    = inventory(l)%V    * 1.D-9 ! Given in neV, want eV
            Wsc     = inventory(l)%Wsc  * 1.D-9 ! ^; W_sc
            Wabs    = inventory(l)%Wabs * 1.D-9 ! ^; W_abs
            thick   = inventory(l)%d    / 100   ! Was given in cm, want in m

            if (vel .le. 1.D0) then
                open(unit = 99, file = errorlog)
                vel = 1.
                write(unit = 99, fmt = *) "One value of vel set to 1 due to being too small."
                close(unit = 99)
            end if

            tStep   = thick / ( vel * numSteps )
            
            !ang = 2*(0.5 * atan(0.5 * abs(Dm) * tan(2. * theta0) / (Vopt - Dm)))**2
    ! Generate (2theta^2, 1-2theta^2)
            !psi(1)      = ang
            !psi(2)      = 1 - ang
    ! Generate (1-2theta^2, theta^2)
            !psi(1)      = 1 - ang
            !psi(2)      = ang
            if (only_endpoint) then
                call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, lambda,A,B,&
                                &thick / vel, psi, rho)
            else
                do 200 n = 1, numSteps
                    call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, &
                                    &A, B, n*tStep, psi, rho)
                    print *, "TODO: Specify a unit file to write this to, or remove this whole bit."
                    stop
                    write(1, 57) x, n*tStep*vel, rho(1, 1), rho(2, 1), rho(1, 2), &
                        &rho(2, 2)
                    x = x + tStep * vel
200             end do
            end if


    ! Update psi with the diagonal elements of rho
    ! Used for if multiple materials that interface
!            if (num_lines .gt. 1) then
!                psi(1) = sqrt(rho(1, 1))
!                psi(2) = sqrt(rho(2, 2))
!            end if
        end do

57  format(6ES16.6E3)
    end subroutine iterate_probability

end program stereo_sim
