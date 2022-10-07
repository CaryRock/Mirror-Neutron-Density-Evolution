! August 21st, 2022
! Cary Rock
! The "driver" program for the ILL reactor simulation, specifically towards
! in/validating the results of the Sarrazin STEREO '22 paper's findings 
! regarding mirror neutrons

! U = V +/- Dm +/- uB
! W = Wabs + v * Wsc
! Dm = m_n' - m_n in units of eV
! eps = 

! https://github.com/haniibrahim/f90getopt - ! GH repo for Fortran (>='03) 
! C-style commandline arguments
! 

program stereo_sim
    !use banfor_module       ! Old way - doesn't account for initial phase
    use exact_banfor_module ! Current method - does account for
    use uuid_module         ! Generates a UUID for each file
    use material_list       ! Methods for reading in an inventory list for
                            ! simulating on ("material.list")
    use get_parameters      ! Checks the command line for parameters, or
                            ! reads them from "parameter.list"
    use NNp_file_writing

! === Begin variable declarations ==============================================
    implicit none
    integer, parameter                  :: prec = 16
    character(128)                      :: INFILE_1, INFILE_2, INFILE_3
    integer                             :: j, k, numSteps  
        ! numSteps is the # of desired time steps to take
    !integer                             :: N            
        ! Number of neutron collisions
    real(8)                             :: Dm, vel, theta0, Vopt, Wsc, Wabs
    real(8)                             :: thick, A, B, tStep, ang, x, lambda
    !complex(8)                          :: i
    real(8),        dimension(2, 2)     :: rho  ! rho is the density matrix
    complex(8),     dimension(2)        :: psi  ! initial n-n' state, e.g.,(1,0)
    !character(80),  dimension(3)        :: arg_in
    
    character(prec), dimension(2, 2)    :: rhonn
    character(prec)                     :: uni, dist, thetaChar, massChar, xChar
    character(80)                       :: directoryName
    character(prec)                     :: DmChar, theta0char, thickness
    character(36)                       :: uuid
    character(64)                       :: filename

    integer                             :: num_lines, num_vels
    type(materiallist), dimension(:), allocatable   :: inventory
    logical                             :: only_endpoint, no_scattering

    real(8), dimension(:, :), allocatable   :: rho_results
    real(8), dimension(4)                   :: rho_av, rho_var
    real(8), dimension(:), allocatable      :: vel_list!, Masses, Angles

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
    DmChar          = trim((DmChar))
    theta0char      = trim((theta0char))

    rho(1, 1)   = 1.0
    rho(2, 1)   = 0.0
    rho(1, 2)   = 0.0
    rho(2, 2)   = 0.0


    58 format(F16.11)

    ! For the moment, hard-code location of "parameter.list" and "material.list"
    ! That can be addressed later
    INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"
    INFILE_3 = "velocity.list"

    ! Check the command line for parameters
    !num_lines = get_lines(INFILE_2)
    !allocate(inventory(num_lines))
    call get_params(INFILE_1, INFILE_2, inventory, Dm, &
                    &theta0, psi, vel, filename, num_lines, &
                    &only_endpoint, no_scattering)

    write( DmChar,      58 ) Dm
    write( theta0char,  58 ) theta0

    ! Shuffle this into get_params for easier running
    !call get_materials(inventory, filename, vel, num_lines, INFILE_2)
    write( thickness,   58 ) sum(inventory%d)

    ! Get the number of velocities to average over - this is definitely going
    ! to need to be refactored
    num_vels = get_vel_lines(INFILE_3) - 1
    !print *, "Number of velocities: ", num_vels
    allocate(rho_results(num_vels, 4))
    allocate(vel_list(num_vels))
    call get_velocities(INFILE_3, num_vels, vel_list)

! === End variable assignments =================================================
    
! === Begin initial file I/O ===================================================
    ! Write each file to a directory of the appropriate mass
    ! Use the file_writing module
    if (psi(1) == cmplx(1.0, 0.0, 8) .and. psi(2) == cmplx(0.0, 0.0, 8)) then
        directoryName = "./data/N/" // trim(filename) // "-" // &
            &trim(adjustl(thickness)) // "/" // trim(adjustl(DmChar)) // "/"
    else if (psi(1) == cmplx(0.0, 0.0, 8) .and. &
                                            &psi(2) == cmplx(1.0, 0.0, 8)) then
        directoryName = "./data/Np/" // trim(filename) // "-" // &
            &trim(adjustl(thickness)) // "/" // trim(adjustl(DmChar)) // "/"
    else
        directoryName = "./data/" // trim(filename) // "-" // &
            &trim(adjustl(thickness)) // "/" // trim(adjustl(DmChar)) // "/"
    end if
    call system('mkdir -p ' // trim (adjustl(directoryName)))

    open(unit = 1, file = trim(adjustl(directoryName)) // trim(filename) // &
            &"_Dm-" // trim(adjustl(DmChar)) // "_theta0-" // &
            &trim(adjustl(theta0char)) // "_data_" // uuid // ".dat", &
            &status = "unknown")

! TODO: Finish converting this section to use the file_writing module
    !mfp = create_msd_struct(psi, filename, thickness, DmChar, theta0char, uuid,&
    !                        &massChar, Dm, theta0, xChar, dist, rhonn, x, rho, &
    !                        &no_scattering, prec)
    !call coord_file_prepare(mfp, directoryName, file, element, u)

    write(1, *) '#' // massChar, thetaChar
    write(1, '(2F16.11)') Dm, theta0
    write(1, 56) '#' // xChar, dist, rhonn(1, 1), rhonn(2, 1), &
                                                    &rhonn(1, 2), rhonn(2, 2)
    write(1, 57) x, 0.0, rho(1, 1), rho(2, 1), rho(1, 2), rho(2, 2)

! === End initial file I/O =====================================================

! === Begin main loop ==========================================================
    ! Since there may be multiple materials, and there may be multiple runs
    ! over multiple parameters for those multiple materials, each material needs
    ! a full treatment - maybe not a million points each, but a full treatment
    ! nonetheless. 
    ! Loop over the number of materials in `inventory`, then loop from 1 to
    ! numSteps in each material.

    do 100 j = 1, num_lines
        vel = vel_list(1)   ! Kludge - fix later
        Vopt    = inventory(j)%V    * 1.D-9 ! Given in neV, want eV
        Wsc     = inventory(j)%Wsc  * 1.D-9 ! ^; W_sc
        Wabs    = inventory(j)%Wabs * 1.D-9 ! ^; W_abs
        thick   = inventory(j)%d    / 100   ! Was given in cm, want in m

        tStep   = thick / ( vel * numSteps )

        print *, "tStep, vel, numSteps, Vopt, Wsc, Wabs, thick: ", tStep, vel, &
                &numSteps, Vopt, Wsc, Wabs, thick
        
        !ang = 2*(0.5 * atan(0.5 * abs(Dm) * tan(2. * theta0) / (Vopt - Dm)))**2
! Generate (2theta^2, 1-2theta^2)
        !psi(1)      = ang
        !psi(2)      = 1 - ang
! Generate (1-2theta^2, theta^2)
        !psi(1)      = 1 - ang
        !psi(2)      = ang
! Generate (1, 0)
        !psi(1)      = cmplx(1.0, 0.0, 8)
        !psi(2)      = cmplx(0.0, 0.0, 8)

        if (only_endpoint) then
            do 400 k = 1, num_vels
                call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, A, B,&
                                &thick/vel, psi, rho)
            
                rho_results(k, 1) = rho(1, 1)
                rho_results(k, 2) = rho(1, 2)
                rho_results(k, 3) = rho(2, 1)
                rho_results(k, 4) = rho(2, 2)
400         end do
            
            do k = 1, 4
                rho_av(k) = sum(rho_results(:,k)) / num_vels
                rho_var(k)= compute_variance(rho_results(:,k))
            end do

            ! write to averaged file here
            !write(1, 59) rho_av(1), rho_av(2), rho_av(3), rho_av(4), &
            !    &rho_var(1), rho_var(2), rho_var(3), rho_var(4)
            write(1, 57) x, numSteps * tStep, rho(1, 1), rho(2, 1), rho(1, 2),&
                        &rho(2, 2)
        else
            do 200 k = 1, numSteps
                call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, lambda, A, B,&
                                &k*tStep, psi, rho)
                write(1, 57) x, k*tStep*vel, rho(1, 1), rho(2, 1), rho(1, 2), &
                            &rho(2, 2)
                x = x + tStep * vel
200         end do
        end if

    ! Update psi with the diagonal elements of rho
    psi(1)  = sqrt(rho(1, 1))
    psi(2)  = sqrt(rho(2, 2))
100 end do
! === End main loop ============================================================

56 format(A19, 5A16)
57 format(6ES16.6E3)
59 format(8ES16.6E3)
    close(unit = 1)

contains
    real(8) function compute_variance(array) result (var)
        implicit none
        real(8), dimension(:), intent(in)   :: array
        real(8)                             :: average
        integer                             :: i

        average = sum(array)
        var = 0.D0

        do i = 1, size(array)
            var = var + (array(i) - average)**2
        end do
        var = var / size(array)
    end function compute_variance
end program stereo_sim
