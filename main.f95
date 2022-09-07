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

! === Begin variable declarations ==============================================
    implicit none
    integer, parameter                  :: prec = 16
    character(128)                      :: INFILE_1, INFILE_2
    integer                             :: j, k, numSteps  
        ! numSteps is the # of desired time steps to take
    !integer                             :: N            
        ! Number of neutron collisions
    real(8)                             :: Dm, vel, theta0, Vopt, Wsc, Wabs,thick
    real(8)                             :: A, B, tStep, ang, x
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

    integer                             :: num_lines
    type(materiallist), dimension(:), allocatable   :: inventory
        

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


    58 format(F16.11)

    ! For the moment, hard-code location of "parameter.list" and "material.list"
    ! That can be addressed later
    INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"

    ! Check the command line for parameters
    num_lines = get_lines(INFILE_2)
    allocate(inventory(num_lines))
    call get_params(INFILE_1, Dm, theta0, psi)

    !print *, "Dm        :", Dm
    !print *, "theta0    :", theta0

    write( DmChar,      58 ) Dm
    write( theta0char,  58 ) theta0

    ! Shuffle this into get_params for easier running
    call get_materials(inventory, filename, vel, num_lines, INFILE_2)
    write( thickness,   58 ) sum(inventory%d)
! === End variable assignments =================================================
    
! === Begin initial file I/O ===================================================
    ! Write each file to a directory of the appropriate mass
    if (psi(1) == cmplx(1.0, 0.0, 8) .and. psi(2) == cmplx(0.0, 0.0, 8)) then
        directoryName = "./data/N/" // trim(filename) // "-" // trim(adjustl(thickness)) // &
        &"/" // trim(adjustl(DmChar)) // "/"
    else if (psi(1) == cmplx(0.0, 0.0, 8) .and. psi(2) == cmplx(1.0, 0.0, 8)) then
        directoryName = "./data/Np/" // trim(filename) // "-" // trim(adjustl(thickness)) // &
        &"/" // trim(adjustl(DmChar)) // "/"
    else
        directoryName = "./data/" // trim(filename) // "-" // trim(adjustl(thickness)) // &
        &"/" // trim(adjustl(DmChar)) // "/"
    end if
    call system('mkdir -p ' // trim (adjustl(directoryName)))

    open(unit = 1, file = trim(adjustl(directoryName)) // trim(filename) // "_Dm-" &
            &// trim(adjustl(DmChar)) // "_theta0-" // &
            &trim(adjustl(theta0char)) // "_data_" // uuid // ".dat", &
            &status = "unknown")
    write(1, *) '#' // massChar, thetaChar
    write(1, '(2F16.11)') Dm, theta0
    write(1, 56) '#' // xChar, dist, rhonn(1, 1), rhonn(2, 1), rhonn(1, 2), rhonn(2, 2)
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
        Vopt    = inventory(j)%V    * 1.D-9 ! Given in neV, want eV
        Wsc     = inventory(j)%Wsc  * 1.D-9 ! ^; W_sc
        Wabs    = inventory(j)%Wabs * 1.D-9 ! ^; W_abs
        thick   = inventory(j)%d    / 100   ! Was given in cm, want in m

        tStep   = thick / ( vel * numSteps )

        print *, "Vopt      :", Vopt
        print *, "Wsc       :", Wsc
        print *, "Wabs      :", Wabs
        print *, "thick     :", thick
        print *, "numSteps  :", numSteps
        
        ang = 2*(0.5 * atan(0.5 * abs(Dm) * tan(2. * theta0) / (Vopt - Dm)))**2
! Generate (2theta^2, 1-2theta^2)
        !psi(1)      = ang
        !psi(2)      = 1 - ang
! Generate (1-2theta^2, theta^2)
        !psi(1)      = 1 - ang
        !psi(2)      = ang
! Generate (1, 0)
        !psi(1)      = cmplx(1.0, 0.0, 8)
        !psi(2)      = cmplx(0.0, 0.0, 8)

        do 200 k = 1, numSteps
            call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, thick, A, B, &
                &k*tStep, psi, rho)
            write(1, 57) x, k*tStep*vel, rho(1, 1), rho(2, 1), rho(1, 2), &
                &rho(2, 2)
!            write(6, 57) x, k*tStep*vel, rho(1, 1), rho(2, 1), rho(1, 2), &
!                &rho(2, 2)
            x = x + tStep * vel
200     end do

    ! Update psi with the diagonal elements of rho
    psi(1)  = sqrt(rho(1, 1))
    psi(2)  = sqrt(rho(2, 2))
100 end do
! === End main loop ============================================================

56 format(A19, 5A16)
57 format(6ES16.6E3)
    close(unit = 1)
end program stereo_sim