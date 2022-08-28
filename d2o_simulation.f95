! This is a section of code that is to simulation the flight of a neutron
! D2O, colliding N times. It is a more-simple case, B = 0 T, but that
! extension would be simple to implement.

! U = V +/- Dm +/- uB
! W = W1 + v * W2
! Dm = m_n' - m_n in units of eV
! eps = 

! https://github.com/haniibrahim/f90getopt - ! GH repo for Fortran (>='03) 
! C-style commandline arguments
! Can either read in commandline arguments somehow (like from this), or devise
! a way to launch parallel (for the full simulation) from a file (somehow)

program d2o
    use banfor_module       ! Old way - doesn't account for initial phase
    use exact_banfor_module ! Current method - does account for
    use uuid_module         ! Generates a UUID for each file
    use material_list       ! Methods for reading in an inventory list for
                            ! simulating on ("material.list")
    use get_parameters      ! Checks the command line for parameters, or
                            ! reads them from "parameter.list"

! === Begin variable declarations ==============================================
    implicit none
    integer, parameter                  :: prec = 16
    integer                             :: j, numSteps  
        ! numSteps is the # of desired time steps to take
    integer                             :: N            
        ! Number of neutron collisions
    real(8)                             :: Dm, vel, theta0, Vopt, W1, W2, thick
    real(8)                             :: A, B, tStep
    complex(8) i
    real(8),        dimension(2, 2)     :: rho  ! rho is the density matrix
    complex(8),     dimension(2)        :: psi  ! initial n-n' state, e.g.,(1,0)
    character(80),  dimension(3)        :: arg_in
    
    character(prec), dimension(2, 2)    :: rhonn
    character(prec)                     :: uni, dist, thetachar, massChar
    character(80)                       :: directoryName
    character(prec)                     :: DmChar, theta0char, thickness
    character(36)                       :: uuid
    character(64)                       :: filename

    integer                             :: num_lines
    type(materiallist), dimension(:), allocatable   :: inventory
        

! === End variable declarations ================================================

! === Begin variable assignments ===============================================

58  format(F16.11)

    ! Check the command line for parameters
    call get_params(Dm, theta0)

    ! Check the material.list file for inventory
    num_lines = get_lines("material.list")
    call get_materials(inventory, filename, vel, num_lines, "material.list")

    read( arg_in(1), 58) Dm
    read( arg_in(2), 58) theta0
    write( DmChar,      58 ) Dm
    write( theta0char,  58 ) theta0
    
    
    
    
    
    
    
    
    read( arg_in(3), 58) thick  ! This is material-dependent. Update this.
    write( thickness,   58 ) thick
    

    theta0 = theta0 * 4.D0 * atan(1.D0) / 180
    !theta0  = 1.0E-2
    !Dm      = 300.0E-9      ! Mass difference n' - n; in eV
    !Vopt    = 5.877E-8      ! Vopt for Cd
    !W1      = 9.914E-15     ! This is W_sc, to be multiplied by vel
    !W2      = 8.4558E-9     ! This is W_abs, measured at 2200 m/s. Constant.
    Vopt    = 1.992E-7      ! Vopt for B4C
    W1      = 2.397E-14     ! Velocity dependent part of B4C
    W2      = 6.102E-9      ! Constant part of B4C
    !thick   = 0.8e-3        
    !thick   = 3.2e-5        ! Thickness of B4C in m
    !thick   = 3.5e-6        ! Thickness of Cd in m
    A       = 0.0
    B       = 0.0
    vel     = 2318.0        ! in m/s
    numSteps= int(1e3)      ! Easier to specify a specific number of steps
    tStep   = thick / vel / numSteps
    N       = 40    ! Arbitrarily set for the moment - later going to be an
                    ! exponential distribution around some value, or, number of
                    ! neutrons to simulate

    i       = cmplx(0.0, 1.0)
    
    rho(1, 1) = 1.0
    rho(2, 1) = 0.0
    rho(1, 2) = 0.0
    rho(2, 2) = 0.0
    psi(1)          = 2*(0.5 * atan(0.5 * abs(Dm) * tan(2. * theta0) / (Vopt - Dm)))**2 !cmplx(1.0, 0.0)
    psi(2)          = 1 - 2*psi(1)!cmplx(0.0, 0.0)
    
    ! This is just for nice formatting
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

! === End variable assignments =================================================


! === Begin initial file I/O  ==================================================
    ! Write each file to a directory of the appropriate mass
    directoryName = (("./data/m/" // trim(adjustl(DmChar))))
    call system('mkdir -p '// trim(adjustl(directoryName)))

    open(1, file = trim(adjustl(directoryName)) // "/B4C_Dm-" // trim(adjustl(DmChar)) // &
        &"_theta0-" // trim(adjustl(theta0char)) // "_data_" // uuid // ".dat", &
        &status = "unknown")
    write(1, *) '#' // massChar, thetaChar
    write(1, '(2F16.11)') Dm, theta0
    write(1, 56) '#' // dist, rhonn(1, 1), rhonn(2, 1), rhonn(1, 2), rhonn(2, 2)
    write(1, 57) 0.0, rho(1, 1), rho(2, 1), rho(1, 2), rho(2, 2)

! === End initial file I/O =====================================================
    

! === Begin main loop ==========================================================
    ! The "main" loop - iterates over the material thickness
    do 100 j = 1, numSteps
        call exactBanfor(Dm, vel, theta0, Vopt, W1, W2, thick, A, B, j*tStep, psi, rho)
        
        ! These go between material layers - the "out" of one is the "in" of the next
        ! psi(1) = sqrt(rho(1, 1))
        ! psi(2) = sqrt(rho(2, 2))

        write(1, 57) j*tStep*vel, rho(1, 1), rho(2, 1), rho(1, 2), rho(2, 2)

100 end do

! === End main loop ============================================================

56  format(A19, 4A16)
57  format(F16.14, 6ES16.6E3)
    close(1)
end program d2o
