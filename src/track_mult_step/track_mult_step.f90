program track_mult_step
    use uuid_module
    use f90getopt
    use material_list
    use exact_banfor_module
    use get_parameters
    use NNp_file_writing
    use distributions
    use NNp_Loop
    !$ use OMP_LIB
    implicit none

! === Purpose ==================================================================

    ! This program is meant to step through both materials and layers in 
    ! increments of the material's scattering length as a modification of 
    ! the density function approach to mimic the Lindblad formalism. That is,
    ! W_sc = 0, but the spatial step size through the material is that of the 
    ! scattering length.

    ! The end goal is to simulate the passing of neutrons -> mirror neutrons
    ! through a material, as well as being initially for the checking of the
    ! STEREO collaboration's results.

! ==============================================================================

! === Begin Variable Initialization ============================================

    type(materiallist), dimension(:), allocatable    :: inventory
    integer, parameter                              :: prec = 18

    ! Arrays
    complex, dimension(2), parameter    ::psiN = &
        &(/ cmplx(1.0, 0.0), cmplx(0.0, 0.0) /)
!    complex, dimension(2), parameter    ::psiM = &
!        &(/ cmplx(0.0, 0.0), cmplx(1.0, 0.0) /)
    complex, dimension(2):: psi

    real, dimension(2, 2)               :: rho, rhoReset
    real, dimension(:,:),allocatable    :: rho_n_result, rho_m_result
    real, dimension(2, 2), parameter    :: rhoN = &
        &reshape((/1.0, 0.0, 0.0, 0.0/), shape(rho))
    real, dimension(2, 2), parameter    :: rhoM = &
        &reshape((/0.0, 0.0, 0.0, 1.0/), shape(rho))
    real, dimension(:, :), allocatable  :: meshNavg, meshMavg
    real, dimension(:, :), allocatable  :: meshNvar, meshMvar
    
    real, dimension(:), allocatable     :: vel_list, Masses, Angles
    real, dimension(:), allocatable     :: matAveN, matAveM
    real, dimension(4)                  :: O, P

    ! Scalars
    real    :: PN, PM
    real    :: dlina, Dm, theta0, vel, Vopt, Wabs, Wsc, A, B, tStep, xStep
    real    :: x, lambda, Navg, Nvar, Mavg, Mvar, PNvels, PMvels, vAvg, temp
    ! JANK
    real    :: remThick, matD, elLambda, dx

    integer :: i, j, k, l, m, numSteps, nMaterials, nMasses, nAngles, nVels
    integer :: kmin, kmax, totSteps, accumulatedSteps

    ! Some of these are only for lazy compatibility
    logical :: only_endpoint, no_scattering, no_absorption, N_initial

    ! Characters and strings
    type(mesh_file_data)                :: mfp
    character(prec), dimension(2, 2)    :: rhonn
    character(prec) :: uni, dist, thetaChar, massChar, xChar, DmChar
    character(prec) :: theta0Char, thickness
    character(2)    :: ff58
    character(3)    :: dat_type
    character(128)  :: prog_type
    character(20)   :: f58
    character(36)   :: uuid
    character(32), dimension(6) :: header
    !character(32)   :: core
    character(64)   :: filename
    character(128)  :: INFILE_1, INFILE_2, INFILE_3, INFILE_4, INFILE_5
    character(256)  :: directoryName, errorlog, file1, file2!, file3, file4

! === End Variable Initialization ==============================================

! === Begin Variable Assignment ================================================
    xChar = "Gloabl X"
    rhonn(1, 1) = "pnn"
    rhonn(2, 1) = "pmn"
    rhonn(1, 2) = "pnm"
    rhonn(2, 2) = "pmm"
    thetaChar = "theta"
    uni = "unitar"
    massChar = "mass"
    dist = "distance"
    uuid = generate_UUID()
    DmChar = trim(DmChar)   ! TODO: Is this necessary? Why did I do this?
    theta0Char = trim(theta0Char)! TODO: Same here

    write(ff58, "(I2)") prec
    f58 = "(F" // trim(adjustl(ff58)) // ".11)"

    dat_type  = "ave"
    prog_type = "track_mult_step"

    !INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"
    INFILE_3 = "velocity.list"
    !INFILE_4 = "deltaMs"
    !INFILE_5 = "thetas"

    errorlog = "error-" // uuid // ".log"

    write( DMChar,      trim(adjustl(f58)) ) Dm
    write( theta0Char,  trim(adjustl(f58)) ) theta0

! === End Variable Assignment ==================================================

! === Begin Preparation/Instantiations =========================================

    call get_vel_params(INFILE_2, INFILE_4, INFILE_5, inventory, psi, &
        &filename, nMaterials, only_endpoint, no_scattering, &
        &no_absorption, Masses, Angles, nMasses, nAngles, nVels)

    !nVels = get_vel_lines(INFILE_3)

    !allocate(vel_list(nVels))
    !call get_velocities(INFILE_3, nVels, vel_list)

    if (psi(1) .eq. cmplx(1.0, 0.0)) then
        N_initial = .true.
        rhoReset = rhoN
    else if (psi(2) .eq. cmplx(1.0, 0.0)) then
        N_initial = .false.
        rhoReset = rhoM
    else
        print *, "Psi    = ", psi
        print *, "Currently, this program cannot handle initial mixed states."
        print *, "Feel free to complain to the maintainer."
        stop
    end if

    ! Note: this counts how many values of the velocity to average over, not how
    ! many velocities there are in total - that is, for averaging, N = nVels
    !nVels = 100  ! TODO: Make this an option later
    allocate(vel_list(nVels))

    PNvels = 0.0
    PMvels = 0.0

    write( thickness, f58 ) sum(inventory%d)

    mfp = create_msd_struct(psi, filename, thickness, DmChar, &
    &theta0char, uuid, massChar, thetaChar, Dm, theta0, xChar, dist, &
    &rhonn, x, rho, no_scattering, no_absorption, prec)

    call file_prepare(mfp, trim(adjustl(prog_type)), dat_type, &
        &directoryName, file1, rhonn(1, 1), 1)
    call file_prepare(mfp, trim(adjustl(prog_type)), dat_type, &
        &directoryname, file2, rhonn(2, 2), 2)

    allocate(meshNavg(nMasses, nAngles))
    allocate(meshMavg(nMasses, nAngles))

    totSteps = sum(inventory%steps)

    allocate(matAveN(nMaterials))
    allocate(matAveM(nMaterials))

    matAveN = 0.
    matAveM = 0.

    x   = 0.0
    rho = rhoReset
    vAvg = 0.0

    header(1) = "Core"
    header(2) = "Material"
    header(3) = "Material Thickness"
    header(4) = "Global Distance"
    header(5) = "PN"
    header(6) = "PM"

    !Wsc = 0.0
    temp = 342.0  ! This can be read in from the material list cards

!          # Mtrl MtrlThck G.Dist PN PM

    print '(A2, A35, A28, A30, A30, A27)', "# ", header(2), header(3), &
      &header(4), header(5), header(6)
    print *, adjustl(header(1)), ",", 20.0, ",", 20.0,",", 1.0,",", 0.0

    !$OMP PARALLEL DO Private(i, j, k, l, m, kmin, kmax, vel, Dm, theta0, &
    !$OMP& tStep, rho, psi, O, P, PN, PM, PNvels, PMvels, Navg, Mavg, &
    !$OMP& dlina, Vopt, Wabs, Wsc, numSteps, x, vAvg, vel_list, &
    !$OMP& matAveN, matAveM)
massLoop:   do i = 1, nMasses
        Dm = Masses(i)

      ! Since the program is split along the masses (for the CPU case),
      ! there are nMasses number of RNGs. Currently handled by GCC.
anglLoop:   do j = 1, nAngles
            theta0  = Angles(j)
            
            ! Generate the velocities to use
            do k = 1, nVels
              call neutron_MW_dist(temp, vel_list(k), 0., 7153.416954)
            end do
            
            vAvg = vAvg + sum(vel_list)
            vAvg = vAvg / nVels
            
            !O       = 0.0
            !P       = 0.0
            !PN      = 1.0
            !PM      = 0.0
            PNvels  = 0.0
            PMvels  = 0.0

            ! Generate all the velocities to be used
            !do k = 1, nVels
            !  call neutron_MW_dist(temp, vel_list(k), 0., 7153.416954)
            !end do

            ! Can this be made more abstract? As in, some common 
            ! interface that will handle the passing of the variables,
            ! but not using a COMMON block?
            ! Lookup: interface, abstract interface, procedure,
            !         function pointers
            call main_vel_loop(inventory, rhoReset, rho, O, P, &
              &vel_list, PNvels, PMvels, Dm, theta0, nMaterials,&
              &nVels, header, N_initial, matAveN, matAveM, no_scattering)
            
            meshNavg(i, j) = PNvels / float(nVels)
            meshMavg(i, j) = PMvels / float(nVels)
        end do anglLoop
    end do massLoop
    !$OMP END PARALLEL DO

! === Begin File Writing and Output ============================================

    ! Generate a random velocity list
    do k = 1, nVels
      call neutron_MW_dist(temp, vel_list(k), 0., 7153.416954)
    end do

    open(unit = 1, file = trim(adjustl(directoryName)) &
        &// "masses.txt", status = "unknown")
    open(unit = 2, file = trim(adjustl(directoryName)) &
        &// "angles.txt", status = "unknown")
    open(unit = 3, file = trim(adjustl(directoryName)) &
        &// "velocities.txt", status = "unknown")
    open(unit = 4, file = trim(adjustl(directoryName)) &
        &// "log.txt", status = "unknown")
    write(unit = 1, fmt = 51) "#", "Mass Delta"
    write(unit = 2, fmt = 51) "#", "Angle"
    write(unit = 3, fmt = 53) "#", "Velocities (cm/s)"

    do i = 1, nMasses
        write(unit = 1, fmt = 50) Masses(i)
    end do
    close(unit = 1)

    do i = 1, nAngles
        write(unit = 2, fmt = 50) Angles(i)
    end do
    close(unit = 2)

    do i = 1, size(vel_list)
        write(unit = 3, fmt = 50) vel_list(i)
    end do
    close(unit = 3)
  
    !write(unit = 4, fmt = '(A32, )') asdf
    write(unit = 4, fmt = 54) header(1), ",", 20., &
      &",", 1.0, ",", 0.0

    do i = 1, nMasses
        write(unit = 10, fmt = 52) meshNavg(i, :)
        write(unit = 20, fmt = 52) meshMavg(i, :)
    end do
    close(unit = 10)
    close(unit = 20)

! === BEGIN JANK ===============================================================
    ! Write a function that'll do this for all the desired 
    ! For outputting the desired probabilities as a function of mass and angle
    ! and material while stepping through the code 
    i = 1
    j = 1
    Dm = Masses(i)
    theta0 = Angles(j)
    !call get_log_data(inventory, rhoReset, rho, O, P, &
    !&vel_list, PNvels, PMvels, Dm, theta0, nMaterials, nVels, &
    !&header, N_initial, matAveN, matAveM)
    O = 0.
    P = 0.
    vAvg = vAvg + sum(vel_list)
    vAvg = vAvg / nVels
    
    O       = 0.0
    P       = 0.0
    PN      = 1.0
    PM      = 0.0
    PNvels  = 0.0
    PMvels  = 0.0

vl: do k = 1, nVels
      if (N_initial) then
        PN = 1.0
        PM = 0.0
      else
        PN = 0.0
        PM = 1.0
      end if

      vel = vel_list(k)
      rho = rhoReset

nMat: do l = 1, nMaterials
        matAveN(l) = 0.0
        matAveM(l) = 0.0
        matD = inventory(l)%d           ! m
        numSteps = inventory(l)%steps   !
        Vopt = inventory(l)%V           ! eV
        Wabs = inventory(l)%Wabs        ! eV
        elLambda = inventory(l)%elscatl ! m

        remThick = matD
        m = 1

        if (N_initial) then
          PN = 1.0
          PM = 0.0
        else
          PN = 0.0
          PM = 1.0
        end if

mStep:  do while(remThick .gt. 0)
          rho = rhoReset
          call scttStep(elLambda, remThick, dx)
          tStep = dx / vel
          call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep, rho)
          O(1) = rho(1, 1)
          O(2) = rho(1, 2)
          O(3) = rho(2, 1)
          O(4) = rho(2, 2)

          if (l .eq. 1) then
            if (N_initial) then
              P(1) = PN * O(1)
              P(2) = PN * O(2)
              P(3) = 0.0
              P(4) = 0.0
            else
              P(1) = 0.0
              P(2) = 0.0
              P(3) = PM * O(3)
              P(4) = PM * O(4)
            end if
          end if

          PN = P(1) + P(3)
          PM = P(2) + P(4)

          P(1) = PN * O(1)
          P(2) = PN * O(2)
          P(3) = PM * O(3)
          P(4) = PM * O(4)

          m = m + 1
        end do mStep  ! matSteps

        ! For later outputting the desired probabilities as a function of material
        ! for a given mass and angle.
        matAveN(l) = matAveN(l) + PN
        matAveM(l) = matAveM(l) + PM
        print *, matAveN(l), matAveM(l)
      end do nMat     ! nMaterials

      PNvels = PNvels + PN
      PMvels = PMvels + PM
    end do vl   ! nVels
    matAveN = matAveN / nVels
    matAveM = matAveM / nVels


! === END JANK =================================================================
    do i = 1, nMaterials
        write(unit = 4, fmt = 54) inventory(i)%matName, ",", inventory(i)%d, &
          &",", matAveN(i), ",", matAveM(i)
    end do

    ! Since I'm basically making a csv, I know this can be done much more simply
    write(unit = 4, fmt = '(A4, A1, ES17.8E3, A1, A5, A1, ES17.8E3, A1, A5, A1, ES17.8E3)') &
      &"mass", ",", Dm, ",", "theta", ",", theta0, ",", "vel", ",", vel_list(1)
    close(unit = 4)


50  format(ES17.8E3)
51  format(A1, A15)
52  format(1608ES17.8E3)
53  format(A1, A17)
54  format(A32, A1, ES17.8E3, A1, ES17.8E3, A1, ES17.8E3)    
! === End File Writing and Output ==============================================

! === End Preparation/Instantiations ===========================================

! ==============================================================================
end program track_mult_step