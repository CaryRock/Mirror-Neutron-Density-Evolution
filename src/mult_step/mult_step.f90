program mult_step
    use uuid_module
    use f90getopt
    use material_list
    use exact_banfor_module
    use get_parameters
    use NNp_file_writing
    use distributions
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
    real, dimension(4)                  :: O, P

    ! Scalars
    real    :: PN, PM
    real    :: dlina, Dm, theta0, vel, Vopt, Wabs, Wsc, A, B, tStep, xStep
    real    :: x, lambda, Navg, Nvar, Mavg, Mvar, PNvels, PMvels, vAvg, temp

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
    prog_type = "mult_step"

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

    x   = 0.0
    rho = rhoReset
    vAvg = 0.0

    header(1) = "Core"
    header(2) = "Material"
    header(3) = "Material Thickness"
    header(4) = "Global Distance"
    header(5) = "PN"
    header(6) = "PM"

    Wsc = 0.0
    temp = 342.0

!          # Mtrl MtrlThck G.Dist PN PM
    print '(A2, A35, A28, A30, A30, A27)', "# ", header(2), header(3), header(4), header(5), header(6)
    print *, adjustl(header(1)), ",", 20.0, ",", 20.0,",", 1.0,",", 0.0
    !$OMP PARALLEL DO Private(i, j, k, l, m, kmin, kmax, vel, Dm, theta0, &
    !$OMP& tStep, rho, psi, O, P, PN, PM, PNvels, PMvels, Navg, Mavg, &
    !$OMP& dlina, Vopt, Wabs, numSteps, x, vAvg, vel_list)
massLoop:   do i = 1, nMasses
        Dm = Masses(i)

      ! Since the program is split along the masses (for the CPU case),
      ! there are nMasses number of RNGs. Currently handled by GCC.
anglLoop:   do j = 1, nAngles
            theta0  = Angles(j)
            
            ! Generate the velocities to use
            do k = 1, nVels
              call YK_MAXW(temp, vel_list(k))
            end do
            
            vAvg = vAvg + sum(vel_list)
            vAvg = vAvg / nVels
            
            !O       = 0.0
            !P       = 0.0
            !PN      = 1.0
            !PM      = 0.0
            PNvels  = 0.0
            PMvels  = 0.0

velsLoop:   do k = 1, nVels!kmin, kmax
                PN = 1.0
                PM = 0.0

                x = 20.0  ! Global Coordinate, cm

                vel = vel_list(k) ! m/s
                rho = rhoReset

mtrlLoop:       do l = 1, nMaterials
                  numSteps = inventory(l)%steps
                  !dlina = inventory(l)%elscatl * numSteps / 100. ! m
                  tStep = inventory(l)%elscatl / 100. / vel    ! s
                  Vopt = inventory(l)%V     * 1.E-9 ! eV
                  Wabs = inventory(l)%Wabs  * 1.E-9 ! eV
                    
! Set-up probabilities for initially starting as a neutron
                  rho = rhoReset
                  call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, &
                    &tStep, psi, rho)
                  O(1) = rho(1, 1)
                  O(2) = rho(1, 2)
                  O(3) = rho(2, 1)
                  O(4) = rho(2, 2)

                  if (l .eq. 1) then
                    if (N_initial) then
                      P(1) = PN * O(1)
                      P(2) = PN * O(2)
                      P(3) = 0.0!PM * O(3)
                      P(4) = 0.0!PM * O(4)
                    else
                      P(1) = 0.0
                      P(2) = 0.0
                      P(3) = PM * O(3)
                      P(4) = PM * O(4)
                    end if
                  !else
                    !P(1) = PN * O(1)
                    !P(2) = PN * O(2)
                    !P(3) = PM * O(3)
                    !P(4) = PM * O(4)
                  end if
                  
                  !!$OMP CRITICAL
                  !print *, "theta0  : ", theta0
                  !print *, "dm      : ", Dm
                  !print *, "vel     : ", vel
                  !print *, "tStep   : ", tStep
                  !print *, "O       : ", O
                  !print *, "rhoReset: ", rhoReset
                  !print *, ""
                  !stop
                  !!$OMP END CRITICAL

matSteps:         do m = 1, numSteps
                    PN = P(1) + P(3)
                    PM = P(2) + P(4)

                    P(1) = PN * O(1)
                    P(2) = PN * O(2)
                    P(3) = PM * O(3)
                    P(4) = PM * O(4)
                  end do matSteps

!
                  !!$OMP CRITICAL
                  !if (PMvels + PM .lt. 0) then
                  !  print *, ""
                  !  print *, "Dumping: "
                  !  print *, "dM (eV), theta (rad), vel (m/s), numSteps: ", Dm, theta0, vel, numSteps
                  !  print *, "# velocities averaged over: ", nVels
                  !  print *, "Average velocity: ", vAvg
                  !  print *, "Particle path length: ", inventory(l)%elscatl * inventory(l)%steps
                  !  print *, "i, j, k: ", i, j, k
                  !  print *, "vel_list: ", vel_list
                  !  print *, "O: ", O
                  !  print *, "P: ", P
                  !  print *, "PN, PM: ", PN, PM
                  !  print *, "PNvels, PMvels: ", PNvels, PMvels
                  !  print *, ""
                  !  stop
                  !end if
!
                  ! angles(101) = 0.001, mass(47) = 199.5...E-9 eV, mass(55) = 501.2...E-9 eV
                  !if (i .eq. nMasses .and. j .eq. 101 .and. k .eq. nVels ) then
                  if (j .eq. 101) then
                    if (k .eq. nVels) then
                      if (i .eq. 47) then
                        x = x + inventory(l)%d
                        print *, inventory(l)%matName, ",",inventory(l)%d, ",",&
                          &x,",", (PNvels + PN)/nVels,",", (PMvels + PM)/nVels
                        write(11, *) inventory(l)%matName, ",",inventory(l)%d, &
                          &",", x,",", (PNvels + PN)/nVels,",", (PMvels + PM)/nVels
                        if (l .eq. nMaterials) then
                          print *, ""
                          write(11, *) ""
                          print *, "dM (eV), theta (rad), vel (m/s), numSteps: ",&
                            &Dm, theta0, vel, numSteps
                          write(11, *) "dM (eV), theta (rad), vel (m/s), numSteps: ", &
                            &Dm, theta0, vel, numSteps
                          print *, "# velocities averaged over: ", nVels
                          write(11, *) "# velocities averaged over: ", nVels
                          print *, "Average velocity: ", vAvg
                          write(11, *) "Average velocity: ", vAvg
                          print *, "Particle path length: ", inventory(l)%elscatl * inventory(l)%steps
                          write(11, *) "Particle path length: ", inventory(l)%elscatl * inventory(l)%steps
                          print *, "O: ", O
                          write(11, *) "O: ", O
                          print *, "numSteps: ", numSteps
                          write(11, *) "numSteps: ", numSteps
                          print *, ""
                        end if
                      end if

                      if (i .eq. 55) then
                        x = x + inventory(l)%d
                        print *, inventory(l)%matName, ",",inventory(l)%d, ",", &
                          &x,",", (PNvels + PN)/nVels,",", (PMvels + PM)/nVels
                        write(12, *) inventory(l)%matName, ",",inventory(l)%d, &
                          &",", x,",", (PNvels + PN)/nVels,",", (PMvels + PM)/nVels
                        if (l .eq. nMaterials) then
                          print *, ""
                          write(12, *) ""
                          print *, "dM (eV), theta (rad), vel (m/s), numSteps: ", &
                            &Dm, theta0, vel, numSteps
                          write(12, *) "dM (eV), theta (rad), vel (m/s), numSteps: ", &
                            &Dm, theta0, vel, numSteps
                          print *, "# velocities averaged over: ", nVels
                          write(12, *) "# velocities averaged over: ", nVels
                          print *, "Average velocity: ", vAvg
                          write(12, *) "Average velocity: ", vAvg
                          print *, "Particle path length: ", inventory(l)%elscatl * inventory(l)%steps
                          write(12, *) "Particle path length: ", inventory(l)%elscatl * inventory(l)%steps
                          print *, "O: ", O
                          write(12, *) "O: ", O
                          print *, "numSteps: ", numSteps
                          write(12, *) "numSteps: ", numSteps
                          print *, ""
                        end if
                      end if
                    end if
                  end if
!                  !$OMP END CRITICAL

                end do mtrlLoop

                PNvels = PNvels + PN
                PMvels = PMvels + PM
            end do velsLoop

            meshNavg(i, j) = PNvels / float(nVels)
            meshMavg(i, j) = PMvels / float(nVels)
        end do anglLoop
    end do massLoop
    !$OMP END PARALLEL DO

! === Begin File Writing and Output ============================================

    ! Generate a random velocity list
    do k = 1, nVels
      call YK_MAXW(temp, vel_list(k))
    end do

    open(unit = 1, file = trim(adjustl(directoryName)) &
        &// "masses.txt", status = "unknown")
    open(unit = 2, file = trim(adjustl(directoryName)) &
        &// "angles.txt", status = "unknown")
    open(unit = 3, file = trim(adjustl(directoryName)) &
        &// "velocities.txt", status = "unknown")
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

    do i = 1, nMasses
        write(unit = 10, fmt = 52) meshNavg(i, :)
        write(unit = 20, fmt = 52) meshMavg(i, :)
    end do
    close(unit = 10)
    close(unit = 20)

50  format(ES17.8E3)
51  format(A1, A15)
52  format(1608ES17.8E3)
53  format(A1, A17)
    
! === End File Writing and Output ==============================================

! === End Preparation/Instantiations ===========================================

! ==============================================================================
end program mult_step