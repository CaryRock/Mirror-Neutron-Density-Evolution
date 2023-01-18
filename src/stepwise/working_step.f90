program step
    use uuid_module
    use f90getopt
    use material_list
    use exact_banfor_module
    use get_parameters
    use NNp_file_writing
    !$ use OMP_LIB
    implicit none

    ! This program steps through each material "between" each scattering event,
    ! that is, Wsc = 0, but each step size is the elastic scatt. length, and
    ! computes the probability of a neutron or mirror neutron appearing at the
    ! end of the material. At each step, the "initial" state of the particle is
    ! reset to the initial.
    ! Must simulate both neutrons and mirror neutrons to get contribution to
    ! each probability from both starting states.

    ! Custom type
    type(materiallist), dimension(:), allocatable   :: inventory
    integer, parameter                              :: prec = 18
    
!    print *, "THIS CODE IS STILL BEING DEBUGGED."

    ! Arrays
    complex, dimension(2)::psiN =(/ cmplx(1.0, 0.0), &
                                                &cmplx(0.0, 0.0) /)
    complex, dimension(2)::psiM =(/ cmplx(0.0, 0.0), &
                                                &cmplx(1.0, 0.0) /)
    complex, dimension(2):: psi, psi_n_result, psi_m_result
    !complex(8), dimension(:,:), allocatable :: prob_results
        ! Depends on the number of steps taken in the material
    real, dimension(2, 2)            :: rho
    real, dimension(:,:), allocatable:: rho_n_result, rho_m_result
    real, dimension(2, 2), parameter :: rhoN = &
                        &reshape((/1.0, 0.0, 0.0, 0.0/), shape(rho))
    real, dimension(2, 2), parameter :: rhoM = &
                        &reshape((/0.0, 0.0, 0.0, 1.0/), shape(rho))
    real, dimension(:, :), allocatable   :: meshNavg, meshMavg
    real, dimension(:, :), allocatable   :: meshNvar, meshMvar

    ! Dimension is basically nMasses rows x nAngles columns, by 
    ! 4 "layers" (<rho_11>, <rho_22>, var(rho_11), var(rho_22)) 
    real, dimension(:), allocatable  :: vel_list, Masses, Angles
    real, dimension(4)  :: O, P
    real                :: PNvels, PMvels, PN, PM

    ! Reals
    real :: dlina, Dm, theta0, vel, Vopt, Wsc, Wabs, A, B, tStep, xStep
    real :: x, lambda, Navg, Nvar, Mavg, Mvar

    ! Integers
    integer :: i, j, k, l, m, numSteps, num_lines, nMasses, nAngles, nVels
    integer :: kmin, kmax
    !$ integer  :: num_vels_to_run

    ! Logicals
    ! Some of these are only for lazy compatibility
    logical :: only_endpoint, no_scattering, no_absorption, N_initial

    ! Characters
    type(mesh_file_data)    :: mfp
    character(prec), dimension(2, 2)    :: rhonn
    character(prec) :: uni, dist, thetaChar, massChar, xChar
    character(prec) :: DmChar, theta0Char, thickness
    character(2)    :: ff58
    character(3)    :: dat_type
    character(4)    :: prog_type
    character(20)   :: f58
    character(36)   :: uuid
    character(64)   :: filename
    character(128)  :: INFILE_1, INFILE_2, INFILE_3, INFILE_4, INFILE_5
    character(256)  :: directoryName, errorlog
    character(256)  :: file1, file2!, file3, file4

! === Begin variable assignments, etc. =========================================
    xChar = "Global X"
    rhonn(1, 1) = "pnn"
    rhonn(2, 1) = "pmn"
    rhonn(1, 2) = "pnm"
    rhonn(2, 2) = "pmm"
    thetaChar = "theta"
    uni = "unitar"
    massChar = "mass"
    dist = "distance"
    uuid = generate_UUID()
    DmChar = trim(DmChar)
    theta0char = trim(theta0char)

    write(ff58, "(I2)") prec
    f58 = "(F" // trim(adjustl(ff58)) // ".11)"

    dat_type = "ave"
    prog_type= "step"
! === End variable assignments, etc. ===========================================

    !INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"
    INFILE_3 = "velocity.list"
    INFILE_4 = "deltaMs"
    INFILE_5 = "thetas"

    errorlog = "error-" // uuid // ".log"

    write( DmChar,      trim(adjustl(f58)) ) Dm
    write( theta0char,  trim(adjustl(f58)) ) theta0
    
    call get_vel_params(INFILE_2, INFILE_4, INFILE_5, inventory, psi, &
                    &filename, num_lines, only_endpoint, no_scattering, &
                    &no_absorption, Masses, Angles, nMasses, nAngles, nVels)
    nVels = get_vel_lines(INFILE_3)

! TODO: Make this either a command-line parameter or the like to be able to manage
!       this more simply.


    allocate(vel_list(nVels))
    call get_velocities(INFILE_3, nVels, vel_list)

    if (psi(1) .eq. cmplx(1.0, 0.0)) then
        N_initial = .true.
    else if (psi(2) .eq. cmplx(1.0, 0.0)) then
        N_initial = .false.
    else
        print *, "Psi    = ", psi
        print *, "Psi(1) = ", psi(1)
        print *, "Psi(2) = ", psi(2)
        print *, "Currently, this program cannot handle initial mixed states."
        print *, "Please complain to the maintainer."
        stop
    end if


    nVels = 1000! Controls how many values of the velocity to average over, not
                ! how many velocities there are total
    
    PNvels = 0.0
    PMvels = 0.0

    write( thickness, f58 ) sum(inventory%d)
    !print *, "thickness: ", thickness

    mfp = create_msd_struct(psi, filename, thickness, DmChar, theta0char, &
        &uuid, massChar, thetaChar, Dm, theta0, xChar, dist, rhonn, &
        &x, rho, no_scattering, no_absorption, prec)

    !call step_file_prepare_ave(mfp, directoryName, file1)
    call file_prepare(mfp, prog_type, dat_type, directoryName, file1, rhonn(1, 1), 1)
    call file_prepare(mfp, prog_type, dat_type, directoryName, file2, rhonn(2, 2), 2)
!   dat_type = "var"
!    call file_prepare(mfp, prog_type, "var", directoryName, file3, rhonn(1, 1), 3)
!    call file_prepare(mfp, prog_type, "var", directoryName, file4, rhonn(2, 2), 4)

    allocate(meshNavg(nMasses, nAngles))
    allocate(meshMavg(nMasses, nAngles))
!    allocate(meshNvar(nMasses, nAngles))
!    allocate(meshMvar(nMasses, nAngles))

    Wsc  = 0.0
    Vopt = inventory(1)%V    * 1.E-9
    Wabs = inventory(1)%Wabs * 1.E-9
    numSteps = inventory(1)%steps
    
    print *, ""
    print *, "Warning! This program has a fixed number of steps! It is &
        &only intended to be run for D2O at the moment!"
    print *, "Number of steps: ", numSteps
    print *, ""
    dlina = inventory(1)%d / 100    ! Distance in m

    rho = rhoN
    
    !!$set_num_threads(16)

    !$OMP PARALLEL DO Private(i, j, k, kmin, kmax, vel, Dm, theta0, &
    !$OMP& tStep, rho, psi, PN, PM, O, P, PNvels, PMvels, Navg, Mavg)
    do i = 1, nMasses
        Dm = Masses(i)

        do j = 1, nAngles
            theta0 = Angles(j)

            PNvels = 0.0
            PMvels = 0.0

            O = 0.0
            P = 0.0
            PN= 1.0
            PM= 0.0

            !kmin = 1 + (i - 1) * nAngles * nVels
            kmin = 1
            !kmax = i * nAngles * nVels
            kmax = nVels
            !if (i .eq. nMasses .and. j .eq. nAngles) kmax = kmin + nVels - 1
vels:       do k = kmin, kmax
              ! velocity in m/s
              !vel = vel_list((i - 1) * nAngles * nVels + (j - 1) * nVels + k)
              vel = vel_list(k)
              tStep = dlina / numSteps / vel

              rho = rhoM
              psi = psiM
              call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep,&
                &psi, rho)
              O(3) = rho(1, 1)
              O(4) = rho(2, 2)

              rho = rhoN
              psi = psiN
              call exactBanfor(Dm, vel, theta0, Vopt, Wsc, Wabs, tStep,&
                &psi, rho)
              O(1) = rho(1, 1)
              O(2) = rho(2, 2)

              !P(1) = PN * O(1)
              !P(2) = PN * O(2)
              !P(3) = 0.0!PM * O(3)
              !P(4) = 0.0!PM * O(4)

              PN = rho(1, 1)!O(1)
              PM = rho(2, 2)!O(2)

              if (numSteps .gt. 1) then
                PN = rho(1, 1)
                PM = rho(2, 2)

matSteps:       do l = 2, numSteps
                  !PN = rho(1, 1)
                  !PM = rho(2, 2)

                  P(1) = PN * O(1)
                  P(2) = PN * O(2)
                  P(3) = PM * O(3)
                  P(4) = PM * O(4)
                  PN = P(1) + P(3)
                  PM = P(2) + P(4)
                end do matSteps
              end if

!              PN = P(1) + P(3)
!              PM = P(2) + P(4)
              PNvels = PNvels + PN
              PMvels = PMvels + PM
            end do vels

            ! Assign to the corresponding mesh arrays
            meshNavg(i, j) = PNvels / nVels
            meshMavg(i, j) = PMvels / nVels
        end do
    end do
    !$OMP END PARALLEL DO

    ! Results need to be written out to a file, preserving structure, as well as
    ! the files for the masses and angles used

    open(unit = 1, file = trim(adjustl(directoryName)) &
                        &// "masses.txt", status = "unknown")
    open(unit = 2, file = trim(adjustl(directoryName)) &
                        &// "angles.txt", status = "unknown")
    open(unit = 3, file = trim(adjustl(directoryName)) &
                        &// "velocities.txt", status = "unknown")
    write(unit = 1, fmt = 51) "#", "Mass Delta"
    write(unit = 2, fmt = 51) "#", "Angle"
    write(unit = 3, fmt = 51) "#", "Velocities"

    do i = 1, nMasses
        write(unit = 1, fmt = 50) Masses(i)
    end do
    close(unit = 1)

    do i = 1, nAngles
        write(unit = 2, fmt = 50) Angles(i)
    end do
    close(unit = 2)

    do i = 1, nVels
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
    stop
end program step