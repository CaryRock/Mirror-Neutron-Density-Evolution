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
    !real, dimension(:), allocatable  :: PNvels, PMVels, P1, P2, P3
    !real, dimension(:), allocatable  :: P4, PN, PM
    real, dimension(4)  :: O, P
    real                :: PNvels, PMvels, PN, PM

    ! Reals
    real :: dlina, Dm, theta0, vel, Vopt, Wabs, A, B, tStep, xStep
    real :: x, lambda, Navg, Nvar, Mavg, Mvar

    ! Integers
    integer :: i, j, k, l, numSteps, num_lines, nMasses, nAngles, nVels
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
                    &no_absorption, Masses, Angles, nMasses, nAngles)
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


    nVels = 100 ! Controls how many values of the velocity to average over, not
                ! how many velocities there are total
    
    ! Allocate the P_i arrays
    !allocate(PNvels(nVels))
    !allocate(PMvels(nVels))
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

    Vopt = inventory(1)%V    * 1.E-9
    Wabs = inventory(1)%Wabs * 1.E-9
    numSteps = inventory(1)%steps!124
    
    print *, ""
    print *, "Warning! This program has a fixed number of steps! It is &
        &only intended to be run for D2O at the moment!"
    print *, "Number of steps: ", numSteps
    print *, ""
    dlina = inventory(1)%d / 100    ! Distance in m
    
    !lambda = 0.0
    !A = 0.0
    !B = 0.0
    !x = 0.0

    rho = rhoN
    
    !!$set_num_threads(16)

    !$OMP PARALLEL DO Private(i, j, k, kmin, kmax, vel, Dm, theta0, tStep, &
    !$OMP& rho, psi, PN, PM, O, P, PNvels, PMvels, Navg, Mavg)
    do i = 1, nMasses
        Dm = Masses(i)
        
        do j = 1, nAngles
            theta0 = Angles(j)
            
            !kmin = 1 + (i - 1) * nAngles * nVels
            kmin = 1
            !kmax = i * nAngles * nVels
            kmax = nVels

            PNvels = 0.0
            PMvels = 0.0
            !if (i .eq. nMasses .and. j .eq. nAngles) kmax = kmin + nVels - 1
            do k = kmin, kmax
                ! velocity in m/s
                vel = vel_list(k)
                !vel = vel_list((i - 1) * nAngles * nVels + (j - 1) * nVels + k)

                tStep= dlina / numSteps / vel  ! Time step size
                !xStep= dlina / numSteps ! Spatial step size

                ! Initial value for PN and PM need to be computed
                rho = rhoN
                psi = psiN

                call exactBanfor(Dm, vel, theta0, Vopt, 0.0, Wabs, &
                    &tStep, psi, rho)
                
                O(1) = rho(1, 1)
                O(2) = rho(2, 2)

                rho = rhoM
                psi = psiM
                call exactBanfor(Dm, vel, theta0, Vopt, 0.0, Wabs, &
                    &tStep, psi, rho)
                O(3) = rho(1, 1)
                O(4) = rho(2, 2)

!                !$OMP CRITICAL
!                if (i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) then
!                  print *, "Dm, theta0: ", Dm, theta0
!                  print *, "vel: ", vel
!                  print *, "Vopt, Wabs: ", Vopt, Wabs
!                  print *, "rho: ", rho(1, 1), rho(2, 2)
!                  print *, "abs: ", 1 - rho(1, 1) - rho(2, 2)
!                  print *, ""
!                end if
!
!                if (i .eq. 1 .and. k .eq. 1) then
!                  print *, "theta, rho: ", theta0, rho(1, 1), rho(2, 2)
!                  print *, "abs: ", 1 - rho(1, 1) - rho(2, 2)
!                  print *, ""
!                end if
!                !$OMP END CRITICAL

                PN = rho(1, 1)
                PM = rho(2, 2)

                ! Call a subroutine to handle the details of computing the 
                ! probabilities

                if (numSteps .gt. 1) then
                  call compute_probabilities(rho, rhoN, rhoM, psi, psiN, &
                    &psiM, Dm, theta0, vel, tStep, Vopt, Wabs, &
                    &numSteps, O, P, PN, PM, N_initial)
                    !&numSteps, P1, P2, P3, P4, PN, PM, N_initial)
                end if

                ! Store computed values in a holding variable/array
                PNvels = PNvels + PN!PN(numSteps)
                PMvels = PMvels + PM!PM(numSteps)
            end do

            ! Compute the average and variance in the PNvels and PMvels
            Navg = PNvels / nVels!sum(PNvels) / nVels
            Mavg = PMvels / nVels!sum(PMvels) / nVels
            !Nvar = compute_variance(PNvels)
            !Mvar = compute_variance(PMvels)

            ! These values should be assigned to corresponding mesh arrays
            meshNavg(i, j) = Navg
            meshMavg(i, j) = Mavg
            !meshNvar(i, j) = Nvar
            !meshMvar(i, j) = Mvar

!            k = k + 1
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

    do i = 1, size(vel_list)
        !write(unit = 3, fmt = 50) vel_list(i) * 100.
    end do
    close(unit = 3)
    
    do i = 1, nMasses
        write(unit = 10, fmt = 52) meshNavg(i, :)
        write(unit = 20, fmt = 52) meshMavg(i, :)
!        write(unit = 30, fmt = 52) meshNvar(i, :)
!        write(unit = 40, fmt = 52) meshMvar(i, :)
    end do

    close(unit = 10)
    close(unit = 20)
!    close(unit = 30)
!    close(unit = 40)


50  format(ES17.8E3)
51  format(A1, A15)
52  format(1608ES17.8E3)
!52  format(2941ES17.8E3)

contains
    subroutine compute_probabilities(rrho, rrhoN, rrhoM, ppsi, ppsiN, &
        &ppsiM, DDm, ttheta0, vvel, ttStep, VVopt, WWabs, nnumSteps, &
        &OO, PP, PPN, PPM, reset)!&PP1, PP2, PP3, PP4, PPN, PPM, reset)
        implicit none
        real, dimension(2, 2), intent(in)    :: rrhoN, rrhoM
        real, dimension(2, 2)                :: rrho
        complex, dimension(2), intent(in)    :: ppsiN, ppsiM
        complex, dimension(2)                :: ppsi
        real, intent(in)     :: DDm, ttheta0, vvel, ttStep, VVopt, WWabs
        real, dimension(4)   :: OO, PP
        real                 :: PPN, PPM
        !real, dimension(:)   :: PP1, PP2, PP3, PP4, PPN, PPM
        integer, intent(in)  :: nnumSteps

        logical, intent(in)  :: reset
        real, dimension(2, 2):: resetRho
        complex, dimension(2):: resetPsi

        !real                 :: AA, BB, llambda
        integer              :: ii

        ! Reset rho and psi, just to make sure
        if (reset .eqv. .true.) then
            resetRho = rrhoN
            resetPsi = ppsiN
        else
            resetRho = rrhoM
            resetPsi = ppsiM
        end if

        rrho = resetRho
        ppsi = resetPsi

        ! Obtain first element of PN and PM
        !PPN(1) = rrho(1, 1)
        !PPM(1) = rrho(2, 2)
        PPN = rrho(1, 1)
        PPM = rrho(2, 2)

        do 101 ii = 2, nnumSteps
            ! Reset rho and psi for use
            !rrho = resetRho
            !ppsi = resetPsi

            ! Iterate over the probability arrays
            !call exactBanfor(DDm, vvel, ttheta0, VVopt, 0.0, WWabs, &
            !    &ttStep, ppsi, rrho)

            !PP1(ii) = PPN(ii - 1) * rrho(1, 1)
            !PP2(ii) = PPN(ii - 1) * rrho(2, 2)
            !PP3(ii) = PPM(ii - 1) * rrho(1, 1)
            !PP4(ii) = PPM(ii - 1) * rrho(2, 2)
            PP(1) = PPN * OO(1)!rrho(1, 1)
            PP(2) = PPN * OO(2)!rrho(2, 2)
            PP(3) = PPM * OO(3)!rrho(1, 1)
            PP(4) = PPM * OO(4)!rrho(2, 2)

            !PPN(ii) = PP1(ii) + PP3(ii)
            !PPM(ii) = PP2(ii) + PP4(ii)
            PPN = PP(1) + PP(3)
            PPM = PP(2) + PP(4)

            !print *, "PPN(ii - 1): ", PPN(ii - 1)
            !print *, "    PPN(II): ", PPN(ii)
101     end do
        return
    end subroutine compute_probabilities

    real function compute_variance(array) result(var)
        implicit none
        real, dimension(:), intent(in)   :: array
        real :: ave
        integer :: s, m

        ave = 0.D0
        var = 0.D0
        s = size(array)

        ave = sum(array) / s

        do m = 1, s
            var = var + (array(m) - ave)**2
        end do

        var = var / (s - 1)
        return
    end function compute_variance
end program step
