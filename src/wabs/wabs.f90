program wabs_step
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
    !
    ! In addition to the above, it reads in a list of W_abs values from a given
    ! input file and repeats the above "step" process for each value of W_abs.
    ! This can definitely benefit from parallel processing.

    ! Custom type
    type(materiallist), dimension(:), allocatable   :: inventory
    integer, parameter                              :: prec = 18
    
!    print *, "THIS CODE IS STILL BEING DEBUGGED."

    ! Arrays
    complex(8), dimension(2)::psiN =(/ cmplx(1.0, 0.0, 8), cmplx(0.0, 0.0, 8) /)
    complex(8), dimension(2)::psiM =(/ cmplx(0.0, 0.0, 8), cmplx(1.0, 0.0, 8) /)
    complex(8), dimension(2)            :: psi, psi_n_result, psi_m_result
    !complex(8), dimension(:,:), allocatable :: prob_results
        ! Depends on the number of steps taken in the material
    real(8), dimension(2, 2)            :: rho
    real(8), dimension(:,:), allocatable:: rho_n_result, rho_m_result
    real(8), dimension(2, 2), parameter :: rhoN = reshape((/1.D0, 0.D0, 0.D0, 0.D0/), shape(rho))
    real(8), dimension(2, 2), parameter :: rhoM = reshape((/0.D0, 0.D0, 0.D0, 1.D0/), shape(rho))
    real(8), dimension(:, :), allocatable   :: meshNavg, meshMavg
    real(8), dimension(:, :), allocatable   :: meshNvar, meshMvar

    ! Dimension is basically nMasses rows x nAngles columns, by 
    ! 4 "layers" (<rho_11>, <rho_22>, var(rho_11), var(rho_22)) 
    real(8), dimension(:), allocatable  :: vel_list, Masses, Angles, wabs_list
    real(8), dimension(:), allocatable  :: PNvels, PMVels, P1, P2, P3, P4, PN, PM

    ! Reals
    real(8) :: dlina, Dm, theta0, vel, Vopt, Wsc, Wabs, A, B, tStep, xStep, x
    real(8) :: lambda, Navg, Nvar, Mavg, Mvar

    ! Integers
    integer :: i, j, k, l, w, id, Mw, Nw, numSteps, num_lines
    integer :: nMasses, nAngles, nVels, nWabs
    integer :: thread_id

    ! Logicals
    ! Some of these are only for lazy compatibility
    logical :: only_endpoint, no_scattering, no_absorption, N_initial

    ! Characters
    type(mesh_file_data)    :: mfp
    character(prec), dimension(2, 2)    :: rhonn
    character(256), dimension(:), allocatable   :: wabsFileList
    character(prec) :: uni, dist, thetaChar, massChar, xChar
    character(prec) :: DmChar, theta0Char, thickness
    character(2)    :: ff58
    character(3)    :: dat_type
    character(4)    :: prog_type
    character(20)   :: f58, f50, f51, f52
    character(36)   :: uuid
    character(64)   :: filename
    character(128)  :: INFILE_1, INFILE_2, INFILE_3
    character(128)  :: INFILE_4, INFILE_5, INFILE_6
    character(256)  :: directoryName, errorlog
    character(256)  :: file1, string1!, file3, file4

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

!50  format(ES17.8E3)
!51  format(A1, A15)
!52  format(1608ES17.8E3)

    dat_type = "ave"
    prog_type= "wabs"
! === End variable assignments, etc. ===========================================

    !INFILE_1 = "parameter.list"
    INFILE_2 = "material.list"
    INFILE_3 = "velocity.list"
    INFILE_4 = "deltaMs"
    INFILE_5 = "thetas"
    INFILE_6 = "wabs"

    errorlog = "error-" // uuid // ".log"

    write( DmChar,      trim(adjustl(f58)) ) Dm
    write( theta0char,  trim(adjustl(f58)) ) theta0
    
    call get_vel_params(INFILE_2, INFILE_4, INFILE_5, inventory, psi, &
                    &filename, num_lines, only_endpoint, no_scattering, &
                    &no_absorption, Masses, Angles, nMasses, nAngles)
    nVels = get_vel_lines(INFILE_3)
    allocate(vel_list(nVels))
    call get_velocities(INFILE_3, nVels, vel_list)


    ! Get the wabs_list values - this should be shuffled into the 
    ! get_parameters_module.f90 code
    nWabs = get_lines(INFILE_6) + 2
    allocate(wabs_list(nWabs))
    allocate(wabsFileList(2*nMasses))
    open(unit = 1, file = INFILE_6, status = "old")
    do w = 1, nWabs
        read(1, '(F37.30)') wabs_list(w)
    end do
    close(unit = 1)

    if (psi(1) .eq. cmplx(1.D0, 0.D0, 8)) then
        N_initial = .true.
    else if (psi(2) .eq. cmplx(1.D0, 0.D0, 8)) then
        N_initial = .false.
    else
        print *, "Psi    = ", psi
        print *, "Psi(1) = ", psi(1)
        print *, "Psi(2) = ", psi(2)
        print *, "Currently, this program cannot handle initial mixed states."
        print *, "Please complain to the maintainer."
        stop
    end if

    ! Allocate the P_i arrays
    allocate(PNvels(nVels))
    allocate(PMvels(nVels))

    write( thickness, f58 ) sum(inventory%d)
    !print *, "thickness: ", thickness

    mfp = create_msd_struct(psi, filename, thickness, DmChar, theta0char, &
        &uuid, massChar, thetaChar, Dm, theta0, xChar, dist, rhonn, &
        &x, rho, no_scattering, no_absorption, prec)

! TODO: UPDATE FILE PREPARATION TO HANDLE WABS COMPUTATIONS - INDEX BY M-INDEX
    do w = 1, nMasses
        call wabs_file_prepare_ave(mfp, directoryName, &
            &wabsFileList(2*(w - 1) + 1), rhonn(1, 1), 1, w - 1)
        call wabs_file_prepare_ave(mfp, directoryName, &
            &wabsFileList(2*(w - 1) + 2), rhonn(2, 2), 2, w - 1)
    end do

    allocate(meshNavg(nWabs, nAngles))
    allocate(meshMavg(nWabs, nAngles))
    !allocate(meshNavg(nMasses, nAngles))
    !allocate(meshMavg(nMasses, nAngles))
!    allocate(meshNvar(nMasses, nAngles))
!    allocate(meshMvar(nMasses, nAngles))

    Vopt = inventory(1)%V    * 1.D-9
    !Wabs = inventory(1)%Wabs * 1.D-9

c     This code currently does not handle numSteps = 1 since the method of
c     stepping requires numSteps >= 2; it should be able to handle numSteps = 1
c     in principle, but the current implementation does not.

c     TODO: IMPLEMENT ACCEPTING VALUE OF numStep = 1 FOR COMPUTATION - SHOULD BE
c     FUNCTIONALLY IDENTICAL TO SETTING OF "-S" FLAG IN MESH CODE
    numSteps = 170
c     If numSteps = 1, then basically run the velmesh algorithm
    !numSteps = 1
    numSteps = 2


    ! TODO: Consider: reading numSteps from the same file as scattering length
    ! is read from, which should be the same file as Vopt is read from.
    allocate(P1(numSteps))
    allocate(P2(numSteps))
    allocate(P3(numSteps))
    allocate(P4(numSteps))
    allocate(PN(numSteps))
    allocate(PM(numSteps))
    !print *, ""
    !print *, "Warning! This program has a fixed number of steps of 170! It is &
    !    &only intended to be run for D2O at the moment!"
    !print *, ""
    dlina = inventory(1)%d / 100
    !dlina = 2.2D-2
    dlina = 4.4D-2
    print *, "Thickness: ", dlina

    lambda = 0.D0
    A = 0.D0
    B = 0.D0
    x = 0.D0
    Wsc = 0.D0
    rho = rhoN
    
    id = 10 ! If paralell, this will be overwritten later. If not, this is fine
    print *, "nMasses = ", nMasses
    print *, "nAngles = ", nAngles
    print *, "nVf     = ", nWabs

    !$OMP PARALLEL DO PRIVATE(l)
    do l = 1, nMasses
        Dm = Masses(l)
        ! Just a big, fat, slow reset of everything
        rho = mfp%rho
        psi = mfp%psi
        P1 = 0.D0
        P2 = 0.D0
        P3 = 0.D0
        P4 = 0.D0
        PN = 0.D0
        PM = 0.D0
        PNvels = 0.D0
        PMvels = 0.D0
        Navg = 0.D0
        Mavg = 0.D0
        meshNavg = 0.D0
        meshMavg = 0.D0

        ! wabs_list loop - parallelize
        do w = 1, nWabs
            Wabs = wabs_list(w)

            do j = 1, nAngles
                theta0 = Angles(j)

                call compute_over_velocities(nVels, vel_list, dlina, rho, rhoN,&
                    &psi, psiN, rhoM, psiM, Dm, theta0, Vopt, Wsc, Wabs, &
                    &numSteps, P1, P2, P3, P4, PN, PM, N_initial, PNvels,PMvels)

                ! Compute the average and variance in the PNvels and PMvels
                Navg = sum(PNvels) / nVels
                Mavg = sum(PMvels) / nVels
                !Nvar = compute_variance(PNvels)
                !Mvar = compute_variance(PMvels)

                ! Assign to the corresponding mesh arrays
                meshNavg(w, j) = Navg
                meshMavg(w, j) = Mavg
                !meshNvar(w, j) = Nvar
                !meshMvar(w, j) = Mvar
            end do
            !f52 = "(2703ES17.8E3)"
            f52 = "(3E17.8E3)"
            ! Critical section
            !!$OMP CRITICAL !write_mesh
            !$ id = 10*OMP_GET_THREAD_NUM()
!            print *, wabsFileList(2*(l - 1) + 1)
            call writeMesh(wabsFileList(2*(l - 1) + 1), meshNavg, w, f52, id)
            call writeMesh(wabsFileList(2*(l - 1) + 2), meshMavg, w, f52, id)
            !!$OMP END CRITICAL !write_mesh

        end do
    end do
    !$OMP END PARALLEL DO
    f50 = "(ES17.8E3)"
    f51 = "(A1, A15)"

    file1 = "masses.txt"
    string1="Mass Delta"
    call writeParameterUsedFile(directoryName, file1, string1, &
            &Masses, nMasses, f51, f50)
    file1 = "angles.txt"
    string1="Angle"
    call writeParameterUsedFile(directoryName, file1, string1, &
            &Angles, nAngles, f51, f50)
    file1 = "wabs.txt"
    string1="W_abs"
    call writeParameterUsedFile(directoryName, file1, string1, &
            &wabs_list, nWabs, f51, f50)

!50  format(ES17.8E3)
!51  format(A1, A15)
!52  format(1608ES17.8E3)

contains
    subroutine writeMesh(file_to_write_to, data_to_write, iter, format, id)
        implicit none
        character(256)          :: file_to_write_to
        real(8), dimension(:, :):: data_to_write
        integer                 :: iter, id, ii
        character(20)           :: format

        logical :: is_open

        inquire(unit = id, opened = is_open)
        if(is_open) then
            !print *, "it's open - track it down"
            close(unit = id)
            !stop
        end if

        open(unit = id, file = file_to_write_to, position = 'append', status = 'old')
        write(unit = id, fmt = trim(adjustl(format))) data_to_write(iter, :)
        close(unit = id)
    end subroutine writeMesh

    subroutine writeParameterUsedFile(directory, fName, header, data, &
        &limit, formatChar, formatFloat)
        implicit none
        real(8), dimension(:)       :: data
        character(256)              :: directory, fName, header
        integer                     :: limit, ii
        character(20)               :: formatChar, formatFloat

        open(unit = 1, file = trim(adjustl(directory)) // trim(adjustl(fName)),&
                &status = "unknown")
        write(unit = 1, fmt = trim(adjustl(formatChar))) "#", trim(adjustl(header))
        do ii = 1, limit
            write(unit = 1, fmt = trim(adjustl(formatFloat))) data(ii)
        end do
        close(unit = 1)
    end subroutine writeParameterUsedFile

    real(8) function computeP_i(P_A, rhojj) result(P_i)
        implicit none
        real(8), intent(in) :: P_A, rhojj
        P_i = P_a * rhojj
        return
    end function computeP_i

    real(8) function computeP_NM(P_a, P_b, rhojj_N, rhojj_M) result(P_NM)
        implicit none
        real(8), intent(in) :: P_a, P_b, rhojj_N, rhojj_M
        P_NM = P_a * rhojj_N + P_b * rhojj_M
        return
    end function computeP_NM

    subroutine compute_probabilities(rrho, rrhoN, rrhoM, ppsi, ppsiN, &
        &ppsiM, DDm, ttheta0, vvel, ttStep, VVopt, WWsc, WWabs, nnumSteps, PP1, PP2, &
        &PP3, PP4, PPN, PPM, reset)
        implicit none
        real(8), dimension(2, 2), intent(in)    :: rrhoN, rrhoM
        real(8), dimension(2, 2)                :: rrho
        complex(8), dimension(2), intent(in)    :: ppsiN, ppsiM
        complex(8), dimension(2)                :: ppsi
        real(8), intent(in)     :: DDm, ttheta0, vvel, ttStep, VVopt, WWsc, WWabs
        real(8), dimension(:)   :: PP1, PP2, PP3, PP4, PPN, PPM
        integer, intent(in)     :: nnumSteps

        logical, intent(in)     :: reset
        real(8), dimension(2, 2):: resetRho
        complex(8), dimension(2):: resetPsi

        real(8)                 :: AA, BB, llambda
        integer                 :: ii

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

        AA = 0.D0
        BB = 0.D0
        llambda = 0.D0

        ! Obtain first element of PN and PM
        PPN(1) = rrho(1, 1)
        PPM(1) = rrho(2, 2)

        do 101 ii = 2, nnumSteps
            ! Reset rho and psi for use
            rrho = resetRho
            ppsi = resetPsi

            ! Iterate over the probability arrays
            call exactBanfor(DDm, vvel, ttheta0, VVopt, WWsc, WWabs, llambda, &
                &AA, BB, ttStep, ppsi, rrho)

            PP1(ii) = PPN(ii - 1) * rrho(1, 1)
            PP2(ii) = PPN(ii - 1) * rrho(2, 2)
            PP3(ii) = PPM(ii - 1) * rrho(1, 1)
            PP4(ii) = PPM(ii - 1) * rrho(2, 2)

            PPN(ii) = PP1(ii) + PP3(ii)
            PPM(ii) = PP2(ii) + PP4(ii)
101     end do
    end subroutine compute_probabilities

    subroutine compute_over_velocities(nVells, vell_list, dllina, rrho, rrhoN, &
    &ppsi, ppsiN, rrhoM, ppsiM, DDm, ttheta0, VVopt, WWsc, WWabs, &
    &nummSteps, PP1, PP2, PP3, PP4, PPN, PPM, N_iinitial, PPNvels, PPMvels)
        implicit none

        complex(8), dimension(2):: ppsi
        complex(8), dimension(2), intent(in)    :: ppsiN, ppsiM
        real(8), dimension(2, 2):: rrho
        real(8), dimension(2, 2), intent(in)    :: rrhoN, rrhoM
        real(8), dimension(:)   :: vell_list, PPN, PPM, PP1, PP2, PP3
        real(8), dimension(:)   :: PP4, PPNvels, PPMvels
        integer:: nVells, nummSteps, kk
        real(8):: dllina, vell, ttStep, llambda, AA, BB
        real(8), intent(in) ::  DDm, ttheta0, VVopt, WWsc, WWabs
        logical, intent(in) ::N_iinitial

        llambda = 0.D0
        AA = 0.D0
        BB = 0.D0

        do kk = 1, nVells
            vell = vell_list(kk)

            ttStep= dllina / vell

            rrho = rrhoN
            ppsi = ppsiN

            call exactBanfor(DDm, vell, ttheta0, VVopt, WWsc, WWabs, llambda, &
                &AA, BB, ttStep, ppsi, rrho)

            PPN(1) = rrho(1, 1)
            PPM(1) = rrho(2, 2)

            call compute_probabilities(rrho, rrhoN, rrhoM, ppsi, ppsiN, ppsiM, &
                &DDm, ttheta0, vell, ttStep, VVopt, WWsc, WWabs, nummSteps, PP1, PP2, &
                &PP3, PP4, PPN, PPM, N_iinitial)

            PPNvels(kk) = PPN(nummSteps)
            PPMvels(kk) = PPM(nummSteps)
        end do
    end subroutine compute_over_velocities

    real(8) function compute_variance(array) result(var)
        implicit none
        real(8), dimension(:), intent(in)   :: array
        real(8) :: ave
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
end program wabs_step