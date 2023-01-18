module NNp_file_writing
    type mesh_file_data
        complex,     dimension(2)    :: psi
        real,        dimension(2, 2) :: rho
        character(:), dimension(:, :), allocatable   :: rhonn
        character(:), allocatable    :: dist, thetaChar, massChar, xChar
        character(:), allocatable    :: DmChar, theta0char, thickness
        character(64)                   :: filename
        character(36)                   :: uuid
        real                         :: Dm, theta0, x
        logical                         :: scattering, absorption
        
        integer                         :: prec
        
        ! Formatting
        character(20)   :: f56! = "A16,5A16"
        character(20)   :: f57! = "6ES16.6E3"
    end type mesh_file_data

contains
    type(mesh_file_data) function create_msd_struct(psi, filename, thickness, &
        &DmChar, theta0char, uuid, massChar, thetaChar, Dm, theta0, xChar, &
        &dist, rhonn, x, rho, Wsc, Wabs, prec) result(mfp)
        implicit none
        !integer, parameter              :: prec = 18
        integer, intent(in)             :: prec
        complex,     dimension(2)    :: psi
        real,        dimension(2, 2) :: rho
        character(prec),dimension(2, 2) :: rhonn
        character(prec)                 :: dist, thetaChar, massChar, xChar
        character(prec)                 :: DmChar, theta0char, thickness
        character(64)                   :: filename
        character(36)                   :: uuid
        real                         :: Dm, theta0, x
        logical                         :: Wsc, Wabs

        character(2)                    :: ff56, ff57

        mfp%psi             = psi
        mfp%rho             = rho
        mfp%rhonn           = rhonn
        mfp%dist            = dist
        mfp%thetaChar       = thetaChar
        mfp%massChar        = massChar
        mfp%xChar           = xChar
        mfp%DmChar          = DmChar
        mfp%theta0char      = theta0char
        mfp%thickness       = thickness
        mfp%filename        = filename
        mfp%uuid            = uuid
        mfp%Dm              = Dm
        mfp%theta0          = theta0
        mfp%x               = x
        mfp%scattering      = Wsc
        mfp%absorption      = Wabs

        ! Get the formatting to update with the rest of the program
        write(ff56, "(A2)") prec
        write(ff57, "(A2)") prec
        mfp%f56             = "A" // ff56 // ",5A6"
        mfp%f57             = "6ES" // ff57 // ".6E3"
        return
    end function create_msd_struct

    character(256) function set_directory(mfp, prog_type) result(directoryName)
        implicit none
        type(mesh_file_data)    :: mfp
        character(4)            :: prog_type
        
        directoryName = "./data/" // prog_type // "/"

        ! if psi = (1, 0)
        !   directoryName = directoryName // "N/"
        ! else if psi = (0, 1)
        !   directoryName = directoryName // "Mn/"
        if (mfp%psi(1) == cmplx(1.0, 0.0) .and. mfp%psi(2) == &
                &cmplx(0.0, 0.0)) then
            directoryName = trim(adjustl(directoryName)) // "N/"
        else if (mfp%psi(1) == cmplx(0.0, 0.0) .and. mfp%psi(2) == &
                    &cmplx(1.0, 0.0)) then
            directoryName = trim(adjustl(directoryName)) // "Mn/"
        else
            directoryName = trim(adjustl(directoryName)) // "mix/"
        end if
        
        ! if scattering ( == true)
        !   directoryName = directoryName // "with_sc/"
        ! else
        !   directoryName = directoryName // "no_sc/"
        if (mfp%scattering) then
            directoryName = trim(adjustl(directoryName)) // "no_sc/"
        else
            directoryName = trim(adjustl(directoryName)) // "with_sc/"
        end if

        directoryName = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) // "/"
    end function set_directory

    subroutine mesh_file_prepare_ave(mfp, directoryName, file, element, u)
        implicit none
        type(mesh_file_data)    :: mfp
        character(256)          :: directoryName, file
        character(3)            :: element
        character               :: unt
        integer                 :: u

        write(unt, '(I1)') u

        directoryName = set_directory(mfp, "mesh")
        
        call system('mkdir -p ' // trim(adjustl(directoryName)))

        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_Mesh_" // unt // "_ave_data_" // mfp%uuid // ".dat"
        
        open(unit = u*10, file = trim(adjustl(file)), status = "unknown")
        write(u*10, "(A16)") "#avg(" // trim(element) // ")"
        
        return
    end subroutine mesh_file_prepare_ave

    subroutine step_file_prepare_ave(mfp, directoryName, file, element, u)
        implicit none
        type(mesh_file_data)    :: mfp
        character(256)          :: directoryName, file
        character(3)            :: element
        character               :: unt
        integer                 :: u

        write(unt, '(I1)') u

        directoryName = set_directory(mfp, "step")

        call system('mkdir -p ' // trim(adjustl(directoryName)))

        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_Step_" // unt // "_ave_data_" // mfp%uuid // ".dat"

        open(unit = u*10, file = trim(adjustl(file)), status = "unknown")
        write(u*10, "(a16)") "#avg(" // trim(element) // ")"

        return
    end subroutine step_file_prepare_ave

    subroutine wabs_file_prepare_ave(mfp, directoryName, file, element, u, m)
        implicit none
        type(mesh_file_data)    :: mfp
        character(256)          :: directoryName, file
        character(4)            :: element, mnt
        character               :: unt
        integer                 :: u, m
        
        write(unt, '(I1)') u
        write(mnt, '(I4.4)') m

        directoryName = set_directory(mfp, "wabs")
        call system('mkdir -p ' // trim(adjustl(directoryName)))

        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_Wabs_" // unt //"_" // mnt // "_ave_data_" // mfp%uuid // ".dat"
        
        open(unit = u, file = trim(adjustl(file)), status = "unknown")
        write(u, "(A16)") "#ave(" // trim(adjustl(element)) // ")"
        close(unit = u)
    end subroutine wabs_file_prepare_ave

    subroutine mesh_file_prepare_var(mfp, directoryName, file, element, u)
        implicit none
        type(mesh_file_data)    :: mfp
        character(256)          :: directoryName, file
        character(3)            :: element
        character               :: unt
        integer                 :: u

        write(unt, '(I1)') u

        directoryName = set_directory(mfp, "mesh")

        call system('mkdir -p ' // trim(adjustl(directoryName)))

        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_Mesh_" // unt // "_var_data_" // mfp%uuid // ".dat"

        open(unit = u, file = trim(adjustl(file)), status = "unknown")
        write(u, "(A16)") "#var(" // trim(adjustl(element)) // ")"
        
        return
    end subroutine mesh_file_prepare_var

    subroutine step_file_prepare_var(mfp, directoryName, file, element, u)
        implicit none
        type(mesh_file_data)    :: mfp
        character(256)          :: directoryName, file
        character(3)            :: element
        character               :: unt
        integer                 :: u

        write(unt, '(I1)') u

        directoryName = set_directory(mfp, "step")

        call system('mkdir -p ' // trim(adjustl(directoryName)))

        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_Step_" // unt // "_var_data_" // mfp%uuid // ".dat"

        open(unit = u*10, file = trim(adjustl(file)), status = "unknown")
        write(u*10, "(A16)") "#var(" // trim(adjustl(element)) // ")"

        return
    end subroutine step_file_prepare_var

    subroutine file_prepare(mfp, prog_type, dat_type, &
                                        &directoryName, file, element, u)
        implicit none
        type(mesh_file_data)    :: mfp
        character(4)            :: prog_type
        character(3)            :: dat_type
        character(256)          :: directoryName, file
        character(3)            :: element
        character               :: unt
        integer                 :: u

        write(unt, '(I1)') u

        directoryName = set_directory(mfp, prog_type)

        call system('mkdir -p ' // trim(adjustl(directoryName)))
        
        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_" // trim(adjustl(prog_type)) // "_" // unt // "_" &
            &// trim(adjustl(dat_type)) // "_data_" // mfp%uuid // ".dat"

        !print *, file

        open(unit = u*10, file = trim(adjustl(file)), status = "unknown")
        write(u*10, "(A16)") "#"// dat_type // "(" // &
                                            &trim(adjustl(element)) // ")"

        return
    end subroutine file_prepare

    subroutine mesh_file_write(rho_av, rho_vr)
        implicit none
        real, dimension(4), intent(in)   :: rho_av, rho_vr

        write(1, "(8ES16.6E3)") rho_av(1), rho_av(2), rho_av(3), rho_av(4), &
            &rho_vr(1), rho_vr(2), rho_vr(3), rho_vr(4)
    
    end subroutine mesh_file_write

    subroutine coord_file_prepare()
        implicit none
        type(mesh_file_data)    :: mfp    ! "mesh_file_prep"
        character(256)          :: directoryName, file

        directoryName = set_directory(mfp, "coor")
        call system('mkdir -p ' // trim (adjustl(directoryName)))
        
        file = trim(adjustl(directoryName)) // trim(adjustl(mfp%filename)) &
            &// "_Dm-" // trim(adjustl(mfp%DmChar)) // "_theta0-" // &
            &trim(adjustl(mfp%theta0char)) // "_data_" // mfp%uuid &
            &// ".dat"

        open(unit = 1, file = trim(adjustl(directoryName)) // &
                &trim(mfp%filename) // "_Dm-" // trim(adjustl(mfp%DmChar)) &
                &// "_theta0-" // trim(adjustl(mfp%theta0Char)) // &
                &"_data_" // mfp%uuid // ".dat", status = "unknown")
        write(1, *) '#' // mfp%massChar, mfp%thetaChar
        write(1, '(2F16.11)') mfp%Dm, mfp%theta0
        write(1, mfp%f56) '#' // mfp%xChar, mfp%dist, mfp%rhonn(1, 1), &
            &mfp%rhonn(2, 1), mfp%rhonn(1, 2), mfp%rhonn(2, 2)
        write(1, mfp%f57) mfp%x, 0.0, &
            &mfp%rho(1, 1), mfp%rho(2, 1), mfp%rho(1, 2), mfp%rho(2, 2)

    end subroutine coord_file_prepare

end module NNp_file_writing