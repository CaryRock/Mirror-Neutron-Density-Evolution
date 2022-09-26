! This is the actual implementation of f90test.f95 that will get the commandline
! options for the program - either from the file "parameters.list", or from the
! command line itself.


! Much of this can probably be ripped from f90test.f95
module get_parameters
contains
! TODO: Implement an option to pass specific material files by command line
    subroutine get_params(INFILE_1, INFILE_2, inventory, mass, &
                            &angle, psi, vel, fName, nLines, skip_mid, skip_Wsc)
        use f90getopt
        use material_list
        implicit none
        
        character(32), parameter    :: version = "Version 0.1.4"
        character(128)              :: INFILE_1, INFILE_2
        character(64), intent(in)   :: fName
        integer                     :: nLines
        type(materiallist), dimension(:), allocatable   :: inventory
        character(256)              :: loc
        real(8)                     :: mass, angle, vel
        complex(8), dimension(2)    :: psi
        logical                     :: i_opt, L_opt, m_opt, t_opt, N_opt
        logical                     :: MN_opt, V_opt
        logical                     :: skip_mid, skip_Wsc

        ! Will 
        type(option_s)  :: opts(11)
        opts(1) = option_s("help",      .false.,    "h")
        opts(2) = option_s("version",   .false.,    "v")
        opts(3) = option_s("mass",      .true. ,    "m")
        opts(4) = option_s("theta",     .true. ,    "t")
        opts(5) = option_s("input_file",.true. ,    "i")
        opts(6) = option_s("initial_n", .true. ,    "N")
        opts(7) = option_s("initial_np",.true. ,    "M")
        opts(8) = option_s("material",  .true. ,    "L")
        opts(9) = option_s("velocity",  .true. ,    "V") 
        opts(10)= option_s("skip-mid",  .false.,    "S")
        opts(11)= option_s("skip_Wsc",  .false.,    "c")

        if(command_argument_count() .eq. 0) then
            print *, "WARNING: No options entered. Please run with '-h' for &
                &the help and to see what options are available."
            STOP
        end if
        
        i_opt   = .false.
        m_opt   = .false.
        t_opt   = .false.
        N_opt   = .false.
        MN_opt  = .false.
        V_opt   = .true.
        skip_mid= .false.
        skip_Wsc= .false.
        L_opt   = .false.

        do
!            select case(getopt("ci:NL:Mm:St:V:hv", opts))
            select case(getopt("ci:NL:Mm:St:hv", opts))
                case(char(0))
                    exit

                case("c")   ! option -c --skip_Wsc
                    skip_Wsc = .true.   ! DO exclude W_scattering

                case("i")   ! option -i --input-file
                    i_opt = .true. ! i excludes m and t
                    !print  *, opts
                    if(isnum(trim(optarg)) == 0) then
                        loc = trim((optarg)) // "/" // trim((INFILE_1))
                        open(unit = 2, file = trim((loc)))
                        !print *, "Reading parameter list from file " &
                        !    &// trim(optarg) // "/" // trim(INFILE_1)
                        read(2, *)  ! Mulch the header
                        read(2, *) mass, angle, vel
                    end if

                case("m")   ! option -m --mass
                    m_opt = .true. ! m mandataory if not using parameter.list
                    if(isnum(trim(optarg)) > 0) then ! Check for a number in "optarg"
                        !print *, "In mass: " // optarg
                        read(optarg, *) mass
                    else
                        print *, trim(optarg)
                        print *, "ERROR: Option -m or --mass: not a number"
                        STOP
                    end if

                case("t")   ! option -t --theta
                    t_opt = .true.  ! t required if using m/not using parameter.list
                    if(isnum(trim(optarg)) > 0) then    ! Check for a number in "optarg"
                        !print *, "In theta: " // optarg
                        read(optarg, *) angle
                    else
                        print *, trim(optarg)
                        print *, "ERROR: Option -t or --theta: not a number"
                        STOP
                    end if

                case("L")   ! option -L --material
                    L_opt = .true.
                    INFILE_2 = trim(optarg)
                    !write(6, *) "Will read inventory from " // INFILE_2

                case("N")   ! option -N --initial_n
                    N_opt   = .true.  ! Conflicts with option M
                    psi(1)  = cmplx(1.0, 0.0)
                    psi(2)  = cmplx(0.0, 0.0)

                case("M")   ! option -M --initial_np
                    MN_opt  = .true.  ! Conflicts with option N
                    psi(1)  = cmplx(0.0, 0.0)
                    psi(2)  = cmplx(1.0, 0.0)

                case("S")   ! option -S --skip-mid
                    skip_mid= .true.

!                case("V")   ! option -V --velocity
!                    V_opt   = .true.
!                    if(isnum(trim(optarg)) > 0) then    ! Check for a number in "optarg"
!                        read(optarg, *) vel
!                    else
!                        print *, trim(optarg)
!                        print *, "ERROR: Option -V or --velocity: not a number"
!                        STOP
!                    end if

                case("h")   ! help output
                    write(*, '(6(A/),/,4(A/))')&
                        "Usage: velmesh [options] ...",&
                        "Options:",&
                        ! -c
                        "   -i  --input-file    Use the file 'parameter.list' &
                            &for the input parameters of mass (-m), angle (-t),&
                            & and velocity (-V).",&
                        "   -m  --mass          Mass difference between n and &
                            &n'. Input a value between between 10 neV. &
                            &and 10 keV (in eV)",&
                        "   -t  --theta angle   Desired mixing angle (in &
                            &radians). Angles > pi/4 will be set to pi/4",&
!                        "   -V  --velocity      Velocity of initial &
!                            &particles.",&
                        "   -h  --help          Print this help screen",&
                        "   -v  --version       Print version information"

                case("v")   ! version information
                    print *, version
            end select
        end do
    
        ! Process input options to check for conflicts - -m -t and -i, etc.
        if (m_opt .and. t_opt .and. V_opt .and. .not. i_opt) then
        else if (i_opt .and. .not. m_opt .and. .not. t_opt .and. .not. V_opt) then
        else if (m_opt .neqv. t_opt .or. m_opt .neqv. V_opt) then
            print *, "Using -m/--mass requires using -t/--theta and &
                &-V/--velocity and not using -i/--input-file"
    !    else if (t_opt .and. .not. m_opt) then
    !        print *, "Using -t/--theta requires using -m/--mass and not using &
    !            &-i/--input-file"
            STOP
        else if (i_opt .and. m_opt .or. i_opt .and. t_opt .or. i_opt .and. V_opt) then
            print *, "Using -i/--input-file disallows using -m/--mass, &
                &-t/--theta, and -V/--velocity."
            STOP
        else
            print *, "Unhandled error case. Please run with '-h' to see the &
                &available options for this program."
            STOP
        end if
        
        if (N_opt .and. MN_opt) then
            print *, "Please specify the initial starting state as being either&
            &purely neutron ('N') or purely mirror neutron ('M'). It makes no &
            &sense to be both." 
            STOP
        else if (.not. N_opt .and. .not. MN_opt) then
            ! Just default to starting as neutron. Figure out specifying 
            ! starting angle later.
            psi(1)  = cmplx(1.0, 0.0)
            psi(2)  = cmplx(0.0, 0.0)
        end if

        if (.not. V_opt) then
            ! Read velocity from data file
            continue
        end if

        nLines = 0
        nLines = get_lines(INFILE_2)
        allocate(inventory(nLines))
        call get_materials(inventory, fName, nLines, INFILE_2)
        
        if (skip_Wsc) then
            inventory%Wsc = 0.D0
        end if
    end subroutine get_params
end module get_parameters
