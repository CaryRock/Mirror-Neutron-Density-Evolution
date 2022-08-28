! This is the actual implementation of f90test.f95 that will get the commandline
! options for the program - either from the file "parameters.list", or from the
! command line itself.


! Much of this can probably be ripped from f90test.f95
module get_parameters
contains
    subroutine get_params(INFILE_1, mass, angle)
        use f90getopt
        implicit none
        
        character(128), intent(in)  :: INFILE_1
        character(256)              :: loc
        real(8)                     :: mass, angle
        logical                     :: i_opt, m_opt, t_opt
        
        ! Will 
        type(option_s)  :: opts(5)
        opts(1) = option_s("help",      .false.,    "h")
        opts(2) = option_s("version",   .false.,    "v")
        opts(3) = option_s("mass",      .true. ,    "m")
        opts(4) = option_s("theta",     .true. ,    "t")
        opts(5) = option_s("input_file",.true. ,    "i")

        if(command_argument_count() .eq. 0) then
            print *, "WARNING: No options entered. Assuming '-i'."
        !    opts = 
        end if
        
        i_opt = .false.
        m_opt = .false.
        t_opt = .false.
        
        do
            select case(getopt("i:m:t:hv", opts))
                case(char(0))
                    exit
                case("i")   ! option -i --input-file
                    i_opt = .true. ! i excludes m and t
                    !print  *, opts
                    if(isnum(trim(optarg)) == 0) then
                        loc = trim((optarg)) // "/" // trim((INFILE_1))
                        open(unit = 2, file = trim((loc)))
                        print *, "Reading parameter.list from location " &
                            &// trim(optarg) // "/" // trim(INFILE_1)
                        read(2, *)  ! Mulch the header
                        read(2, *) mass, angle
                    end if
                
                case("m")   ! option -m --mass
                    m_opt = .true. ! m mandataory if not using parameter.list
                    if(isnum(trim(optarg)) > 0) then ! Check for a number in "optarg"
                        !print *, "In mass: " // optarg
                        read(optarg, *) mass
                    else
                        print *, trim(optarg)
                        print *, "ERROR: Option -m or --mass: not a number"
                        stop
                    end if
                case("t")   ! option -t --theta
                    t_opt = .true.  ! t required if using m/not using parameter.list
                    if(isnum(trim(optarg)) > 0) then    ! Check for a number in "optarg"
                        !print *, "In theta: " // optarg
                        read(optarg, *) angle
                    else
                        print *, trim(optarg)
                        print *, "ERROR: Option -t or --theta: not a number"
                        stop
                    end if
                case("h")   ! help output
                    write(*, '(6(A/),/,4(A/))')&
                        "Usage: test.exe [options] ...",&
                        "Options:",&
                        "   -i  --input-file    Use the file 'parameter.list' for &
                            &input parameters.",&
                        "   -m  --mass delta    Mass difference - between 10 neV &
                            &and 10 keV (in eV)",&
                        "   -t  --theta angle   Desired mixing angle (in radians). &
                            &Angles > pi/4 will be set to pi/4",&
                        "   -h  --help          Print this help screen",&
                        "   -v  --version       Print version information"
                case("v")   ! version information
                    print *, "Version 0.1"
            end select
        end do
    
        ! Process input options to check for conflicts - -m -t and -i, etc.
        if (m_opt .and. t_opt .and. .not. i_opt) then
        else if (i_opt .and. .not. m_opt .and. .not. t_opt) then
        else if (m_opt .neqv. t_opt) then
            print *, "Using -m/--mass requires using -t/--theta and not using &
                &-i/--input-file"
    !    else if (t_opt .and. .not. m_opt) then
    !        print *, "Using -t/--theta requires using -m/--mass and not using &
    !            &-i/--input-file"
        else if (i_opt .and. m_opt .or. i_opt .and. i_opt) then
            print *, "Using -i/--input-file disallows using -m/--mass and &
                &-t/--theta."
            stop
        else
            print *, "Unhandled error case. Please run with '-h' to see the &
                &available options for this program."
            stop
        end if
    
    end subroutine get_params
end module get_parameters