program f90test
    use f90getopt
    implicit none

    real(8)         :: Dm, theta0
    logical         :: i_present, m_present, t_present

    type(option_s)  :: opts(5)
    opts(1) = option_s("help",      .false.,    "h")
    opts(2) = option_s("version",   .false.,    "v")
    opts(3) = option_s("mass",      .true. ,    "m")
    opts(4) = option_s("theta",     .true. ,    "t")
    opts(5) = option_s("input_file",.true. ,    "i")

    if(command_argument_count() .eq. 0) then
        print *, "ERROR: This program has required options. Please check them &
            &out with '-h'"
        stop
    end if
    
    i_present = .false.
    m_present = .false.
    t_present = .false.
    
    do
        select case(getopt("i:m:t:hv", opts))
            case(char(0))
                exit
            case("i")   ! option -i --input-file
                i_present = .true. ! i excludes m and t
                if(isnum(trim(optarg)) == 0) then
                    print *, "Reading parameter.list..."
                    print *, "In input: " // optarg
                end if
            
            case("m")   ! option -m --mass
                m_present = .true. ! m mandataory if not using parameter.list
                if(isnum(trim(optarg)) > 0) then ! Check for a number in "optarg"
                    print *, "In mass: " // optarg
                    read(optarg, *) Dm
                else
                    print *, "ERROR: Option -m or --mass: not a number"
                    stop
                end if
            case("t")   ! option -t --theta
                t_present = .true.  ! t required if using m/not using parameter.list
                if(isnum(trim(optarg)) > 0) then    ! Check for a number in "optarg"
                    print *, "In theta: " // optarg
                    read(optarg, *) theta0
                else
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
    if (m_present .and. t_present .and. .not. i_present) then
        print *, "Cool and good"
    else if (i_present .and. .not. m_present .and. .not. t_present) then
        print *, "Also cood and good"
    else if (m_present .neqv. t_present) then
        print *, "Using -m/--mass requires using -t/--theta and not using &
            &-i/--input-file"
!    else if (t_present .and. .not. m_present) then
!        print *, "Using -t/--theta requires using -m/--mass and not using &
!            &-i/--input-file"
    else if (i_present .and. m_present .or. i_present .and. i_present) then
        print *, "Using -i/--input-file disallows using -m/--mass and &
            &-t/--theta."
    else
        print *, "Unhandled error case."
        stop
    end if
end program f90test