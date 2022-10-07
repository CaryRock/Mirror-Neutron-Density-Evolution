program generate
    implicit none
    ! Short, sweet, simple program to generate the deltaM and theta values
    ! for this NNprime code

    ! Just gonna do this dumb, can be extended later.

    ! This program will generate data from dM = 0.5D-9 to 10.0D3 eV, and theta
    ! values from 0.5D-7 to 0.785D0 radians, with 20 steps per decade.

    !real, parameter     :: mMMax = 9.5, mMMin = 1.0, tMMax = 9.5, tMMin = 1.0
    integer, parameter  :: mEMin = -9, mEMax = 3, tEMax = 0, tEMin = -8
    real, parameter  :: mDivsPerDecade = 20, tDivsPerDecade = 20
    real             :: i, j, pow

50  format(ES16.8E3)

    open(unit = 1, file = "deltaMs", status = "unknown")
    j = mEMin

    do j = mEMin, mEMax
        pow = 10.0**j
        i = 1. / mDivsPerDecade
        do while (i .lt. 10.0)
            write(unit = 1, fmt = 50) i * pow
            i = i + 0.05
        end do
    end do
    close(unit = 1)

    open(unit = 1, file = "thetas", status = "unknown")
    j = tEMin
    do j = tEMin, tEMax
        pow = 10.0**j
        i = 1. / tDivsPerDecade
        do while (i .lt. 10.0)
            if (j .eq. tEMax .and. floor(i*100) .eq. 80) then
                i = 0.785
                write(unit = 1, fmt = 50) i * pow
                goto 001
            else 
                write(unit = 1, fmt = 50) i * pow
                i = i + 0.05
            end if
        end do
    end do

001 close(unit = 1)
    stop
end program generate
