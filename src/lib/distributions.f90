module constants
    real, parameter :: kB = 8.617E-5    ! eV / K
    real, parameter :: pi = 4.0 * atan(1.0)
    real, parameter :: c  = 2.99792E8   ! m/s
    real, parameter :: mn = 939.5654133E6 / c**2    ! Neutron mass, MeV / c^2
end module constants

module distributions
contains
    real function MBDist(vel, x) result(pdfx)
        implicit none
        real, intent(in):: vel, x

        pdfx = x * vel*vel*vel * exp(-vel*vel*x/2.0)
    end function MBDist
    
    subroutine setup_MBDist(num_samples, temp, velocities, PDF, CDF, vmin, vmax)
        ! Idea: Evaluate 1e6 points along the MB Distribution, using a normal
        ! distribution as the "starting point"/overlay to weight which values to
        ! start sampling from. Somthing like importance sampling and all that,
        ! but way more dumb.
        use constants
        implicit none
        integer, intent(in) :: num_samples
        real, intent(in)    :: temp, vmin, vmax

        real, parameter     :: r2opi = sqrt(2/pi)
        real                :: x, dv

        real, dimension(:)  :: velocities, PDF, CDF
        integer                         :: i

        x = mn / (kB * temp)
        PDF = 0.0
        CDF = 0.0
        dv = (vmax - vmin) / num_samples

        ! Establish a range of values to evaluate distribution over
        
        do i = 1, num_samples + 1
            velocities(i) = vmin + (i - 1)* dv   ! Set velocity
            PDF(i) = r2opi * MBDist(velocities(i), x) ! Evaluate at velocity
        end do

        ! Accumulate values in the CDF
        CDF(1) = PDF(1)
        do i = 2, num_samples
            CDF(i) = CDF(i - 1) + PDF(i)
        end do
    end subroutine setup_MBDist

    real function sample_MBDist(velocities, CDF, num_samples) result(sampled)
        implicit none
        real, dimension(:), intent(in)  :: CDF, velocities
        real                            :: random

        integer, intent(in)             :: num_samples
        integer(16)                     :: index

        call random_number(random)

        print *, random
        ! TODO: implement a method to translate random number calls to velocities
        ! from the CDF -> PDF inversion(s).

        ! We know how many elements are in the CDF and PDF, so utilize
        index = int(random * num_samples) + 1
        index = int(cdf(index) * num_samples) + 1

        sampled = velocities(index)
        
    end function sample_MBDist
end module distributions