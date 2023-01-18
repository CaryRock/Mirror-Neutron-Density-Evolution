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

    recursive subroutine YK_MAXW(TSOUR, V)
      implicit none
      real, parameter :: AK2M = 1.648656E4
      real, parameter :: AK3M = 2.472984E4

      real, intent(in):: TSOUR
      real            :: V, VPFT, TEMPE, V3MP2, V2MP2, A, AF, AMAX, VMP
      real            :: random

      VPFT = 6637.43        ! For 342 K thermal source

      if (TSOUR .lt. 273.15) VPFT = 1.E8  ! For colder source
!     VPFT - m/s - two regions of thermal flux parameterization

!     If the temperature of the source is constant:
      TEMPE = TSOUR
!     Temperature is a radial function:
!     (for cold source) NOT INITIATED
!     REMIS = sqrt(X0*X0 + Y0*Y0)
!     TEMPE = TVSR(REMIS)   ! TVSR isn't in the CERNLIB documentation.

!     V3MP2 - most probably velocity, squared - (m/s)^2
      V3MP2 = AK3M * TEMPE
      V2MP2 = AK2M * TEMPE
      VMP = sqrt(V3MP2)
      AMAX = 45000. * exp(-1.5)/VMP

      ! Start sampling the distribution
      AF = 1.0
      A = 10.0
      do while(A .gt. AF)
        call random_number(random)
        V = 10. * VMP * random
        call random_number(random)
        A = AMAX * random
        if (V .lt. VPFT) then
          AF = 20000. * exp(-V * V / V2MP2) * V * V * V / (V2MP2 * V2MP2)
        else
          AF = 493.597 / V
        end if
      end do

      return
    end subroutine YK_MAXW
end module distributions