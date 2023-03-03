program test_MW
  use distributions
  implicit none
  real  :: temp, vel, minVel, maxVel
  integer, parameter  :: nVels = 10000
!  real  :: velocities(nVels)
  integer :: i

  temp = 342.
  minVel = 200.
  maxVel = 7153.416954

  do i = 1, nVels
    call neutron_MW_dist(temp, vel, minVel, maxVel)
    write(1, '(1ES17.8E3)') vel
  end do

end program test_MW
