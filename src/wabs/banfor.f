      PROGRAM banfor
c
c-----April 04, 2022 
c-----Barezhiani ANalitical FORmulas 
c-----from paper "N-N' in absorbing material" 
c-----Symmetry 2022, 1, 0. 
c-----https://doi.org/10.3390/sym1010000
c
c --> U = V +/- Dm  +/- muB   will treat that as single real number U=V
c --> W=W1+v[m/s]*W2 in eV
c --> Dm=Delta_m =mn'-mn = e.g. 2.E-7 eV  
c --> eps = epsilon e.g. 2.E-9 eV e.g. in vacuum for theta_0=1.E-2
c 
c --> units: m, sec, eV
c
      Real hbar,nmass,Qe,Dm,vel,s2ze,TOF,theta0
      Real theta,omega,Vopt,V,W,W1,W2,A,B,U1,U2,WW1,WW2,DE,dlina
      COMPLEX(8) pnn8,pmn8
      REAL(8) pnn,pmn,unitar, ll
      COMPLEX i,zeta,D,zeta2,arg,cze,cze2,sze,sze2,H1,H2,arg2
      COMPLEX czestar,szestar,cze2star,sze2star,uroo,roo
      COMPLEX(8) snn,smn,H1C,H2C,snnC,smnC
      COMPLEX(8) aee1,aee2,ee1,ee2,ee1C,ee2C,ee1CC,ee2CC
      INTEGER nti,ntf,ntdiff
c
c --> Added by Cary
      real temp
c
c --> constants and parameters:
      do 500 jo=1,9
      do 501 jn = 1, 20
c      theta0=jn*10.**(-7+jo)
      Dm = 165.9E-9!log10(10**(float(jn)/20.)) * 10**(-15. + float(jo))
      theta0 = 0.001!0.70795
      if(theta0.ge.0.8)theta0=0.785398163
c      write(6,*)theta0
c     
c      Dm=100.0E-9           ! delta_m in 
      hbar=6.582119569E-16 ! precise eV*s hbar
      Vopt=1.659E-7
      V=Vopt - Dm          ! V=200 neV corresponding to B4C
c      W1=6.102E-9          ! constant part of W for B4C
c      W2=2.397E-14         ! part of W to be multiplied by velocity in m/s B4C
      W1=10**(-15. + float(jo * 20 + jn) / 20.)!log10(10**(float(jn)/20.)) * 10**(-15. + float(jo))!2.959E-15!8.4558E-9         ! constant part of W for Cd
      W2=0.!2.137E-16!9.914E-15         ! part of W to be multiplied by velocity in m/s Cd
      eps=0.5*ABS(Dm)*tan(2.*theta0)! eps for givem Dm and theta0 
      dlina = 0.022
c      eps=0.0            
c
      i=CMPLX(0.0,1.0)     !  sqrt(-1.) 
c
      vel=1000.            !  m/s
      !call MAXWV(temp, vel)
c --> start printing list of parameters and tests: 
c      write(6,*)'absorber Cd 3.5 mm'
c      write(6,*)'absorber Cd 3.5 mm'  
c      write(6,*)'neutron velocity          = ',vel
c      write(6,*)'theta_0                   = ',theta0
c      write(6,*)'Delta m                   = ',Dm
c      write(6,*)'epsilon                   = ',eps
c      write(6,*)'Vopt (Re) Fermi potential = ',Vopt
c
      W=W1+vel*W2
c
c      write(6,*)'W (Im) fermi potential    = ',W
c      write(6,*)'Eff. V = Vopt - Dm        = ',V
c
      arg=2.*eps/(V-i*W)
      zeta2=atan(arg) 
      zeta=zeta2/2.
c
c      write(6,*)'tan 2*zeta value in (7)   = ',arg
c      write(6,*)'2*zeta in (7)             = ',zeta2
c      write(6,*)'zeta from (7)             = ',zeta
c      write(6,*) 
c
      cze=cos(zeta)
      sze=sin(zeta)
      cze2=cze**2
      sze2=sze**2
c
c      write(6,*)'cos(zeta) - complex       = ',cze
c      write(6,*)'sin(zeta) - complex       = ',sze
c      write(6,*)'cos^2(zeta) - complex     = ',cze2
c      write(6,*)'sin^2(zeta) - complex     = ',sze2
c
      czestar=CONJG(cze)
      szestar=CONJG(sze)
      cze2star=CONJG(cze2)
      sze2star=CONJG(sze2)
c      
c      write(6,*)'cos(zeta) - conjugated    = ',czestar
c      write(6,*)'sin(zeta) - conjugated    = ',szestar
c      write(6,*)'cos^2(zeta) - conjugated  = ',cze2star
c      write(6,*)'sin^2(zeta) - conjugated  = ',sze2star
c      write(6,*)
c
      uroo=(V-i*W)**2 + 4.*eps**2
      roo =sqrt(uroo)
      if(REAL(zeta).lt.0.0)roo=-roo
      H1=0.5*(V-i*W+roo)
      H2=0.5*(V-i*W-roo)
      H1C=CONJG(H1)
      H2C=CONJG(H2)
c
c      write(6,*)'(V-i*W)^2 + 4.*eps^2      = ',uroo
c      write(6,*)'sqrt[(V-i*W)^2+4.*eps^2]  = ',roo
c      write(6,*)'H1 eigenvalue complex     = ',H1
c      write(6,*)'H2 eigenvalue complex     = ',H2
c      write(6,*)'H1_conjugated             = ',H1C
c      write(6,*)'H2_conjugated             = ',H2C
c
      U1=REAL(H1)
      U2=REAL(H2)
      DE=U1-U2
      WW1=AIMAG(H1)
      WW2=AIMAG(H2)
c
c      write(6,*)'REAL part of H1           = ',U1
c      write(6,*)'REAL part of H2           = ',U2
c      write(6,*)'dE = Re(H1) - Re(H2)      = ',DE
c      write(6,*)'IMAGINARY part of H1      = ',WW1
c      write(6,*)'IMAGINARY part of H2      = ',WW2
c
      time=dlina/vel       !   Cd absorber
      TOF=time/hbar
c      write(6,*)'time over absorber, sec   = ',time
c      write(6,*)'TOF =  time/hbar eV^-1    = ',TOF   
c 
      aee1=-i*H1*TOF
      ee1=cdexp(aee1) 
      aee2=-i*H2*TOF    
      ee2=cdexp(aee2)
      ee1C=CONJG(ee1)
      ee2C=CONJG(ee2)
c
c      write(6,*)'exp(-i*H1*TOF) complex    = ',ee1
c      write(6,*)'exp(-i*H2*TOF) complex    = ',ee2
c      write(6,*)'exp(-i*H1*TOF) conjugate  = ',ee1C
c      write(6,*)'exp(-i*H2*TOF) conjugate  = ',ee2C
c     write(6,*)
c
      snn =cze2*ee1+sze2*ee2
      snnC=cze2star*ee1C+sze2star*ee2C
      smn =cze*sze*(ee1-ee2)
      smnC=czestar*szestar*(ee1C-ee2C)
      pnn=snn*snnC
      pmn=smn*smnC
      unitar=pnn+pmn
c
c      write(6,*)'snn                       = ',snn
c      write(6,*)'snnC                      = ',snnC
c      write(6,*)'smn                       = ',smn
c      write(6,*)'smnC                      = ',smnC
c      write(6,*)'pnn                       = ',pnn
c      write(6,*)'pmn                       = ',pmn
c      write(6,*)'unitar                    = ',unitar
      ll = (1 - exp(-2. * dlina * W1 / (hbar * vel)))
      write(6,57)W1,pnn,pmn, ll
c  57  format(F18.12,1P2E16.6)
c  57  format(1P2E18.7, 1P2E18.6)
  57  format(4ES18.8E3)
 501  enddo
 500  enddo
c
      STOP
      end

      subroutine MAXWV(temp, vel)
      implicit none
      parameter (AK2M = 1.648656E+4)    ! 2*k_B / m_n
      parameter (AK3M = 2.472984E+4)    ! 3*k_B / m_n

      real VPFT, TEMPE, V2MP2, V3MP2, VMP, AMAX, A, AF

      VPFT = 6637.43                    ! For 342K thermal source
      if(temp .lt. 273.)VPFT = 1.E+8    ! For cold(er) source
c --> VPFT - [m/sec] two regions of thermal flux parameterization
c --> if the temperature of the source is constant:
      TEMPE = temp
c --> Temperature is a radial function:
c --> (for cold source) NOT INITIATED
c     REMIS = sqrt(X0**2 + Y0**2)
c     TEMPE = TVSR(REMIS)
c
c --> Used FLUX spectrum is divided by factor 10^11;
c     V3MP2 - most probable velocity square (m/s)**2
      V3MP2 = AK3M * TEMPE
      V2MP2 = AK2M * TEMPE
      VMP = sqrt(V3MP2)
      AMAX = 45000.*EXP(-1.5)/VMP
 1    continue
      V = 10.*VMP*RNDM(1.)
      A = AMAX*RNDM(1.)
      if(vel.lt.VPFT)then
      AF = 20000.*exp(-vel**2/V2MP2)*V**3 / V2MP2**2
                     else
      AF = 493.597/vel
                     endif
      if(A.gt.AF)go to 1

      return
      end subroutine MAXWV