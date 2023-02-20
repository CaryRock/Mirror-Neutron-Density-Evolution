      PROGRAM banforb4c
      use exact_banfor_module
c
c-----January 29, 2023
c-----Berezhiani ANalitical FORmulas 
c-----from paper "N-N' in absorbing material" 
c-----Symmetry 2022, 1, 0. 
c-----https://doi.org/10.3390/sym1010000
c-----wide range of theta_0 vs wide range dm
c
c --> U = V +/- Dm +/- muB   will treat that as single real number U=V
c --> W=W1+v[m/s]*W2 in eV
c --> Dm=Delta_m =mn'-mn = e.g. 2.E-7 eV  
c --> eps = epsilon e.g. 2.E-9 eV e.g. in vacuum for theta_0=1.E-2
c 
c --> units: m, sec, eV
c
      
      real rho(2, 2)
      complex psi(2)

      Real theta0,TOF
      Real V,W!,U1,U2,WW1,WW2,DE
      REAL     Dm!,vre12,vim12,vre21,vim21,vre11,vre22
      REAL     xdm(260),yt0(158),array(260,158,4)!,probave(4)
      REAL    specn(300),speca(200)
      REAL     p1,p2,p3,p4,PP1,PP2,PP3,PP4
c     array(Dm,theta0,1=ro11,ro22,Re ro12,Im ro 12)
      COMPLEX i!,zeta,zeta2,arg,cze,cze2,sze,sze2,H1,H2
      !COMPLEX czestar,szestar,cze2star,sze2star,uroo,roo
      !COMPLEX     snn,smn,smm,H1C,H2C,snnC,smnC,smmC
      !COMPLEX     rho11,rho22,rho12,rho21
      !COMPLEX     psi1t,psi2t,psi1tC,psi2tC
      !COMPLEX     aee1,aee2,ee1,ee2,ee1C,ee2C
c
      common /input/Dm,theta0,V,W,del,nelsc,TOF
c      common /output/vre12,vim12,vre21,vim21,vre11,vre22
      common /output/p1,p2,p3,p4
      common /accu/PP1,PP2,PP3,PP4
c
      !COMMON/PAWC/H(50000)
      !CALL HLIMIT(50000)
      !CALL HBOOK1(15,'velocity spectrum m/s',300,0.,12000.,0.)
      !CALL HBOOK1(16,'wavelength in AA',200,0.,20.,0.)
c
c --> initial parameters for 8 mm B4C layer
      i=CMPLX(0.0,1.0)      ! sqrt(-1.) 
c     TSOUR=340.            ! Kelvin
      dlina=0.008           ! m thickness of B4C layer = 8
      hbar=6.582119569E-16  ! precise eV*s hbar
c      Vopt=1.659E-7         ! for heavy water
      Vopt=1.992E-7         ! for B4C
c     Vopt=5.877E-8         ! for Cd
c      W1=2.9588E-15         ! for heavy water constant part (absorption)
c      W1=0.0                ! for test
      W2=0.0                ! correction
c     W2=2.1367E-14         ! for heavy water to be x by velocity
      W1=6.102E-9           ! constant part of W for B4C
c     W2=2.397E-14          ! part of W to be x by velocity in m/s B4C
c     W1=8.4558E-9          ! constant part of W for Cd
c     W2=9.914E-15          ! part of W to be x by velocity in m/s Cd
      NSTAT=1               ! 10000 for approx 3 hours D2O
      nelsc=1               !  number of elastic segments 117 for D2O 
      write(6,613) TSOUR,dlina,NSTAT,nelsc,Vopt,W1,W2
  613 format('Compare B4C 8 mm 1steps 2300 m/s with Cary',/
     *       'thermal spectrum at TSOUR = ', F10.3,/
     *       'thickness in m,           = ', F10.3,/
     *       'NSTAT - averaging vel     = ', I10,/
     *       '# of elastic scat in layer= ', I10,/
     *       'V_optical Real part, eV   = ', 1PE12.4,/
     *       'W1=W_abs Imaginary, eV    = ', E12.4,/
     *       'W2=W_scat Imaginary, eV   = ', E12.4/)
c -------------------------------------------------------------------- 
      Vop=Vopt       !  Re optical potential in eV  
      del=dlina      !  elastic scattering length or layer thickness, m

c --> for del=2.05 cm neslsc=126 (Cary) in D2O
c
c --> begining of the process
c --> loop on Dm 260 points: from 1 neV to 8912.5 eV
c
      Dm=2.E-7
      V=Vop-Dm
      Vel=2300.
      W=W1
      write(6,17)Vel,Dm
  17  format('calculation for velocity = ',F12.3,' and Dm = ',1PE12.3)
c
c --> loop on theta_0 158 points: from 1E-8 rad to <0.7853 rad
      DO 210 jp=1,8
      do 210 jm=1,20
      xt20=0.05*(jm-1)
      theta0=10**(xt20+jp-9.)
      if(theta0 .gt. 0.785398)go to 210
      indy=jm+20*(jp-1)
      yt0(indy)=theta0
c
      TOF=del/(Vel*hbar)
c     
      rho = 0.0
      rho(1, 1) = 1.0
      call exactBanfor(Dm,vel,theta0,Vopt,W2,W1,del/vel,psi,rho)
      p1 = rho(1, 1)
      p2 = rho(1, 2)
      p3 = rho(2, 1)
      p4 = rho(2, 2)
   !   call SIBAF   ! calculation of probabilities per layer
c
      write(6,413)theta0,V,Vopt,p1,p2,p3,p4
  413 format(1P7E13.5)
c --> loop for material segments (117 for D2O)
c
  210 continue 
c
      STOP
c --> outputs
      WRITE(21,97)yt0
      WRITE(22,97)yt0
   97 FORMAT(12x,1P158E12.4)
      do 91 m=1,260
      WRITE(21,98)xdm(m),(array(m,l,1),l=1,158)
      WRITE(22,98)xdm(m),(array(m,l,2),l=1,158)  
   98 FORMAT(1P159E12.4)
   91 enddo
c
c      nol=0
c      do 513 mm=1,260
c      do 513 kt=1,158
c      v=array(mm,kt,1)
c      if(v.ge.9.E-5 .and. v.le.1.E-4)then
c      if(v.ge.0.018 .and. v.le.0.020)then
c      nol=nol+1
c      write(6,415)nol,xdm(mm),yt0(kt),v
c                                       endif
c  515 FORMAT(I8,1P3E13.5)
c  513 enddo
c
      CALL HUNPAK(15,specn,CH,1)
      dx=40.
      do 34 ihi=1,300
      xval=dx*(ihi-0.5)
      write(6,228)xval,specn(ihi)
 228  format(1P2E16.6)
  34  enddo
      call HPRINT(15)
c
      CALL HUNPAK(16,speca,CH,1)
      dx=0.1
      do 35 ihi=1,200
      xval=dx*(ihi-0.5)
      write(6,228)xval,speca(ihi)
  35  enddo
      call HPRINT(16)
c
      write(6,*)
      call cpu_time(runtime)
      WRITE(6,*)'run time (sec) = ',runtime
      stop
      end
c
c######################################################################
c
      subroutine SIBAF
c-----October 6, 2022 
c-----Berezhiani ANalitical FORmulas 
c-----from paper "N-N' in absorbing material" 
c-----Symmetry 2022, 1, 0. 
c-----https://doi.org/10.3390/sym1010000
c-----propagation of entangled nn' state in material with V and W
c-----wide range of theta_0 vs wide range dm
c
c --> U = V +/- Dm +/- muB will treat that as single real number U=V
c --> W=W1+v[m/s]*W2 in eV
c --> Dm=Delta_m =mn'-mn = e.g. 2.E-7 eV  
c --> eps = epsilon e.g. 2.E-9 eV e.g. in vacuum for theta_0=1.E-2
c 
c --> units: m, sec, eV
c
      Real theta0,TOF
      Real V,W,U1,U2,WW1,WW2,DE
      REAL     p1,p2,p3,p4
      REAL     Dm!,vre12,vim12,vre21,vim21,vre11,vre22
      COMPLEX i,zeta,zeta2,arg,cze,cze2,sze,sze2,H1,H2
      COMPLEX czestar,szestar,cze2star,sze2star,uroo,roo
      COMPLEX     snn,smn,smm,H1C,H2C,snnC,smnC,smmC
      !COMPLEX     rho11,rho22,rho12,rho21
      !COMPLEX     psi1t,psi2t,psi1tC,psi2tC
      COMPLEX     aee1,aee2,ee1,ee2,ee1C,ee2C

      common /input/Dm,theta0,V,W,del,nelsc,TOF
      common /output/p1,p2,p3,p4
c
      i=CMPLX(0.0,1.0)      !  sqrt(-1.) 
c
      eps=0.5*ABS(Dm)*tan(2.*theta0)! eps for givem Dm and theta0 
      arg=2.*eps/(V-i*W)
      zeta2=atan(arg) 
      zeta=zeta2/2.
c
      cze=cos(zeta)
      sze=sin(zeta)
      cze2=cze**2
      sze2=sze**2
c
      czestar=CONJG(cze)
      szestar=CONJG(sze)
      cze2star=CONJG(cze2)
      sze2star=CONJG(sze2)
c
      uroo=(V-i*W)**2 + 4.*eps**2
      roo =sqrt(uroo)
      if(REAL(zeta).lt.0.0)roo=-roo
      H1=0.5*(V-i*W+roo)
      H2=0.5*(V-i*W-roo)
      H1C=CONJG(H1)
      H2C=CONJG(H2)
c
      U1=REAL(H1)
      U2=REAL(H2)
      DE=U1-U2
      WW1=AIMAG(H1)
      WW2=AIMAG(H2)
c 
      aee1=-i*H1*TOF
      ee1=cexp(aee1) 
      aee2=-i*H2*TOF    
      ee2=cexp(aee2)
      ee1C=CONJG(ee1)
      ee2C=CONJG(ee2)
c
      snn =cze2*ee1+sze2*ee2
      snnC=cze2star*ee1C+sze2star*ee2C
      smn =cze*sze*(ee1-ee2)
      smnC=czestar*szestar*(ee1C-ee2C)
      smm =sze2*ee1+cze2*ee2
      smmC=sze2star*ee1C+cze2star*ee2C
c
c      psi1t =snn*psi10+smn*psi20
c      psi1tC=snnC*psi10+smnC*psi20
c      psi2t =smn*psi10+smm*psi20
c      psi2tC=smnC*psi10+smmC*psi20
c
      p1=REAL(snn*snnC)
      p2=REAL(smn*smnC)
      p3=p2
      p4=REAL(smm*smmC)
c
c      rho11=psi1t*psi1tC
c      rho22=psi2t*psi2tC
c      rho12=psi1t*psi2tC
c      rho21=psi2t*psi1tC
c
c      vre12=REAL(rho12)
c      vim12=AIMAG(rho12)
c      vre21=REAL(rho21)
c      vim21=AIMAG(rho21)
c      vre11=REAL(rho11)
c      vre22=REAL(rho22)
c
      return
      end  
c
c#######################################################################
!c
      !subroutine MAXW(TSOUR,V)
!c     ---------------------------------------
!c     Maxwellian distribution of velocities V at temperature TSOUR:
!c     output: in m/sec
!c
      !parameter (AK2M=1.648656E+4)
      !parameter (AK3M=2.472984E+4)
!c --> AK2M = 2*K_boltzman/m_n
!c --> AK3M = 3*K_boltzman/m_n
!c
      !VPFT=6637.43                 !  for 342K thermal source
      !if(TSOUR.lt.273.)VPFT=1.E+8  !  for cold(er) source
!c --> VPFT - [m/sec] two regions of thermal flux parametrization
!c --> if the temperature of the source is constant:
      !TEMPE=TSOUR
!c --> temperature is a radial function:
!c --> (for cold source) NOT INITIATED
!c     REMIS=sqrt(X0**2+Y0**2)
!c     TEMPE=TVSR(REMIS)
!c
!c --> Used FLUX spectrum is divided by factor 10^11;
!c     V3MP2 - most probable velocity square (m/s)**2
      !V3MP2=AK3M*TEMPE
      !V2MP2=AK2M*TEMPE
      !VMP=sqrt(V3MP2)
      !AMAX=45000.*EXP(-1.5)/VMP
 !1    continue
      !V=10.*VMP*RNDM(1.)
      !A=AMAX*RNDM(1.)
      !if(V.lt.VPFT)then           ! #1
      !AF=20000.*exp(-V**2/V2MP2)*V**3/V2MP2**2
                   !else           ! #1
      !AF=493.597/V
                   !endif          ! #1
      !if(A.GT.AF)go to 1
!c
      !return
      !end

