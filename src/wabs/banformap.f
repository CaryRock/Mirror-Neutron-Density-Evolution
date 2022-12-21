      PROGRAM banformap
c
c-----October 6, 2022 
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
      Real theta0,TOF
      Real V,W,U1,U2,WW1,WW2,DE
      REAL(16) Dm,vre12,vim12,vre21,vim21,vre11,vre22
      REAL(16) probav(4),xdm(260),yt0(158),array(260,158,4)
      REAL(4) specn(300)
c     array(Dm,theta0,1=ro11,2=ro22,3=Re ro12,4=Im ro 12)
      COMPLEX i,zeta,zeta2,arg,cze,cze2,sze,sze2,H1,H2
      COMPLEX czestar,szestar,cze2star,sze2star,uroo,roo
      COMPLEX(16) snn,smn,smm,H1C,H2C,snnC,smnC,smmC
      COMPLEX(16) psi10,psi20,rho11,rho22,rho12,rho21
      COMPLEX(16) psi1t,psi2t,psi1tC,psi2tC
      COMPLEX(8)  aee1,aee2,ee1,ee2,ee1C,ee2C
c
      common /input/Dm,theta0,V,W,TOF,psi10,psi20
      common /output/vre12,vim12,vre21,vim21,vre11,vre22
c
      COMMON/PAWC/H(50000)
      CALL HLIMIT(50000)
      CALL HBOOK1(15,'velocity spectrum m/s',300,0.,12000,0.)
c
c --> initial parameters of the heavy water layer
      i=CMPLX(0.0,1.0)      ! sqrt(-1.) 
      TSOUR=340.            ! Kelvin
      psi10=CMPLX(1.0,0.0)  ! initial neutron state
      psi20=CMPLX(0.0,0.0)  ! initial mirror state
      dlina=0.022           ! in m - thickness of layer
      hbar=6.582119569E-16  ! precise eV*s hbar
      Vopt=1.659E-7         ! for heavy water
c     Vopt=1.992E-7         ! for B4C
c     Vopt=5.877E-8         ! for Cd
      W1=2.9588E-15         ! for heavy water constant part (absorption)
      W2=0.0                ! correction
c     W2=2.1367E-14         ! for heavy water to be x by velocity
c     W1=6.102E-9           ! constant part of W for B4C
c     W2=2.397E-14          ! part of W to be x by velocity in m/s B4C
c     W1=8.4558E-9          ! constant part of W for Cd
c     W2=9.914E-15          ! part of W to be x by velocity in m/s Cd  
      NSTAT=10000
      NSTAT=1   
c
c --> loop on W=W1 values
      DO 1002 iwp=1,9
      DO 1003 iwm=1,9
      W=float(iwm)*10**(-16.+iwp)
c --> loop of Delta_m values
c      DO 200 kp=1,13
c      DO 200 km=1,20
c      xd20=0.05*(km-1)
c      Dm=10**(xd20+kp-10.)
c      indx=km+20*(kp-1)
c      xdm(indx)=Dm
c      V=Vopt-Dm
      DO 200 kp=1,1
      Dm=1.659E-7 
      indx=kp
      xdm(indx)=Dm
      V=Vopt-Dm
c
c --> loop of theta_0 values     
c      DO 210 jp=1,9
c      do 210 jm=1,20
c      xt20=0.05*(jm-1)
c      theta0=10**(xt20+jp-9.)
c      if(theta0 .gt. 0.785398)go to 210
c      indy=jm+20*(jp-1)
c      yt0(indy)=theta0
      DO 210 jp=1,3
      theta0=10**(-1.*jp)
      indy=jp
      yt0(indy)=theta0
c --> here Dm in eV and theta0 are prepared:
c --> now will start averaging over thermal velocity spectrum
      do 175 j0=1,4
      probav(j0)=0.0D0
  175 enddo
c
      do 177 j=1,NSTAT
      call MAXW(TSOUR,Vel)
c     Vel=2318.
      Vel=1000.
      call HF1(15,Vel,1.)
      TOF=dlina/(Vel*hbar)   !  TOF in units of hbar
c --> input parameters: Dm,theta0,V,W,TOF,psi10,psi20
c      W=W1+vel*W2
c      WRITE(6,*)Dm,theta0,V,W,TOF,psi10,psi20
      call SIBAF
c <-- output parameters: vre12,vim12,vre21,vim21,vre11,vre22
c --> array(Dm,theta0,1=ro11,2=ro22,3=Rero12,4=Imro12)
c      WRITE(6,*)vre12,vim12,vre21,vim21,vre11,vre22

      probav(1)=probav(1)+vre11
      probav(2)=probav(2)+vre22
      probav(3)=probav(3)+vre12
      probav(4)=probav(4)+vim12

  177 enddo
c --> array(Dm,theta0,1=ro11,ro22,Re ro12,Im ro 12)
c --> probav(4),xdm(260),yt0(158),array(260,158,4) 
      array(indx,indy,1)=probav(1)/float(NSTAT)           
      array(indx,indy,2)=probav(2)/float(NSTAT) 
      array(indx,indy,3)=probav(3)/float(NSTAT) 
      array(indx,indy,4)=probav(4)/float(NSTAT) 
c
  210 continue
  200 continue 
      absl=0.5*hbar*vel/W
      absp=1.-exp(-dlina/absl)
      write(6,577)W,((array(1,jp,jj),jj=1,2),jp=1,3),absp
  577 FORMAT(1P8E12.5)
 1003 continue
 1002 continue
      STOP
c
c --> output rho 11, rho 22
      WRITE(20,97)yt0
      WRITE(21,97)yt0
c   97 FORMAT(12x,1P158E12.4)
   97 FORMAT(12x,1P3E12.4)
c      do 91 m=1,260
      do 91 m=1,1
c      WRITE(20,98)xdm(m),(array(m,l,1),l=1,158)
c      WRITE(21,98)xdm(m),(array(m,l,2),l=1,158)  
c   98 FORMAT(1P159E12.4)
      WRITE(20,98)xdm(m),(array(m,l,1),l=1,3)
      WRITE(21,98)xdm(m),(array(m,l,2),l=1,3)  
   98 FORMAT(1P3E12.4)
   91 enddo
c
      CALL HUNPAK(15,specn,CH,1)
      dx=40.
      do 34 ihi=1,300
      xval=dx*(ihi-0.5)
      write(6,128)xval,specn(ihi)
 128  format(1P2E16.6)
  34  enddo
      call HPRINT(15)
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
c-----wide range of theta_0 vs wide range dm
c
c --> U = V +/- Dm +/- muB will treat that as single real number U=V
c --> W=W1+v[m/s]*W2 in eV
c --> Dm=Delta_m =mn'-mn = e.g. 2.E-7 eV  
c --> eps = epsilon e.g. 2.E-9 eV e.g. in vacuum for theta_0=1.E-2
c 
c --> units: m, sec, eV
      common /input/Dm,theta0,V,W,TOF,psi10,psi20
      common /output/vre12,vim12,vre21,vim21,vre11,vre22
c
      Real theta0,TOF
      Real V,W,U1,U2,WW1,WW2,DE
      REAL(16) Dm,vre12,vim12,vre21,vim21,vre11,vre22
      COMPLEX i,zeta,zeta2,arg,cze,cze2,sze,sze2,H1,H2
      COMPLEX czestar,szestar,cze2star,sze2star,uroo,roo
      COMPLEX(16) snn,smn,smm,H1C,H2C,snnC,smnC,smmC
      COMPLEX(16) psi10,psi20,rho11,rho22,rho12,rho21
      COMPLEX(16) psi1t,psi2t,psi1tC,psi2tC
      COMPLEX(8)  aee1,aee2,ee1,ee2,ee1C,ee2C
c
      i=CMPLX(0.0,1.0)      !  sqrt(-1.) 
c
      print *, "W, tht: ", W, theta0
      !print *, "Psis 1: ", psi10, psi20
      print *, "Rhos 1: ", vre11, vre22
      print *, ""

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
      ee1=cdexp(aee1) 
      aee2=-i*H2*TOF    
      ee2=cdexp(aee2)
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
      psi1t =snn*psi10+smn*psi20
      psi1tC=snnC*psi10+smnC*psi20
      psi2t =smn*psi10+smm*psi20
      psi2tC=smnC*psi10+smmC*psi20
c
      rho11=psi1t*psi1tC
      rho22=psi2t*psi2tC
      rho12=psi1t*psi2tC
      rho21=psi2t*psi1tC
c
      vre12=REAL(rho12)
      vim12=AIMAG(rho12)
      vre21=REAL(rho21)
      vim21=AIMAG(rho21)
      vre11=REAL(rho11)
      vre22=REAL(rho22)

      !print *, "Psis 2: ", psi10, psi20
      print *, "Rhos 2: ", vre11, vre22
      print *, ""
c
      !print *, "W, V: ", W, V
      
      if (int(theta0 * 100) .eq. 10 .and. int(vre11 * 100) .eq. 89) stop
      return
      end  
c
c#######################################################################
c
      subroutine MAXW(TSOUR,V)
c     ---------------------------------------
c     Maxwellian distribution of velocities V at temperature TSOUR:
c     output: in m/sec
c
      parameter (AK2M=1.648656E+4)
      parameter (AK3M=2.472984E+4)
c --> AK2M = 2*K_boltzman/m_n
c --> AK3M = 3*K_boltzman/m_n
c
      VPFT=6637.43                 !  for 342K thermal source
      if(TSOUR.lt.273.)VPFT=1.E+8  !  for cold(er) source
c --> VPFT - [m/sec] two regions of thermal flux parametrization
c --> if the temperature of the source is constant:
      TEMPE=TSOUR
c --> temperature is a radial function:
c --> (for cold source) NOT INITIATED
c     REMIS=sqrt(X0**2+Y0**2)
c     TEMPE=TVSR(REMIS)
c
c --> Used FLUX spectrum is divided by factor 10^11;
c     V3MP2 - most probable velocity square (m/s)**2
      V3MP2=AK3M*TEMPE
      V2MP2=AK2M*TEMPE
      VMP=sqrt(V3MP2)
      AMAX=45000.*EXP(-1.5)/VMP
 1    continue
      V=10.*VMP*RNDM(1.)
      A=AMAX*RNDM(1.)
      if(V.lt.VPFT)then           ! #1
      AF=20000.*exp(-V**2/V2MP2)*V**3/V2MP2**2
                   else           ! #1
      AF=493.597/V
                   endif          ! #1
      if(A.GT.AF)go to 1
c
      return
      end

