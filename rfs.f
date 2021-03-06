
      program rfs
c
c   ***********************************************************************
c   this program calculates the complex reflection coefficient and field
c   profile of an axisymmetric rf structure, in response to a constant
c   amplitude wave of known frequency continuously injected from the
c   left boundary into the rf structure.
c
c   the rf structure is formed of multiple sections of uniform and
c   linearly tapered waveguides with resistive walls. examples are the
c   open cavity of a gyromonotron and the severed interaction structure
c   of a gyro-twt.
c
c   by assumption, the radius of the (circular cross-section) rf
c   structure is either constant or a slowly varying function of the
c   axial position z. the left end section of the rf structure is a
c   propagating waveguide of uniform cross-section and resistivity.
c   the right end section is either a propagating or a cut-off waveguide,
c   also of uniform cross-section and resistivity.
c
c   a single te(m,n) mode is assumed throughout the structure.
c
c   the program is based on the formalism in "spectral domain analysis
c   of open cavities" (by k. r. chu, et al.), int. j. infrared and
c   millimeter waves, vol. 13, no. 10, pp. 1571-1598, 1992.
c
c   note: comment statements in the main program beginning with ccc 
c         call for action by the user.
c
c   date of first version: april, 1991
c   date of this version: March 1, 1999
c   ************************************************************************
c
      implicit real (a,b,d-h,j-z), complex (c)
c   in the main program and subprograms cbc, difeq, radius, rho, and 
c   closs, all variables beginning with a,b,d-h,j-z are real numbers, 
c   all variables beginning with i are integers, and all variables 
c   beginning with c are complex numbers.
      dimension axmn(9,8),cg(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(100001),arw(100001),anrho(100001),afamp(100001),
     &afphse(100001)
      dimension afreq(1001),afmax(1001),agamp(1001)
      common/const/ci,pi,realc
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,xmn,im
      common/diag/az,arw,anrho,afamp,afphse,idiag,icont
      external cbc
      data axmn/
     & 3.832e0, 1.841e0, 3.054e0, 4.201e0, 5.318e0, 6.416e0,
     & 7.501e0, 8.578e0, 9.647e0,
     & 7.016e0, 5.331e0, 6.706e0, 8.015e0, 9.282e0,10.520e0,
     &11.735e0,12.932e0,14.116e0,
     &10.174e0, 8.536e0, 9.970e0,11.346e0,12.682e0,13.987e0,
     &15.268e0,16.529e0,17.774e0,
     &13.324e0,11.706e0,13.170e0,14.586e0,15.964e0,17.313e0,
     &18.637e0,19.942e0,21.229e0,
     &16.471e0,14.864e0,16.348e0,17.789e0,19.196e0,20.576e0,
     &21.932e0,23.268e0,24.587e0,
     &19.616e0,18.016e0,19.513e0,20.973e0,22.401e0,23.804e0,
     &25.184e0,26.545e0,27.889e0,
     &22.760e0,21.164e0,22.672e0,24.145e0,25.590e0,27.010e0,
     &28.410e0,29.791e0,31.155e0,
     &25.904e0,24.311e0,25.826e0,27.310e0,28.768e0,30.203e0,
     &31.618e0,33.015e0,34.397e0/
c axmn(iabs(im)+1,in) is the n-th non-zero root of jm'(x)=0 up to im=8 and in=8
c
c   universal constants
      ci=cmplx(0.0,1.0)
      pi=3.1415926
      realc=2.99792e10
c
ccc instruction for plotting within do loop 100 (1:plot, 0:do not plot)
      iplot=1
c
ccc instruction for writing within do loop 100 (1:write, 0:do not write)
      iwrite=1

ccc specify no of points on the z-axis to be used to mark the positions
c   of the left/right boundaries and the junctions between sections. 
      izmark=4
c   
c   rf structure dimension arrays zmark rwl,and rwr specified below are
c   to be passed to the subprograms through the common block. the array 
c   elements must be in unit of cm.
c
c   zmark is a z-coordinate array to mark, from left to right, the posi- 
c   tions of the left boundary (zmark(1)), the junctions between sections 
c   (zmark(2),...,zmark(izmark-1)), and the right boundary (zmark(izmark)).
c
c   zmark(1) must be located within the left end section and zmark(izmark)  
c   must be located within the right end section. the end section (on the 
c   left or right) is assumed to be either a progagating or a cut-off 
c   waveguide of uniform cross-section and resistivity.
c
c   in theory, lengths of the uniform end sections do not affect the results.
c   in practice, set the length of the cut-off end section (if any) suffi- 
c   ciently short to avoid numerical difficulties due to the exponential 
c   growth or attenuation of f with z.

c   rwl(i) is the wall radius immediately to the left of z=zmark(i). 
c   rwr(i) is the wall radius immediately to the right of z=zmark(i). 
c   wall radius between zmark(i) and zmark(i+1) will be linearly
c   interpolated between rwr(i) and rwl(i+1) by function radius(z).
c
ccc specify rf structure dimension arrays zmark, rwl, and rwr in cm. 
      r=0.9
      r2=1.1
      l=11.7
      l2=5.0
      theta=10.0
      zmark(1)=0.0
      zmark(2)=zmark(1)+l
      zmark(3)=zmark(2)+(r2-r)/tan(theta*pi/180.0)
      zmark(4)=zmark(3)+l2
      rwl(1)=r
      rwr(1)=r
      rwl(2)=r
      rwr(2)=r
      rwl(3)=r2
      rwr(3)=r2
      rwl(4)=r2
      rwr(4)=r2

c   wall resistivity arrays rhol and rhor specified below are to be
c   passed to the subprograms through the common block. array elements
c   must be in mks unit of ohm-m (1 ohm-m = 100*ohm-cm).
c   (example: resistivity of copper at room temperature=1.72e-8 ohm-m)
c
c   rhol(i) is the wall resistivity immediately to the left of z=zmark(i).
c   rhor(i) is the wall resistivity immediately to the right of z=zmark(i).
c   wall resistivity between zmark(i) and zmark(i+1) will be linearly
c   interpolated between rhor(i) and rhol(i+1) by function rho(z).
c
ccc specify wall resistivity arrays rhol and rhor in ohm-m. 
      rhocu=1.72e-8
      rhomks=0.0
      rhol(1)=rhomks
      rhor(1)=rhomks
      rhol(2)=rhomks
      rhor(2)=rhomks
      rhol(3)=rhomks
      rhor(3)=rhomks
      rhol(4)=rhomks
      rhor(4)=rhomks
c   write rf structure dimensions and resistivity
      write(*,2) r,r2,l,l2,theta
    2 format(/' rf structure dimensions (length in cm):',
     &/' r=',1pe10.3,', r2=',1pe10.3,', l=',1pe10.3,
     &', l2=',1pe10.3/' theta=',1pe10.3,' degree')
      write(*,3) izmark
    3 format(/' rf structure dimension arrays: (izmark=',i4,')')
      do 50 i=1,izmark,6
      if(izmark.ge.(i+5)) imax=i+5
      if(izmark.lt.(i+5)) imax=izmark
      write(*,4) (zmark(ii),ii=i,imax)
    4 format(/'   zmark(cm)=',6(1pe10.3,','))
      write(*,5) (rwl(ii),ii=i,imax)
    5 format('     rwl(cm)=',6(1pe10.3,','))
      write(*,6) (rwr(ii),ii=i,imax)
    6 format('     rwr(cm)=',6(1pe10.3,','))
   50 continue
      write(*,7) izmark
    7 format(/' wall resistivity arrays: (izmark=',i4,')')
      do 60 i=1,izmark,6
      if(izmark.ge.(i+5)) imax=i+5
      if(izmark.lt.(i+5)) imax=izmark
      write(*,8) (zmark(ii),ii=i,imax)
    8 format(/'   zmark(cm)=',6(1pe10.3,','))
      write(*,9) (rhol(ii),ii=i,imax)
    9 format(' rhol(ohm-m)=',6(1pe10.3,','))
      write(*,10) (rhor(ii),ii=i,imax)
   10 format(' rhor(ohm-m)=',6(1pe10.3,','))
   60 continue
      if(iplot.eq.1) go to 53
      l=zmark(izmark)-zmark(1)
      delz=l/float(1000)
      az(1)=zmark(1)
      arw(1)=radius(az(1))
      anrho(1)=rho(az(1))/1.72e-8
      rwmax=arw(1)
      nrhomx=anrho(1)
      do 54 i=1,1000
      az(i+1)=zmark(1)+delz*dfloat(i)
      arw(i+1)=radius(az(i+1))
      anrho(i+1)=rho(az(i+1))/1.72e-8
      if(arw(i+1).gt.rwmax) rwmax=arw(i+1)
      if(anrho(i+1).gt.nrhomx) nrhomx=anrho(i+1)
   54 continue
      call sscale(zmark(1),zmark(izmark),0.0e0,rwmax)
      call splot(az,arw,1001,'rr  rw(cm) vs z(cm)$')
      if(nrhomx.eq.0.0e0) go to 53
      call sscale(zmark(1),zmark(izmark),0.0e0,nrhomx)
      call splot(az,anrho,1001,'**  rho(ohm-m)/1.72e-8 vs z(cm)$')
   53 continue

ccc specify no of steps for z-integration. 
c   (for a new problem, always check convergence with respect to izstep).
      izstep=2000

ccc specify m and n of te(m,n) mode
      im=1
      in=1
      xmn=axmn(iabs(im)+1,in)
      write(*,11) im,in,xmn
   11 format(/' mode: im=',i3,', in=',i3,', xmn=',1pe11.4/)

      gampmx=0.0
      fampmx=0.0
ccc specify incident wave frequency.  
      iwmax=3
      do 100 iw=1,iwmax
      if(iw.eq.1) fhz=9.83e9
      if(iw.eq.2) fhz=10.05e9
      if(iw.eq.3) fhz=10.39e9
      afreq(iw)=fhz
      w=2.0*pi*fhz
c   angular frequency of the incident wave (in radian/sec, denoted by w) 
c   is passed to subprograms through the common block.

ccc guess complex reflection coefficient  
      cguess=cmplx(0.5,0.5)
      cg(1)=cguess
c   the complex reflection coefficient (denoted by variable cg(1)) has been
c   formulated to be the root of the equation cbc(cg)=0. cg(1) defined above 
c   is a guessed value to be supplied to subroutine muller as a starting 
c   value for the search of the correct reflection coefficient (referred to 
c   the z=0 plane).

c   calculate reflection coefficient (referred to the z=0 plane)
      ep1=1.0e-6
      ep2=ep1
      iroot=1
      imaxit=50
      icont=0
      idiag=0
      call muller(0,iroot,cg,imaxit,ep1,ep2,cbc,.false.)
      idiag=1
      value=cabs(cbc(cg(1)))
c   muller is called to solve equation cbc(cg)=0 for the complex root cg.
c   'idiag' is an instruction for diagnostics. if (and only if) idiag=1,
c   function cbc will store the profiles of wall radius, resistivity, and
c   rf field, and pass them through the common block to the main program
c   for plotting.
c   'icont' monitors the no of times that function cbc is called for root 
c   searching ('icont' goes up by one upon each call to cbc).
c   'value' monitors the accuracy of the root returned by muller.
      gr=cg(1)
      gi=-ci*cg(1)
      gamp=cabs(cg(1))
      gphase=atan2(gi,gr)
      agamp(iw)=gamp
      if(agamp(iw).gt.gampmx) gampmx=agamp(iw)

c   find rwmax, rhomax, and fmax
      rwmax=arw(1)
      rhomax=anrho(1)
      fmax=afamp(1)
      do 70 iz=1,izstep
      if(arw(iz+1).gt.rwmax) rwmax=arw(iz+1)
      if(anrho(iz+1).gt.rhomax) rhomax=anrho(iz+1)
      if(afamp(iz+1).gt.fmax) fmax=afamp(iz+1)
  70  continue
      afmax(iw)=fmax
      if(afmax(iw).gt.fampmx) fampmx=afmax(iw)

c   plot rf structure shape, wall resistivity, and rf field profile
      if(iplot.ne.1) go to 72
      call sscale(zmark(1),zmark(izmark),0.0,rwmax)
      call splot(az,arw,izstep+1,'rr  rw(cm) vs z(cm)$')
      if(rhomax.eq.0.0) go to 71
      call sscale(zmark(1),zmark(izmark),0.0,rhomax)
      call splot(az,anrho,izstep+1,'**  rho(ohm-m)/1.72e-8 vs z(cm)$')
   71 continue
      call sscale(zmark(1),zmark(izmark),0.0,fmax)
      call splot(az,afamp,izstep+1,'ff  amplitude of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),-pi,pi)
      call splot(az,afphse,izstep+1,'pp  phase(radian) of f vs z(cm)$')
      write(*,12)
   12 format(/)
   72 continue

      if(iwrite.ne.1) go to 99
c   write results
      write(*,13) fhz,gamp,gphase,icont,value,fmax,afamp(izstep+1),
     &izstep,ep1
   13 format(' freq=',1pe11.4,' hz, gamp=',1pe11.4,
     &', gphase(radian)=',1pe11.4,'(',i3,',',1pe8.1,')'/
     &' fmax=',1pe10.3,', famp at right boundary=',1pe10.3,
     &', izstep=',i5,', ep1=',1pe8.1)
      write(*,14)
   14 format('**********************************************************
     &**********************')
   99 continue
  100 continue
      if(iwmax.eq.1) go to 101
      call sscale(afreq(1),afreq(iwmax),0.0e0,gampmx)
      call splot(afreq,agamp,iwmax,
     &'gg  amplitude of reflection coefficient vs freq.(hz)$')
      call sscale(afreq(1),afreq(iwmax),0.0e0,fampmx)
      call splot(afreq,afmax,iwmax,
     &'mm  spatial maximum of famp vs freq.(hz)$')
  101 continue
      stop
      end
c ******************************************************************
      function cbc(cg)
      implicit real (a,b,d-h,j-z), complex (c)
      dimension y(20),dy(20),q(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(100001),arw(100001),anrho(100001),afamp(100001),
     &afphse(100001)
      common/const/ci,pi,realc
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,xmn,im
      common/diag/az,arw,anrho,afamp,afphse,idiag,icont
      external difeq
      z1=zmark(1)
      rw=radius(z1)
      if((w/realc)**2.le.xmn**2/rw**2) write(*,1)
      if((w/realc)**2.le.xmn**2/rw**2) stop 100
    1 format(' *** w below cutoff in first section ***')
c
c initialize rf field array y at left boundary.
      ckz2=(w/realc)**2-xmn**2*closs(z1)/rw**2
      ckz=csqrt(ckz2)
      cdum1=cexp(ci*ckz*z1)
      cdum2=cexp(-ci*ckz*z1)
      cdum3=cdum1+cg*cdum2
      cdum4=cdum1-cg*cdum2
      y(1)=cdum3
      y(2)=-ci*cdum3
      y(3)=ci*ckz*cdum4
      y(4)=ckz*cdum4
      y(5)=z1
c store radius, resistivity (normalized to that of copper), and rf field at 
c left boundary
      az(1)=z1
      arw(1)=radius(z1)
      anrho(1)=rho(z1)/1.72e-8
      afamp(1)=sqrt(y(1)**2+y(2)**2)
      afphse(1)=atan2(y(2),y(1))
c integrate differential equations for rf field array y from left boundary   
c to right boundary.
      do 10 i=1,5
   10 q(i)=0.0
      l=zmark(izmark)-z1
      delz=l/float(izstep)
      do 100 iz=1,izstep
c advance rf field array y by one step.
      call rkint(difeq,y,dy,q,1,5,delz)
      if(idiag.ne.1) go to 99
c store radius, resistivity (normalized to that of copper), and rf field 
c as functions of z
c (arrays az, arw, anrho, afamp, and afphse are passed to the main program
c through the common block).
      z=y(5)
      az(iz+1)=z
      arw(iz+1)=radius(z)
      anrho(iz+1)=rho(z)/1.72e-8
      afamp(iz+1)=sqrt(y(1)**2+y(2)**2)
      afphse(iz+1)=atan2(y(2),y(1))
   99 continue
  100 continue
c calculate cbc(cg).
      fr=y(1)
      fi=y(2)
      fpr=y(3)
      fpi=y(4)
      cf=cmplx(fr,fi)
      cfp=cmplx(fpr,fpi)
      zf=y(5)
      rwf=radius(zf)
      if((w/realc)**2.le.xmn**2/rwf**2) go to 202
      ckz2=(w/realc)**2-xmn**2*closs(zf)/rwf**2
      ckz=csqrt(ckz2)
      cbc=cfp-ci*ckz*cf
      go to 201
  202 continue
      ckapa2=xmn**2*closs(zf)/rwf**2-(w/realc)**2
      ckapa=csqrt(ckapa2)
      cbc=cfp+ckapa*cf
  201 continue
      icont=icont+1
      return
      end
c ******************************************************************
      subroutine difeq(y,dy,ieqfst,ieqlst)
c this subroutine is called by subroutine rkint to calculate the 
c derivatives of rf field array y (defined as array dy). subroutine
c rkint is called by function cbc to integrate the differential
c equations for the rf field array.
      implicit real (a,b,d-h,j-z), complex (c)
      dimension y(20),dy(20)
      common/const/ci,pi,realc
      common/mode/w,xmn,im
      z=y(5)
      rw=radius(z)
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      kz2r=ckz2
      kz2i=-ci*ckz2
      dy(1)=y(3)
      dy(2)=y(4)
      dy(3)=-kz2r*y(1)+kz2i*y(2)
      dy(4)=-kz2i*y(1)-kz2r*y(2)
      dy(5)=1.0
      return
      end
c **********************************************************************
      function radius(z)
c radius(z) is the wall radius at position z.
      implicit real (a,b,d-h,j-z), complex (c)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      imax=izmark-1
      do 100 i=1,imax
      if(z.ge.zmark(i).and.z.le.zmark(i+1))
     &radius=rwr(i)+(rwl(i+1)-rwr(i))*(z-zmark(i))/
     &(zmark(i+1)-zmark(i))
  100 continue
c default values
      if(z.lt.zmark(1)) radius=rwl(1)
      if(z.gt.zmark(izmark)) radius=rwr(izmark)
      return
      end
c **********************************************************************
      function rho(z)
c rho(z) is the wall resistivity at position z.
      implicit real (a,b,d-h,j-z), complex (c)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      imax=izmark-1
      do 100 i=1,imax
      if(z.ge.zmark(i).and.z.le.zmark(i+1))
     &rho=rhor(i)+(rhol(i+1)-rhor(i))*(z-zmark(i))/(zmark(i+1)-zmark(i))
  100 continue
c default values
      if(z.lt.zmark(1)) rho=rhol(1)
      if(z.gt.zmark(izmark)) rho=rhor(izmark)
      return
      end
c **********************************************************************
      function closs(z)
c closs(z) is a factor to account for wall losses of the te(m,n) mode
c at position z. closs(z)=1.0 for zero wall resistivity.
      implicit real (a,b,d-h,j-z), complex (c)
      common/const/ci,pi,realc
      common/mode/w,xmn,im
      delta=sqrt(realc**2*rho(z)/(2.0*pi*w*9.0e9))
      rw=radius(z)
      wcmn=xmn*realc/rw
      cdum1=(1.0+ci)*delta/rw
      cdum2=(float(im**2)/(xmn**2-float(im**2)))*(w/wcmn)**2
      closs=1.0-cdum1*(1.0+cdum2)
      return
      end
c **********************************************************************
      subroutine rkint (derivy, y, dy, q, neqfst, neqlst, dx)
c ref. ralston & wilf "mathematical methods for digital computers",p.117
          real  y(neqlst), dy(neqlst), q(neqlst)
          real  a(4), b(4), c(4)
          real  dx,t
          external derivy
          data a /0.5e0, 0.29289322e0, 1.7071068e0, 0.16666667e0/,
     1         b /2.0e0, 1.0e0, 1.0e0, 2.0e0/,
     2         c /0.5e0, 0.29289322e0, 1.7071068e0, 0.5e0/
          do 1 j = 1, 4
          call derivy (y, dy, neqfst, neqlst)
          do 2 i = neqfst, neqlst
          t = a(j)*(dy(i) - b(j)*q(i))
          y(i) = y(i) + dx*t
          q(i) = q(i) + 3.0e0*t - c(j)*dy(i)
    2     continue
    1     continue
      return
      end
c ********************************************************************
      subroutine muller (kn,n,rts,maxit,ep1,ep2,fn,fnreal)
          implicit complex(a-h,o-z)
          complex num,lambda
           real  ep1,ep2,eps1,eps2
           real  abso
      external fn
           real  aimag
      logical fnreal
      dimension rts(6)
          aimag(x)=  (0.e0,-1.e0)*x
c
c this subroutine taken from elementary numerical analysis:
c                               an algorithmic approach
c                            by: conte and de boor
c                            algorithm 2.11 - muller's method
c
c initialization.
          eps1=amax1(ep1,1.e-12)
          eps2=amax1(ep2,1.e-20)
      ibeg = kn + 1
      iend = kn + n
c
      do 100 i = ibeg, iend
      kount = 0
c compute first three estimates for root as,
c     rts(i) + 0.5 , rts(i) - 0.5 , rts(i).
      abso=cabs(rts(i))
      if (abso) 11,12,11
   11 firsss=rts(i)/100.0e0
      go to 13
   12 firsss=cmplx(1.0e0,0.0e0)
   13 continue
      secdd=firsss/100.0e0
    1 h=.5e0*firsss
      rt = rts(i) + h
      assign 10 to nn
      go to 70
   10 delfpr = frtdef
      rt = rts(i) - h
      assign 20 to nn
      go to 70
   20 frtprv = frtdef
      delfpr = frtprv - delfpr
      rt = rts(i)
      assign 30 to nn
      go to 70
   30 assign 80 to nn
      lambda = -0.5
c compute next estimate for root.
   40 delf = frtdef - frtprv
      dfprlm = delfpr * lambda
      num = - frtdef * ( 1.0 + lambda ) * 2
      g = ( 1.0 + lambda * 2 ) * delf - lambda * dfprlm
      sqr = g * g + 2.0 * num * lambda * ( delf - dfprlm )
      if ( fnreal .and. real(sqr) .lt. 0.0 ) sqr = 0.0
      sqr = csqrt (sqr)
      den = g + sqr
      if ( real(g) * real(sqr) + aimag(g) * aimag(sqr) .lt. 0.0 )
     * den = g - sqr
      if ( cabs (den) .eq. 0.0 ) den = 1.0
      lambda = num / den
      frtprv = frtdef
      delfpr = delf
      h = h * lambda
      rt = rt + h
      if ( kount .gt. maxit ) write(*,1492)
 1492      format(' lack of convergence in muller')
      if ( kount .gt. maxit ) go to 100
c
   70 kount = kount + 1
      frt=fn   (      rt   )
      frtdef = frt
      if ( i .lt. 2 ) go to 75
c
      do 71 j = 2, i
      den = rt - rts( j - 1 )
      if ( cabs (den) .lt. eps2  ) go to 79
   71 frtdef = frtdef / den
   75 go to nn, (10,20,30,80)
   79 rts(i)=rt+ secdd
      go to 1
c check for convergence.
   80 if ( cabs (h) .lt. eps1 * cabs (rt) ) go to 100
      if ( amax1 ( cabs (frt), cabs (frtdef) ) .lt. eps2 ) go to 100
c
c check for divergence.
      if ( cabs (frtdef) .lt. 10.0 * cabs (frtprv) ) go to 40
      h = h / 2.0
      lambda = lambda / 2.0
      rt = rt - h
      go to 70
  100 rts(i) = rt
      return
      end
c *********************************************************************
      subroutine bplot(x, y, npt, tit)
      character   char,charr,tit(4),title(101),arr(51,101)
      character v,h,plus,b,etit
      dimension x(1), y(1), xx(101)
      logical terr
      character*80 teor
      data v,h,plus,b /'|','-','+',' '/, etit/'$'/, nc /0/
      data teor/'-error-    title must be, ((1-101)characters long) and
     1(terminated by $ )'/
      nscale = 0
      if (npt .eq. 0) go to  30
      if(nc.ne.0) go to 20
      nc=1
      call mxmn (x,npt,n1x,n0x)
      call mxmn (y,npt,n1y,n0y)
      call quan(x(n1x), x(n0x), xmax, xmin)
      call quan(y(n1y), y(n0y), ymax, ymin)
      go to 13
c
      entry      bscale (x0,x1, y0,y1)
      nscale = 1
      nc = 1
      xmax = x1
      xmin = x0
      ymax = y1
      ymin = y0
   13 dx = (xmax - xmin)*.01
      rdx = 1.0/dx
      dy = (ymax - ymin)*.02
      rdy = 1.0/dy
      do 15 n = 1,101,20
      xx(n) = xmin + (n-1)*dx
   15 if(abs(xx(n)).le. .001*abs((xmax-xmin)))xx(n) = 0.0
      do 50 i = 1,51
      charr = b
      if(mod(i-1 , 5) .eq. 0)charr = h
      do 55 j = 1,101
      char = charr
      if(mod(j-1 , 10) .ne. 0)go to 55
      char = v
      if(charr .ne. b)char = plus
   55 arr(i,j) = char
   50 continue
      if(nscale.eq.1) return
   20 do 16 n=1,npt
      nx = (x(n) - xmin)*rdx + 1.5
      ny = (ymax - y(n))*rdy + 1.5
      if(nx .lt.   1)go to 16
      if(ny .lt.   1)go to 16
      if(ny .gt.  51)go to 16
      if(nx .gt. 101)go to 16
      arr(ny,nx) = tit(1)
   16 continue
   30 if(tit(2) .eq. b .and. tit(3) .eq. b .and. tit(4) .eq. b)return
      terr = .false.
      do 70 i = 2,102
      ntit = i - 1
      if(tit(i+4) .eq. etit)go to 75
   70 continue
      terr = .true.
      ntit = 73
   75 mspc = (101 - ntit)/2
      nspc = mspc + ntit + 1
      do 80 i=1,101
      title(i) = b
      if(i .le. mspc .or. i .ge. nspc)go to 80
      title(i) = tit(i-mspc+4)
      if(terr)title(i) = teor(i-mspc:i-mspc)
   80 continue
      nc=0
      write(*,1) title
    1 format(1h1,///,15x,101a1)
      do 17 ny=1,51
      if(mod(ny-1 , 5) .ne. 0)go to 65
      ay = ymax-(ny-1)*dy
      if(abs(ay) .le. .001*abs((ymax-ymin)))ay = 0.0
      write(*, 9) ay, (arr(ny,n),n = 1,101)
      go to 17
   65 write(*,10)     (arr(ny,n),n = 1,101)
   17 continue
    9 format(1pe14.4,1x,101a1)
   10 format(15x,101a1)
      write(*,12)  (xx(n),n=1,101,20)
   12 format(1p6e20.4)
      return
      end
c *********************************************************************
      subroutine splot(x, y, npt, tit)
      character   char,charr,tit(4),title(101),arr(51,101)
      character v,h,plus,b,etit
      dimension x(1), y(1), xx(101)
      logical terr
      character*80 teor
      data v,h,plus,b /'|','-','+',' '/, etit/'$'/, nc /0/
      data teor/'-error-   title must be, ((1-nxap)characters long) and
     1(terminated by $ )'/
      nscale = 0
      if (npt .eq. 0) go to  30
      if(nc.ne.0) go to 20
      nc=1
      call mxmn (x,npt,n1x,n0x)
      call mxmn (y,npt,n1y,n0y)
      call quan(x(n1x), x(n0x), xmax, xmin)
      call quan(y(n1y), y(n0y), ymax, ymin)
      go to 13
c
      entry      sscale (x0,x1, y0,y1)
      nxr=10
      nxp=6
      nyr=5
      nyp=4
      nxap=nxr*nxp+1
      nyap=nyr*nyp+1
      nscale = 1
      nc = 1
      xmax = x1
      xmin = x0
      ymax = y1
      ymin = y0
   13 dx = (xmax - xmin)/(nxap-1.0)
      rdx = 1.0/dx
      dy = (ymax - ymin)/(nyap-1.0)
      rdy = 1.0/dy
      nxp2=nxp*2
      do 15 n = 1,nxap,nxp2
      xx(n) = xmin + (n-1)*dx
   15 if(abs(xx(n)).le. .001*abs((xmax-xmin)))xx(n) = 0.0
      do 50 i = 1,nyap
      charr = b
      if(mod(i-1 , nyp) .eq. 0)charr = h
      do 55 j = 1,nxap
      char = charr
      if(mod(j-1 , nxp) .ne. 0)go to 55
      char = v
      if(charr .ne. b)char = plus
   55 arr(i,j) = char
   50 continue
      if(nscale.eq.1) return
   20 do 16 n=1,npt
      nx = (x(n) - xmin)*rdx + 1.5
      ny = (ymax - y(n))*rdy + 1.5
      if(nx .lt.   1)go to 16
      if(ny .lt.   1)go to 16
      if(ny .gt.  nyap)go to 16
      if(nx .gt. nxap)go to 16
      arr(ny,nx) = tit(1)
   16 continue
   30 if(tit(2) .eq. b .and. tit(3) .eq. b .and. tit(4) .eq. b)return
      terr = .false.
      nxap1=nxap+1
      do 70 i = 2,nxap1
      ntit = i - 1
      if(tit(i+4) .eq. etit)go to 75
   70 continue
      terr = .true.
      ntit = 73
   75 mspc = (nxap - ntit)/2
      nspc = mspc + ntit + 1
      do 80 i=1,nxap
      title(i) = b
      if(i .le. mspc .or. i .ge. nspc)go to 80
      title(i) = tit(i-mspc+4)
      if(terr)title(i) = teor(i-mspc:i-mspc)
   80 continue
      nc=0
      write(*,1) (title(it),it=1,nxap)
    1 format(1h0,/,12x,101a1)
      do 17 ny=1,nyap
      if(mod(ny-1 , nyp) .ne. 0)go to 65
      ay = ymax-(ny-1)*dy
      if(abs(ay) .le. .001*abs((ymax-ymin)))ay = 0.0
      write(*, 9) ay, (arr(ny,n),n = 1,nxap)
      go to 17
   65 write(*,10)     (arr(ny,n),n = 1,nxap)
   17 continue
    9 format(1x,1pe11.4,101a1)
   10 format(12x,101a1)
      write(*,12)  (xx(n),n=1,nxap,nxp2)
   12 format(1x,1pe17.4,1p5e12.4)
      return
      end
c ******************************************************************
      subroutine mxmn (y,j ,nmax, nmin)
      real y(1)
      real quant(5)
      data quant/10.,5.,2.5,2.,1./
      nmax=1
      nmin=1
      do 1 n=1,j
      if (y(n) .gt. y(nmax)) nmax=n
      if(y(n).lt.y(nmin)) nmin=n
    1 continue
      return
c
      entry quan(x1, x0, xmax, xmin)
c     provide conveniently quantized limits for graphing
c     fake out the optimizer
      dx = x1 - x0
      if(dx.eq.0.) dx = abs(x1)
      if(dx.eq.0.) dx = 1
    2 ix = alog10(dx)
      if(dx.lt.1.) ix = ix -1
      xch = 10.**ix
c     roundoff could make xman gt 10
      xman = amin1(10., dx/xch)
      do 18 n=1,5
   18 if(xman.le. quant(n)) qx = quant(n)*xch*.1
      a = .5*(x1 + x0)/qx
      if(a.ge.0.) midx = a + .5
      if(a.lt.0.) midx = a - .5
      xmax = (midx + 5)*qx
      xmin = (midx - 5) *qx
      dx = 1.1*dx
      if(x1.gt.xmax+.01*dx .or. x0.lt.xmin-.01*dx) go to 2
      return
      end

