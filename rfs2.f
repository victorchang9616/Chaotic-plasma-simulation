      program rfs2
c
c   (all in double precision)
c    ***********************************************************************
c    this program calculates the reflection/transmission characteristics, 
c    ohmic losses, and axial field profile of an axisymmetric rf structure, 
c    in response to an incident wave of constant power and known frequency
c    continuously injected from the left end into the rf structure.
c
c    'rfs2' is a double precision version of 'rfs' with more diagnostics.
c    for example, at a specified incident frequency and power, profiles of
c    the net wave power, the forward and backward wave powers, and the wall
c    losses are diagnosed and plotte as functions of the axial position z. 
c    return and transmission losses are also printed and plotted as functions
c    of the incident wave frequency.
c
c    the rf structure is formed of multiple sections of uniform and linearly
c    tapered waveguides with resistive walls.
c
c    by assumption, the radius of the (circular cross-section) rf structure
c    is either constant or a slowly varying function of the axial position z.
c    the left end section of the rf structure is a propagating waveguide of
c    uniform cross-section and resistivity. the right end section is either
c    a cut-off or a propagating waveguide, also of uniform crosss-section and
c    resistivity.
c
c    a single te(m,n) mode is assumed throughout the structure.
c
c    the program is based on the formalism in "spectral domain analysis
c    of open cavities" (by k. r. chu, et al.), int. j. infrared and
c    millimeter waves, vol. 13, no. 10, pp. 1571-1598, 1992.
c
c    note: comment statements in the main program beginning with ccc
c          call for action by the user.
c
c    date of first version: february, 1993
c    date of this version: March 1, 1999
c    ************************************************************************
c
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
c    in the main program and subprograms cbc, difeq, radius, rho, and
c    closs, all variables beginning with a,b,d-h,j-z are real numbers,
c    all variables beginning with i are integers, and all variables
c    beginning with c are complex numbers.
      dimension axmn(9,8),cg(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(100001),arw(100001),anrho(100001),afamp(100001),
     &afphse(100001),apnet(100001),apfwd(100001),apbwd(100001),
     &apohmz(100001),agfwd(100001),agbwd(100001),agohm(100001)
      dimension afreq(1001),atldb(1001),arldb(1001),afmax(1001),
     &agamp(1001)
      common/const/ci,pi,realc
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,im
      common/phase/phsin,phsout,delphs
      common/diag/az,arw,anrho,afamp,afphse,apnet,apfwd,apbwd,apohmz,
     &pohmsm,idiag,icont
      external cbc
      data axmn/
     & 3.832d0, 1.841d0, 3.054d0, 4.201d0, 5.318d0, 6.416d0,
     & 7.501d0, 8.578d0, 9.647d0,
     & 7.016d0, 5.331d0, 6.706d0, 8.015d0, 9.282d0,10.520d0,
     &11.735d0,12.932d0,14.116d0,
     &10.174d0, 8.536d0, 9.970d0,11.346d0,12.682d0,13.987d0,
     &15.268d0,16.529d0,17.774d0,
     &13.324d0,11.706d0,13.170d0,14.586d0,15.964d0,17.313d0,
     &18.637d0,19.942d0,21.229d0,
     &16.471d0,14.864d0,16.348d0,17.789d0,19.196d0,20.576d0,
     &21.932d0,23.268d0,24.587d0,
     &19.616d0,18.016d0,19.513d0,20.973d0,22.401d0,23.804d0,
     &25.184d0,26.545d0,27.889d0,
     &22.760d0,21.164d0,22.672d0,24.145d0,25.590d0,27.010d0,
     &28.410d0,29.791d0,31.155d0,
     &25.904d0,24.311d0,25.826d0,27.310d0,28.768d0,30.203d0,
     &31.618d0,33.015d0,34.397d0/
c axmn(iabs(im)+1,in) is the n-th non-zero root of jm'(x)=0 up to im=8 and in=8
c
c    universal constants
      ci=dcmplx(0.0,1.0)
      pi=3.1415926
      realc=2.99792d10
c
ccc instruction for plotting within do loop 100 (1:plot, 0:do not plot)
      iplot=0
c
ccc instruction for writing within do loop 100 (1:write, 0:do not write)
      iwrite=1

ccc  specify no of points on the z-axis to be used to mark the positions
c    of the left/right boundaries and the junctions between sections.
      izmark=10

c    rf structure dimension arrays zmark rwl,and rwr specified below are to
c    be passed to the subprograms through the common block. the array elements
c    must be in unit of cm.
c
c    zmark is a z-coordinate array to mark, from left to right, the posi-
c    tions of the left boundary (zmark(1)), the junctions between sections
c    (zmark(2),...,zmark(izmark-1)), and the right boundary (zmark(izmark)).
c
c    zmark(1) must be located within the left end section and zmark(izmark)
c    must be located within the right end section. the end section (on the
c    left or right) is assumed to be either a progagating or a cut-off
c    waveguide of uniform cross-section and resistivity.
c
c    in theory, lengths of the uniform end sections do not affect the results.
c    in practice, set the length of the cut-off end section (if any) suffi-
c    ciently short to avoid numerical difficulties due to the exponential
c    growth or attenuation of f with z.
c
c    rwl(i) is the wall radius immediately to the left of z=zmark(i).
c    rwr(i) is the wall radius immediately to the right of z=zmark(i).
c    wall radius between zmark(i) and zmark(i+1) will be linearly
c    interpolated between rwr(i) and rwl(i+1) by function radius(z).
c
ccc  specify reference wall radius (rwref) in cm.
      rwref=0.2654
c     rwref=0.2757
ccc  specify rf structure dimension arrays zmark, rwl, and rwr in cm.
      r1=0.3175
      r2=rwref
c     rsev=1.0*rwref
      rsev=0.2787
      r3=rwref
      r4=0.3175
      l1=2.0
      lt1=2.27
      l2=4.5
      lt2=0.17
      lsev=1.66
      lt3=lt2
      l3=l2
      lt4=lt1
      l4=l1

      zmark(1)=0.0
      zmark(2)=zmark(1)+l1
      zmark(3)=zmark(2)+lt1
      zmark(4)=zmark(3)+l2
      zmark(5)=zmark(4)+lt2
      zmark(6)=zmark(5)+lsev
      zmark(7)=zmark(6)+lt3
      zmark(8)=zmark(7)+l3
      zmark(9)=zmark(8)+lt4
      zmark(10)=zmark(9)+l4
      rwl(1)=r1
      rwr(1)=rwl(1)
      rwl(2)=r1
      rwr(2)=rwl(2)
      rwl(3)=r2
      rwr(3)=rwl(3)
      rwl(4)=r2
      rwr(4)=rwl(4)
      rwl(5)=rsev
      rwr(5)=rwl(5)
      rwl(6)=rsev
      rwr(6)=rwl(6)
      rwl(7)=r3
      rwr(7)=rwl(7)
      rwl(8)=r3
      rwr(8)=rwl(8)
      rwl(9)=r4
      rwr(9)=rwl(9)
      rwl(10)=r4
      rwr(10)=rwl(10)

c    wall resistivity arrays rhol and rhor specified below are to be
c    passed to the subprograms through the common block. array elements
c    must be in mks unit of ohm-m (1 ohm-m = 100*ohm-cm).
c    (example: resistivity of copper at room temperature=1.72e-8 ohm-m)
c
c    rhol(i) is the wall resistivity immediately to the left of z=zmark(i).
c    rhor(i) is the wall resistivity immediately to the right of z=zmark(i).
c    wall resistivity between zmark(i) and zmark(i+1) will be linearly
c    interpolated between rhor(i) and rhol(i+1) by function rho(z).
c
ccc  specify wall resistivity arrays rhol and rhor.
      rhocu=(1.72e-8)*1.0
      rhosev=rhocu*2.5e6
      rhol(1)=rhocu
      rhor(1)=rhocu
      rhol(2)=rhocu
      rhor(2)=rhocu
      rhol(3)=rhocu
      rhor(3)=rhocu
      rhol(4)=rhocu
      rhor(4)=rhocu 
      rhol(5)=rhosev
      rhor(5)=rhosev
      rhol(6)=rhosev
      rhor(6)=rhosev
      rhol(7)=rhocu 
      rhor(7)=rhocu
      rhol(8)=rhocu
      rhor(8)=rhocu
      rhol(9)=rhocu
      rhor(9)=rhocu
      rhol(10)=rhocu
      rhor(10)=rhocu

      write(*,2) r1,r2,rsev,r3,r4,
     &l1,lt1,l2,lt2,lsev,lt3,l3,lt4,l4
    2 format(' rf structure dimensions in cm: '/' (upper line: radii of
     &uniform sections, lower line: lengths of all sections)'
     &//1pe9.3,8x,1pe10.3,8x,1pe10.3,8x,1pe10.3,7x,1pe10.3,
     &/1pe8.2,8(1pe9.2))
      write(*,3) izmark
    3 format(/' rf structure dimension arrays: (izmark=',
     &i4,')')
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
      delz=l/dfloat(1000)
      az(1)=zmark(1)
      arw(1)=radius(az(1))
      anrho(1)=rho(az(1))/1.72d-8
      rwmax=arw(1)
      nrhomx=anrho(1)
      do 54 i=1,1000
      az(i+1)=zmark(1)+delz*dfloat(i)
      arw(i+1)=radius(az(i+1))
      anrho(i+1)=rho(az(i+1))/1.72d-8
      if(arw(i+1).gt.rwmax) rwmax=arw(i+1)
      if(anrho(i+1).gt.nrhomx) nrhomx=anrho(i+1)
   54 continue
      call sscale(zmark(1),zmark(izmark),0.0d0,rwmax)
      call splot(az,arw,1001,'rr  rw(cm) vs z(cm)$')
      if(nrhomx.eq.0.0d0) go to 53
      call sscale(zmark(1),zmark(izmark),0.0d0,nrhomx)
      call splot(az,anrho,1001,'**  rho(ohm-m)/1.72d-8 vs z(cm)$')
   53 continue

ccc  specify no of steps for z-integration.
c    (for a new problem, always check convergence with respect to izstep).
      izstep=1000

ccc  specify m and n of te(m,n) mode
      im=1
      in=1

      xmn=axmn(iabs(im)+1,in)
      kcapmn=bj(im,xmn)**2*(1.0-dfloat(im*im)/xmn**2)
      fcref=xmn*realc/(2.0*pi*rwref)
      write(*,11) im,in,xmn,kcapmn
   11 format(/' mode: im=',i3,', in=',i3,', xmn=',1pe11.4,
     &', kcapmn=',1pe10.3)
      write(*,12) rwref,fcref
   12 format(' rwref=reference wall radius=',1pd12.5,' cm'/
     &' fcref=cutoff frequency w.r.t. rwref=',1pd12.5,' hz'/)

ccc specify input (forward) wave power (in watts) and phase (named phsin,
c   in radian) at left end. phase of output wave (named phsout) and phase
c   difference between output and input waves (named delphs) will be
c   returned to the main program by function cbc through the common block.
      pin=1.0
      phsin=0.0

      gampmx=0.0
      tlmin=0.0
      rlmin=0.0
      fampmx=0.0
      mindb=0.0
ccc  specify incident wave frequency.
      iwmax=9
      do 100 iw=1,iwmax
      fhzmin=32.0d9
      fhzmax=38.0d9
      fhzmin=fcref*1.00001
      fhzmax=1.15d0*fcref
      if(iwmax.gt.1) fstep=(fhzmax-fhzmin)/dfloat(iwmax-1)
      if(iwmax.eq.1) fstep=0.0
      fhz=fhzmin+fstep*(iw-1)
      afreq(iw)=fhz
      w=2.0*pi*fhz
c    angular frequency of the incident wave (in radian/sec, denoted by w)
c    is passed to subprograms through the common block.
c
ccc  guess complex reflection coefficient
      cguess=dcmplx(0.5,0.5)
      if(iw.eq.1) cg(1)=cguess
c    the complex reflection coefficient (denoted by variable cg(1)) has been
c    formulated to be the root of the equation cbc(cg)=0. cg(1) defined above
c    is a guessed value to be supplied to subroutine muller as a starting
c    value for the search of the correct reflection coefficient.

c    calculate reflection coefficient (referred to the z=0 plane)
      ep1=1.0d-5
      ep2=ep1
      iroot=1
      imaxit=50
      icont=0
      idiag=0
      call muller(0,iroot,cg,imaxit,ep1,ep2,cbc,.false.)
      idiag=1
      value=cdabs(cbc(cg(1)))
c    muller is called to solve equation cbc(cg)=0 for the complex root cg.
c    'idiag' is an instruction for diagnostics. if (and only if) idiag=1,
c    function cbc will store the profiles of wall radius, resistivity, rf
c    field and powers, and pass them through the common block to the main
c    program for plotting.
c    'icont' monitors the no of times that function cbc is called for root
c    searching ('icont' goes up by one upon each call to cbc).
c    'value' monitors the accuracy of the root returned by muller.
      gr=cg(1)
      gi=-ci*cg(1)
      gamp=cdabs(cg(1))
      gphase=datan2(gi,gr)
      pconsv=(apnet(1)-apnet(izstep+1)-pohmsm)/pin
      tldb=10.0*dlog10(apfwd(izstep+1)/pin)
      rldb=10.0*dlog10(apbwd(1)/pin)
      if(pohmsm.ne.0.0) pohmdb=10.0*dlog10(pohmsm/pin)
      agamp(iw)=gamp
      atldb(iw)=tldb
      arldb(iw)=rldb
      if(agamp(iw).gt.gampmx) gampmx=agamp(iw)
      if(atldb(iw).lt.tlmin) tlmin=atldb(iw)
      if(arldb(iw).lt.rlmin) rlmin=arldb(iw)
      if(atldb(iw).lt.mindb) mindb=atldb(iw)
      if(arldb(iw).lt.mindb) mindb=arldb(iw)

      agfwd(1)=10.0d0*dlog10(apfwd(1)/pin)
      agbwd(1)=10.0d0*dlog10(apbwd(1)/pin)
      if(apohmz(1).ne.0.0d0) agohm(1)=10.0d0*dlog10(apohmz(1)/pin)
      if(apohmz(1).eq.0.0d0) agohm(1)=-1.0d50
      do 80 iz=1,izstep
      agfwd(iz+1)=10.0d0*dlog10(apfwd(iz+1)/pin)
      agbwd(iz+1)=10.0d0*dlog10(apbwd(iz+1)/pin)
      if(apohmz(iz+1).ne.0.0d0)
     &agohm(iz+1)=10.0d0*dlog10(apohmz(iz+1)/pin)
      if(apohmz(iz+1).eq.0.0d0) agohm(iz+1)=-1.0d50
   80 continue

c    find rwmax, nrhomx, fmax, pnetmx, pnetmn, gfwdmx, gfwdmn, gbwdmx, gbwdmn, 
c    gohmmx, gohmmn
      rwmax=arw(1)
      nrhomx=anrho(1)
      fmax=afamp(1)
      gfwdmx=agfwd(1)
      gfwdmn=agfwd(1)
      gbwdmx=agbwd(1)
      gbwdmn=agbwd(1)
      pnetmx=apnet(1)
      pnetmn=apnet(1)
      gohmmx=agohm(1)
      gohmmn=agohm(1)
      do 70 iz=1,izstep
      if(arw(iz+1).gt.rwmax) rwmax=arw(iz+1)
      if(anrho(iz+1).gt.nrhomx) nrhomx=anrho(iz+1)
      if(afamp(iz+1).gt.fmax) fmax=afamp(iz+1)
      if(agfwd(iz+1).gt.gfwdmx) gfwdmx=agfwd(iz+1)
      if(agfwd(iz+1).lt.gfwdmn) gfwdmn=agfwd(iz+1)
      if(agbwd(iz+1).gt.gbwdmx) gbwdmx=agbwd(iz+1)
      if(agbwd(iz+1).lt.gbwdmn) gbwdmn=agbwd(iz+1)
      if(apnet(iz+1).gt.pnetmx) pnetmx=apnet(iz+1)
      if(apnet(iz+1).le.pnetmn) pnetmn=apnet(iz+1)
      if(agohm(iz+1).gt.gohmmx) gohmmx=agohm(iz+1)
      if(agohm(iz+1).le.gohmmn) gohmmn=agohm(iz+1)
   70 continue
      afmax(iw)=fmax
      if(afmax(iw).gt.fampmx) fampmx=afmax(iw)
      if(pnetmn.ge.0.0d0) pnetmn=0.0d0
      if(gfwdmn.ge.0.0d0) gfwdmn=0.0d0
      if(gbwdmn.ge.0.0d0) gbwdmn=0.0d0
      if(gohmmn.ge.0.0d0) gohmmn=0.0d0
      if(gbwdmn.lt.(gbwdmx-30.0d0)) gbwdmn=gbwdmx-30.0d0
      if(gohmmn.lt.(gohmmx-30.0d0)) gohmmn=gohmmx-30.0d0

c    plot rf structure shape, wall resistivity, field and power profiles
      if(iplot.ne.1) go to 72
      call sscale(zmark(1),zmark(izmark),0.0d0,rwmax)
      call splot(az,arw,izstep+1,'rr  rw(cm) vs z(cm)$')
      if(nrhomx.eq.0.0) go to 71
      call sscale(zmark(1),zmark(izmark),0.0d0,nrhomx)
      call splot(az,anrho,izstep+1,'**  rho(ohm-m)/1.72d-8 vs z(cm)$')
   71 continue
      call sscale(zmark(1),zmark(izmark),0.0d0,fmax)
      call splot(az,afamp,izstep+1,'ff  amplitude of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),-pi,pi)
      call splot(az,afphse,izstep+1,'pp  phase(radian) of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),pnetmn,pnetmx)
      call splot(az,apnet,izstep+1,'nn  net wave power(w) vs z(cm)$') 
      write(*,21)
   21 format(/ ' note: forward wave power (hence the plot below) may be
     &inaccurately diagnosed'/'       in regions where w is slightly abo
     &ve cutoff or where the wall is tapered.'/
     &'       it was artificially set to 1.0d-50 watt where w < wcmn.')
      call sscale(zmark(1),zmark(izmark),gfwdmn,gfwdmx)
      call splot(az,agfwd,izstep+1,
     &'ff  10*log(forward wave power/input power) vs z(cm)$')
      if(pohmsm.eq.0.0d0) go to 73
      call sscale(zmark(1),zmark(izmark),gohmmn,gohmmx)
      call splot(az,agohm,izstep+1,
     &'oo  10*log(ohmic loss per cm/input power) vs z(cm)$')
   73 continue
      write(*,22)
   22 format(/ ' note: backward wave power (hence the plot below) may be
     & inaccurately diagnosed'/'       in regions where w is slightly ab
     &ove cutoff or where the wall is tapered.'/
     &'       it was artificially set to 1.0d-50 watt where w < wcmn.')
      call sscale(zmark(1),zmark(izmark),gbwdmn,gbwdmx)
      call splot(az,agbwd,izstep+1, 
     &'bb  10*log(backward wave power/input power) vs z(cm)$') 
      write(*,13)
   13 format(/)
   72 continue

      if(iwrite.ne.1) go to 99
c    write results
      write(*,14) fhz,gamp,gphase,icont,value,pin,fmax,izstep,ep1
   14 format(' freq=',1pe11.4,' hz, gamp=',1pe11.4,
     &', gphase(radian)=',1pe11.4,'(',i3,',',1pe8.1,')'/
     &' pin (incident power)=',1pe10.3,
     &' w, fmax=',1pe10.3,', izstep=',i6,', ep1=',1pe8.1)
      write(*,15) phsin,phsout,delphs
   15 format(' phsin(radian)=',1pd10.3,', phsout(radian)=',1pd10.3,
     &', delphs(radian)=',1pd10.3)
      write(*,16) apnet(1),apnet(izstep+1),apfwd(1),
     &apfwd(izstep+1),apbwd(1),apbwd(izstep+1),pohmsm,pconsv,
     &tldb,rldb
   16 format(' pnet (net wave power):      at left end=',1pe10.3,
     &' w, at right end=',1pe10.3,' w'/
     &' pfwd (forward wave power):  at left end=',1pe10.3,
     &' w, at right end=',1pe10.3,' w'/
     &' pbwd (backward wave power): at left end=',1pe10.3,
     &' w, at right end=',1pe10.3,' w'/
     &' pohmsm (sum of wall losses)=',1pe10.3,
     &' w, power conservation factor=',1pe10.3/
     &' 10*log(pfwd at right end/pin)= ',1pe10.3,
     &' db (transmission loss)'/
     &' 10*log(pbwd at left end/pin)=  ',1pe10.3,
     &' db (return loss)')
      if(pohmsm.ne.0.0) write(*,17) pohmdb
   17 format(' 10*log(pohmsm/pin)= ',1pe10.3,' db (ohmic loss)')
      write(*,18)
   18 format('**********************************************************
     &**********************')
   99 continue
  100 continue
      if(iwmax.eq.1) go to 101
      call sscale(afreq(1),afreq(iwmax),0.0d0,gampmx)
      call splot(afreq,agamp,iwmax,
     &'gg  amplitude of reflection coefficient vs freq.(hz)$')
      call sscale(afreq(1),afreq(iwmax),tlmin,0.0d0)
      call splot(afreq,atldb,iwmax,
     &'tt  transmission loss(db) vs freq.(hz)$')
      call sscale(afreq(1),afreq(iwmax),rlmin,0.0d0)
      call splot(afreq,arldb,iwmax,
     &'rr  return loss(db) vs freq.(hz)$')
      call sscale(afreq(1),afreq(iwmax),mindb,0.0d0)
      call splot(afreq,atldb,iwmax,'t   $')
      call splot(afreq,arldb,iwmax,
     &'rr  transmission and return losses(db) vs freq.(hz)$')
c     call sscale(afreq(1),afreq(iwmax),0.0d0,fampmx)
c     call splot(afreq,afmax,iwmax,
c    &'mm  spatial maximum of famp vs freq.(hz)$')
  101 continue
 1000 continue
      stop
      end
c ******************************************************************
      function cbc(cg)
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(20),dy(20),q(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(100001),arw(100001),anrho(100001),afamp(100001),
     &afphse(100001),apnet(100001),apfwd(100001),apbwd(100001),
     &apohmz(100001)
      common/const/ci,pi,realc
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,im
      common/phase/phsin,phsout,delphs
      common/diag/az,arw,anrho,afamp,afphse,apnet,apfwd,apbwd,apohmz,
     &pohmsm,idiag,icont
      external difeq
      z1=zmark(1)
      rw=radius(z1)
      if((w/realc)**2.le.xmn**2/rw**2) write(*,1)
      if((w/realc)**2.le.xmn**2/rw**2) stop 100
    1 format(' *** w below cutoff in first section ***')
c
c initialize rf field array y at left boundary according to the incident
c forward wave power (pin) specified in the main program, and store radius,
c resistivity (normalized to copper), rf field, and powers at left boundary.
c (wave power in w, wall loss per unit length in w/cm)
      ckz2=(w/realc)**2-xmn**2*closs(z1)/rw**2
      ckz=cdsqrt(ckz2)
      cdum1=cdexp(ci*ckz*z1)
      cdum2=cdexp(-ci*ckz*z1)
      cf0=dcmplx(1.0,0.0)*cdexp(ci*phsin)
      do 20 i=1,2
      cdum3=cf0*(cdum1+cg*cdum2)
      cdum4=cf0*(cdum1-cg*cdum2)
      y(1)=cdum3
      y(2)=-ci*cdum3
      y(3)=ci*ckz*cdum4
      y(4)=ckz*cdum4
      y(5)=0.0
      y(6)=z1
c
      cf=dcmplx(y(1),y(2))
      cfp=dcmplx(y(3),y(4))
      kzr=ckz
      cdum10=cdexp(ci*kzr*z1)
      cdum20=cdexp(-ci*kzr*z1)
      cdum5=ci*kzr*cf
      cdum6=2.0*ci*kzr
      dum17=w*xmn**2*kcapmn/8.0d7
      delta=dsqrt(realc**2*rho(z1)/(2.0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      cffwd=(cdum5+cfp)*cdum20/cdum6
      cfbwd=(cdum5-cfp)*cdum10/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z1)
c       the backward component of cf is cfbwd*cdexp(-ci*kzr*z1)
      az(1)=z1
      arw(1)=radius(z1)
      anrho(1)=rho(z1)/1.72d-8
      afamp(1)=dsqrt(y(1)**2+y(2)**2)
      afphse(1)=datan2(y(2),y(1))
      apfwd(1)=kzr*dum17*cdabs(cffwd)**2
      apbwd(1)=kzr*dum17*cdabs(cfbwd)**2
      cfstar=dcmplx(y(1),-y(2))
      cfpstr=dcmplx(y(3),-y(4))
      cdum30=cf*cfpstr-cfstar*cfp
      apnet(1)=0.5*ci*dum17*cdum30
      apohmz(1)=dum17*dum18*dum19*(y(1)**2+y(2)**2)
      cf0=cf0*dsqrt(pin/apfwd(1))
   20 continue
c integrate differential equations for rf field array y from left boundary
c to right boundary.
      do 10 i=1,6
   10 q(i)=0.0
      l=zmark(izmark)-z1
      delz=l/dfloat(izstep)
      do 100 iz=1,izstep
c advance rf field array y by one step.
      call rkint(difeq,y,dy,q,1,6,delz)
      if(idiag.ne.1) go to 99
c store radius, resistivity (normalized to copper), rf field, and powers 
c as functions of z (wave power in w, wall loss per unit length in w/cm).
c forward wave power (apfwd) and backward wave power (apbwd) are
c artificially set to 1.0d-50 where w < wcmn.
c arrays az, arw, afamp, afphse, apfwd, apbwd, apnet, and apohmz are
c passed to the main program through the common block.
      z=y(6)
      rw=radius(z)
      cf=dcmplx(y(1),y(2))
      cfp=dcmplx(y(3),y(4))
      delta=dsqrt(realc**2*rho(z)/(2.0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      az(iz+1)=z
      arw(iz+1)=radius(z)
      anrho(iz+1)=rho(z)/1.72d-8
      afamp(iz+1)=dsqrt(y(1)**2+y(2)**2)
      afphse(iz+1)=datan2(y(2),y(1))
      cfstar=dcmplx(y(1),-y(2))
      cfpstr=dcmplx(y(3),-y(4))
      cdum30=cf*cfpstr-cfstar*cfp
      apnet(iz+1)=0.5*ci*dum17*cdum30
      apohmz(iz+1)=dum17*dum18*dum19*(y(1)**2+y(2)**2)
      if((w/realc)**2.lt.(xmn/rw)**2) go to 29
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      ckz=cdsqrt(ckz2)
      kzr=ckz
      cdum10=cdexp(ci*kzr*z)
      cdum20=cdexp(-ci*kzr*z)
      cdum5=ci*kzr*cf
      cdum6=2.0*ci*kzr
      cffwd=(cdum5+cfp)*cdum20/cdum6
      cfbwd=(cdum5-cfp)*cdum10/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z)
c       the backward component of cf is cfbwd*cdexp(-ci*kzr*z)
      apfwd(iz+1)=kzr*dum17*cdabs(cffwd)**2
      apbwd(iz+1)=kzr*dum17*cdabs(cfbwd)**2
      if(iz.ne.izstep) go to 39
      ffwdr=cffwd*cdum10
      ffwdi=-ci*cffwd*cdum10
      phsout=datan2(ffwdi,ffwdr)
c     phsout=afphse(iz+1)
      delphs=phsout-phsin
      go to 39
   29 continue
      apfwd(iz+1)=1.0d-50
      apbwd(iz+1)=1.0d-50
      if(iz.ne.izstep) go to 39
      phsout=1.0d10
      delphs=1.0d10
   39 continue
   99 continue
  100 continue
c store sum of wall losses (in w). pohmsm is passed to the main program
c through the common block.
      pohmsm=y(5)
c calculate cbc(cg).
      fr=y(1)
      fi=y(2)
      fpr=y(3)
      fpi=y(4)
      cf=dcmplx(fr,fi)
      cfp=dcmplx(fpr,fpi)
      zf=y(6)
      rwf=radius(zf)
      if((w/realc)**2.le.xmn**2/rwf**2) go to 202
      ckz2=(w/realc)**2-xmn**2*closs(zf)/rwf**2
      ckz=cdsqrt(ckz2)
      cbc=cfp-ci*ckz*cf
      go to 201
  202 continue
      ckapa2=xmn**2*closs(zf)/rwf**2-(w/realc)**2
      ckapa=cdsqrt(ckapa2)
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
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(20),dy(20)
      common/const/ci,pi,realc
      common/mode/w,pin,xmn,kcapmn,im
      z=y(6)
      rw=radius(z)
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      kz2r=ckz2
      kz2i=-ci*ckz2
      delta=dsqrt(realc**2*rho(z)/(2.0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum17=w*xmn**2*kcapmn/8.0d7
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      dy(1)=y(3)
      dy(2)=y(4)
      dy(3)=-kz2r*y(1)+kz2i*y(2)
      dy(4)=-kz2i*y(1)-kz2r*y(2)
      dy(5)=dum17*dum18*dum19*(y(1)**2+y(2)**2)
      dy(6)=1.0
      return
      end
c **********************************************************************
      function radius(z)
c radius(z) is the wall radius at position z.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
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
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
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
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      common/const/ci,pi,realc
      common/mode/w,pin,xmn,kcapmn,im
      delta=dsqrt(realc**2*rho(z)/(2.0*pi*w*9.0d9))
      rw=radius(z)
      wcmn=xmn*realc/rw
      cdum1=(1.0+ci)*delta/rw
      cdum2=(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      closs=1.0-cdum1*(1.0+cdum2)
      return
      end
c **********************************************************************
      subroutine rkint (derivy, y, dy, q, neqfst, neqlst, dx)
c ref. ralston & wilf "mathematical methods for digital computers",p.117
          real*8  y(neqlst), dy(neqlst), q(neqlst)
          real*8  a(4), b(4), c(4)
          real*8  dx,t
          external derivy
          data a /0.5d0, 0.29289322d0, 1.7071068d0, 0.16666667d0/,
     1         b /2.0d0, 1.0d0, 1.0d0, 2.0d0/,
     2         c /0.5d0, 0.29289322d0, 1.7071068d0, 0.5d0/
          do 1 j = 1, 4
          call derivy (y, dy, neqfst, neqlst)
          do 2 i = neqfst, neqlst
          t = a(j)*(dy(i) - b(j)*q(i))
          y(i) = y(i) + dx*t
          q(i) = q(i) + 3.0d0*t - c(j)*dy(i)
    2     continue
    1     continue
      return
      end
c ********************************************************************
      subroutine muller (kn,n,rts,maxit,ep1,ep2,fn,fnreal)
          implicit complex*16 (a-h,o-z)
          complex*16 num,lambda
           real*8  ep1,ep2,eps1,eps2
           real*8  abso
      external fn
           real*8  aimag
      logical fnreal
      dimension rts(6)
          aimag(x)=  (0.d0,-1.d0)*x
c
c this subroutine taken from elementary numerical analysis:
c                               an algorithmic approach
c                            by: conte and de boor
c                            algorithm 2.11 - muller's method
c
c initialization.
          eps1=dmax1(ep1,1.d-12)
          eps2=dmax1(ep2,1.d-20)
      ibeg = kn + 1
      iend = kn + n
c
      do 100 i = ibeg, iend
      kount = 0
c compute first three estimates for root as,
c     rts(i) + 0.5 , rts(i) - 0.5 , rts(i).
      abso=cdabs(rts(i))
      if (abso) 11,12,11
   11 firsss=rts(i)/100.0d0
      go to 13
   12 firsss=dcmplx(1.0d0,0.0d0)
   13 continue
      secdd=firsss/100.0d0
    1 h=.5d0*firsss
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
      sqr = cdsqrt (sqr)
      den = g + sqr
      if ( real(g) * real(sqr) + aimag(g) * aimag(sqr) .lt. 0.0 )
     * den = g - sqr
      if ( cdabs (den) .eq. 0.0 ) den = 1.0
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
      if ( cdabs (den) .lt. eps2  ) go to 79
   71 frtdef = frtdef / den
   75 go to nn, (10,20,30,80)
   79 rts(i)=rt+ secdd
      go to 1
c check for convergence.
   80 if ( cdabs (h) .lt. eps1 * cdabs (rt) ) go to 100
      if ( dmax1 ( cdabs (frt), cdabs (frtdef) ) .lt. eps2 ) go to 100
c
c check for divergence.
      if ( cdabs (frtdef) .lt. 10.0 * cdabs (frtprv) ) go to 40
      h = h / 2.0
      lambda = lambda / 2.0
      rt = rt - h
      go to 70
  100 rts(i) = rt
      return
      end
c *********************************************************************
      subroutine bplot(x, y, npt, tit)
      implicit real*8 (a-h,o-z)
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
   15 if(dabs(xx(n)).le. .001*dabs((xmax-xmin)))xx(n) = 0.0
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
      if(dabs(ay) .le. .001*dabs((ymax-ymin)))ay = 0.0
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
      implicit real*8 (a-h,o-z)
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
   15 if(dabs(xx(n)).le. .001*dabs((xmax-xmin)))xx(n) = 0.0
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
      if(dabs(ay) .le. .001*dabs((ymax-ymin)))ay = 0.0
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
      implicit real*8 (a,b,d-h,o-z), complex*16 (c)
      real*8 y(1)
      real*8 quant(5)
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
      if(dx.eq.0.) dx = dabs(x1)
      if(dx.eq.0.) dx = 1
    2 ix = dlog10(dx)
      if(dx.lt.1.) ix = ix -1
      xch = 10.**ix
c     roundoff could make xman gt 10
      xman = dmin1(10.0d0, dx/xch)
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
c ***********************************************************************
      function bj(in,x)
c
c bj: bessel function j(x) of integer order.
c in: positive or negative order of bj (-40< in <40 ).
c x : positive or negative argument of bj (x=0 or 1.0e-10<x <300 ).
c
      implicit real*8 (a,b,d-h,j-z)
      dimension t(1000)
      isign1=1
      isign2=1
c     sx=sngl(x)
c     modified by c.s.hsue
      xabs=dabs(x)
      inabs=iabs(in)
      call bes(inabs,xabs,0,sbj,t)
      if(in.lt.0) isign1=(-1)**inabs
      if(x.lt.0.0d0) isign2=(-1)**inabs
      bj=sbj*dfloat(isign1*isign2)
      return
      end
c *********************************************************************
      subroutine bes(no,x,kode,result,t)
c       bessel function
cd * * * * * * * * * * * * * mini documentation * * * * * * * * * * * *
cd
cd       program:  bes(n0,x,kode,result,t)
cd       classification:  c3.e
cd       author / major users:  g. gilbert and g. baker
cd           modified by max miller,asc, austin,texas.
cd       originator:  rcc/ focus library
cd       latest mods:  march 19,1976.
cd       reviewer:
cd       location:
cd       date:  march 19,1976.
cd       description / purpose:  to evaluate bessel functions j, i.
cd           for a given r*4 arguments x and a given integer n bes cal-
cd           culates j(x) or i(x) for m = 0,1, . . .,n.
cd       arguments:
cd           no    non negative integer order of the bessel function.
cd            x    argument  of bessel function.
cd            k    a switch.  for k = 0 j(x) is computed.
cd                 for k = 1 i(x) is computed.
cd           result   real  value of the bessel function.
cd            t    temporary storage array.  its value must be the
cd                 maximum  of either no or integer(2*abs(x) + 10.).
cd       language and limitations:  fortran.
cd       entry points:  bes.
cd
cd * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      implicit real*8 (a-h,o-z)
      dimension t(1000)
  107 format(55h1negative order not accepted in bessel function routine)
      klam = 1
      ko = no + 1
      if (x) 6,1,6
    1 if(no) 5,2,3
    2 t(ko) = 1.
      result = 1.
      return
    3 result = 0.
      t(1)=1.0
      do 100 i=1,no
  100 t(i+1)=0.0
      return
    5 write(6,107)
      stop 1
    6 if(no) 5,7,7
    7 if(kode) 8,9,8
    8 klam = klam + 1
    9 continue
      absx=dabs(x)
      newno=no
      if(absx.ge.0.1) go to 201
      s=dlog10(absx)
      is=int(s)
      newno=idint(20.0/dfloat(is))
      newno=iabs(newno)
      if(newno.ge.no) go to 201
      newno1=newno+1
      do 200 i=newno1,no
  200 t(i+1)=0.0
c
  201 jo = 2 * idint(absx)
      mo=min0(newno,no)
      if(mo-jo) 11,12,12
   11 mo = jo
   12 mo = mo + 11
      t(mo) = 0.
      lub = mo-1
      t(lub) = 1.d-25
      go to (23,51),klam
c
c calculate j(x)
   23 f = 2 * lub
      mo = mo - 3
      i2 = mo
   24 f = f - 2.
      t(i2 + 1) = f / x * t(i2+2) - t(i2+3)
      if(i2) 25,26,25
   25 i2 = i2 - 1
      go to 24
   26 sum = t(1)
      do 40 j =  3,mo,2
   40 sum = sum + 2. * t(j)
      f = 1./ sum
      do 50 j = 1,ko
   50 t(j) = t(j) * f
      result = t(ko)
      return
c
c calculate i(x)
   51 f = 2 * lub - 2
      mo = mo - 3
      i2 = mo
  511 t(i2+1) = f / x * t(i2+2) + t(i2+3)
      if(i2) 52,53,52
   52 i2 = i2 - 1
      f = f - 2.
      go to 511
   53 sum = t(1)
      do 70 j = 2,mo
   70 sum = sum + 2.* t(j)
      f = 1. / sum *dexp (x)
      do 80 j = 1,ko
   80 t(j) = t(j) * f
      result = t(ko)
      return
      end

