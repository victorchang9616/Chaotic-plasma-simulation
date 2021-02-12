      program gtsim
c   
c   (all in double precision)
c   ***********************************************************************
c
c   this program calculates the gain, growth rate, output power, efficiency,
c   reflection coefficient, etc. of a single-mode gyro-twt or carm amplifier  
c   by numerically tracing the electron orbits (in both linear and nonlinear 
c   regimes) using either the iterative or noniterative procedure.
c
c   the intercaction structure is formed of a single section or multiple
c   sections of uniform and linearly tapered waveguides with resistive
c   walls.
c
c   by assumption, the radius of the (circular cross-section) interaction
c   structure is either constant or a slowly varying function of the axial
c   position z.
c
c   to set the boundary conditions, either of the end sections is
c   assumed to be connected to a uniform waveguide of the same cross
c   section as that of the end. the left end is assumed to be above
c   cutoff (propagating wave only). the right end can be above or below 
c   cutoff.
c
c   a single te(m,n) mode is assumed throughout the structure.
c   the program looks for a quasi steady-state solution(i.e. everything
c   repeats itself over a finite period). this is done by averaging the 
c   electron performance over a wave period.
c
c   note: 1. the equations and numerical methods are documented in:
c            k.r. chu, h.y. chen, c.l. hung, t.t. chang, l.r. barnett,
c            s.h. chen, t.t. yang, and d.j. dialetis, "theory and
c            experiment of ultrahigh-gain gyrotron traveling wave
c            amplifier", ieee trans. plasma science, vol. 27, pp. 391-
c            404 (1999).
c         2. comment statements in the main program beginning with ccc 
c            call for action by the user.
c
c   date of first version: january, 1993
c   date of this version: november 2, 2000
c   ************************************************************************
c
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
c   all variables beginning with a,b,d-h,j-z are real numbers. all
c   variables beginning with i are integers. all variables beginning 
c   with c are complex numbers.
      dimension axmn(9,8),cg(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(20001),arw(20001),anrho(20001),afamp(20001),
     &afphse(20001),apnet(20001),apfwd(20001),apbwd(20001),
     &apohmz(20001),agfwd(20001),agfwdz(20001),
     &agbwd(20001),agohm(20001)
      dimension atz(101),afreq(101),apin1(101,101),apin2(101,101),
     &apfwd1(101,101),apfwd2(101,101),aefwd1(101,101),aefwd2(101,101),
     &aebwd2(101,101),aeohm2(101,101),agbwd0(101,101),agbwd2(101,101),
     &agfwd0(101,101),agfwd1(101,101),agfwd2(101,101),
     &icont2(101,101),avalu2(101,101),apcnsv(101,101)
c   arc(it,ipz,irc,ipsi,iphi)
      dimension arc(11,41,1,11,11),apsi(11,41,1,11,11),
     &app(11,41,1,11,11),aphi(11,41,1,11,11),apz(11,41,1,11,11),
     &atime(11,41,1,11,11),agamma(11,41,1,11,11),
     &weight(11,41,1,11,11)
      dimension apti(11,41,1,11,11),aptf(11,41,1,11,11),
     &aptcsv(11,41,1,11,11),agmcsv(11,41,1,11,11)
      dimension itot(100)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,is,im
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      common/bext/bz0
      common/diagf/az,arw,anrho,afamp,afphse,apnet,apfwd,apbwd,apohmz,
     &pohmsm,idiag,icont
      common/diagb/arc,apsi,app,aphi,apz,atime,agamma
      common/phase/phsin,phsout,delphs
      common/check/agmcsv,apti,aptf,aptcsv,gmcsv,ptcsv
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
c   universal constants
      ci=dcmplx(0.0d0,1.0d0)
      pi=3.1415926d0
      realc=2.99792d10
      e=4.8032d-10
      me=9.1095d-28

ccc instruction for plotting within do loop 100 (1:plot, 0:do not plot)
      iplot=1
c
ccc instruction for writing within do loop 100 (1:write, 0:do not write)
      iwrite=1
c     
ccc specify procedure of calculation 
      itera=1
c   itera=1: iterative (calculate reflection coefficient at left end by
c            imposing outgoing wave boundary condition at right end)
c   itera=0: noniterative (set reflection coefficient at left end to 0 
c            and ignore boundary condition at right end)

ccc specify no of points on the z-axis to be used to mark the positions
c   of the left/right ends and the junctions between sections. 
      izmark=8
c

c   rf structure dimension arrays zmark rwl,and rwr specified below are
c   to be passed to the subprograms through the common block. the array 
c   elements must be in unit of cm.
c
c   zmark is a z-coordinate array to mark, from left to right, the posi- 
c   tions of the left end (zmark(1)), the junctions between sections 
c   (zmark(2),...,zmark(izmark-1)), and the right end (zmark(izmark)).
c
c   rwl(i) is the wall radius immediately to the left of z=zmark(i). 
c   rwr(i) is the wall radius immediately to the right of z=zmark(i). 
c   wall radius between zmark(i) and zmark(i+1) will be linearly
c   interpolated between rwr(i) and rwl(i+1) by function radius(z).
c
c   in theory, length of a uniform and cutoff end section does not affect
c   the results. in practice, set the length of the cutoff end section
c   (if any) sufficiently short to avoid numerical difficulties due to the 
c   exponential growth or attenuation of f with z.
c
ccc specify reference wall radius (rwref) in cm.
c   (for the calculation of cutoff freq., grazing mag. field, etc.)
      rwref=0.2659
ccc specify rf structure dimension arrays zmark, rwl, and rwr in cm. 
      r1=0.3175
      r2=rwref
      l1=2.0
      lt1=2.27
      l2=7.73d0
      lt2=2.0d0
      l3=7.73d0
      write(*,2) r1,r2,l1,lt1,l2,lt2,l3
    2 format( /' rf structure dimensions in cm:'/
     &' r1=',1pd11.4,' cm, r2=',1pd11.4,' cm, l1=',1pd11.4,
     &' cm, lt1=',1pd11.4,' cm'/' l2=',1pd11.4,' cm, lt2=',
     &1pd11.4,' cm, l3=',1pd11.4,' cm')

      zmark(1)=0.0
      zmark(2)=zmark(1)+l1
      zmark(3)=zmark(2)+lt1
      zmark(4)=zmark(3)+l2
      zmark(5)=zmark(4)+lt2
      zmark(6)=zmark(5)+l3
      zmark(7)=zmark(6)+lt1
      zmark(8)=zmark(7)+l1
      rwl(1)=r1
      rwr(1)=rwl(1)
      rwl(2)=r1
      rwr(2)=rwl(2)
      rwl(3)=r2
      rwr(3)=rwl(3)
      rwl(4)=r2
      rwr(4)=rwl(4)
      rwl(5)=r2
      rwr(5)=rwl(5)
      rwl(6)=r2
      rwr(6)=rwl(6)
      rwl(7)=r1
      rwr(7)=rwl(7)
      rwl(8)=r1
      rwr(8)=rwl(8)

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
      rhocu=1.72e-8
      rho2=1.5d4*rhocu
      rhol(1)=rhocu
      rhor(1)=rhocu
      rhol(2)=rhocu
      rhor(2)=rhocu
      rhol(3)=rhocu
      rhor(3)=rho2
      rhol(4)=rho2
      rhor(4)=rho2
      rhol(5)=rhocu
      rhor(5)=rhocu
      rhol(6)=rhocu
      rhor(6)=rhocu
      rhol(7)=rhocu
      rhor(7)=rhocu
      rhol(8)=rhocu
      rhor(8)=rhocu
c   write rf structure dimensions and resistivity
      write(*,3) izmark
    3 format(/' rf structure dimension arrays: (izmark=',i4,')')
      do 50 i=1,izmark,6
      if(izmark.ge.(i+5)) imax=i+5
      if(izmark.lt.(i+5)) imax=izmark
      write(*,4) (zmark(ii),ii=i,imax)
    4 format(/'   zmark(cm)=',6(1pd10.3,','))
      write(*,5) (rwl(ii),ii=i,imax)
    5 format('     rwl(cm)=',6(1pd10.3,','))
      write(*,6) (rwr(ii),ii=i,imax)
    6 format('     rwr(cm)=',6(1pd10.3,','))
   50 continue
      write(*,7) izmark
    7 format(/' wall resistivity arrays: (izmark=',i4,')')
      do 60 i=1,izmark,6
      if(izmark.ge.(i+5)) imax=i+5
      if(izmark.lt.(i+5)) imax=izmark
      write(*,8) (zmark(ii),ii=i,imax)
    8 format(/'   zmark(cm)=',6(1pd10.3,','))
      write(*,9) (rhol(ii),ii=i,imax)
    9 format(' rhol(ohm-m)=',6(1pd10.3,','))
      write(*,11) (rhor(ii),ii=i,imax)
   11 format(' rhor(ohm-m)=',6(1pd10.3,','))
   60 continue
      l=zmark(izmark)-zmark(1)
      delz=l/dfloat(1000)
      az(1)=zmark(1)
      arw(1)=radius(az(1))
      anrho(1)=rho(az(1))/1.72d-8
      rwmax=arw(1)
      nrhomx=anrho(1)
      do 20 i=1,1000
      az(i+1)=zmark(1)+delz*dfloat(i)
      arw(i+1)=radius(az(i+1))
      anrho(i+1)=rho(az(i+1))/1.72d-8
      if(arw(i+1).gt.rwmax) rwmax=arw(i+1)
      if(anrho(i+1).gt.nrhomx) nrhomx=anrho(i+1)
   20 continue
      if(iplot.eq.1) go to 101
      call sscale(zmark(1),zmark(izmark),0.0d0,rwmax)
      call splot(az,arw,1001,'rr  rw(cm) vs z(cm)$')
      if(nrhomx.eq.0.0d0) go to 102
      call sscale(zmark(1),zmark(izmark),0.0d0,nrhomx)
      call splot(az,anrho,1001,'**  rho(ohm-m)/1.72e-8 vs z(cm)$')
  102 continue
  101 continue

ccc specify m and n of te(m,n) mode
      im=1
      in=1

      xmn=axmn(iabs(im)+1,in)
      kcapmn=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
ccc instruction for cold test of the gyro-twt circuit
c   (icold=1: yes, icold=0: no)
      icold=1
      fhzmin=33.2d9
      fhzmax=38.0d9
      iwmax=81
      if(icold.eq.1) call cldtst(fhzmin,fhzmax,iwmax)
ccc specify cyclotron harmonic number
      is=1
c
ccc specify normalized beam electron guiding center position (nrcavg) and
c   normalized guiding center spread (ndrc). both are normalized to rwref.
c   (rc can not be set to 0 because that would make psi meaningless)
      do 1000 ircavg=1,1
      nrcavg=0.35
      ndrc=0.0
c
ccc specify beam voltage in kv
c   (all electrons are assumed to have the same initial energy)
      do 1000 ivbkv=1,1
      vbkv=99.45
c
ccc specify beam current in amp
      do 1000 irib=1,1
      ribamp=3.0
c
ccc specify beam alpha0
c   (alpha0 can not be set to 0 because that would make phi meaningless)
      do 1000 ialpha=1,1
      alpha0=0.85
c
c   calculate other parameters
      rcavg=nrcavg*rwref
      drc=ndrc*rwref
      nrc1=nrcavg-0.5*ndrc
      nrc2=nrcavg+0.5*ndrc
      rc1=nrc1*rwref
      rc2=nrc2*rwref
      opnrc1=oprc(is,im,in,1)
      opnrc2=oprc(is,im,in,2)
      opnrc3=oprc(is,im,in,3)
      gamma0=1.0d0+vbkv/511.0d0
      theta0=datan2(alpha0,1.0d0)
      angle=datan2(alpha0,1.0d0)*180.0/pi
      p0=dsqrt(gamma0**2-1.0d0)*me*realc
      pp0=p0*dsin(theta0)
      pz0=p0*dcos(theta0)
      v0=p0/(gamma0*me)
      vp0=v0*dsin(theta0)
      vz0=v0*dcos(theta0)
      gamz0=1.0/dsqrt(1.0-(vz0/realc)**2)
      rib=ribamp*3d9
      pbkw=vbkv*ribamp
      wcref=xmn*realc/rwref
      fcref=wcref/(2.0*pi)
      oeg=gamma0*wcref/(gamz0*dfloat(is))
      bzg=oeg*me*realc/e
c   calculate skin depth for rhomax at cutoff frequency (w.r.t. rwref)
      rhomax=nrhomx*1.72d-8
      delta=dsqrt(realc**2*rhomax/(2.0d0*pi*wcref*9.0d9))
c
ccc specify external magnetic field (bz0) in gauss or in ratio of bz0 to bzg 
      do 1000 ib=1,1
c     bz0=1.25525e4
c     bz0=1.25400e4
c     bz0bzg=bz0/bzg
      bz0bzg=0.99
      bz0=bz0bzg*bzg

      oe=bz0bzg*oeg
      oc=oe/gamma0
      fcyc=oc/(2.0*pi)
      rl0=vp0/oc
      nrl0=rl0/rwref
c   calculate points of intersections between waveguide mode (radius=rwref)
c   and beam-wave resonance line.
      if(bz0bzg.gt.1.0) go to 103
      fg=gamz0*fcref
      kzg=gamz0*wcref*vz0/realc**2
      f1=fg
      f2=fg
      kz1=kzg
      kz2=kzg
      go to 104
  103 continue
      dum1=1.0-(wcref/(dfloat(is)*oc*gamz0))**2
      f1=dfloat(is)*oc*gamz0**2*(1.0-(vz0/realc)*dsqrt(dum1))/(2.0*pi)
      f2=dfloat(is)*oc*gamz0**2*(1.0+(vz0/realc)*dsqrt(dum1))/(2.0*pi)
      kz1=dfloat(is)*(oc/realc)*gamz0**2*((vz0/realc)-dsqrt(dum1))
      kz2=dfloat(is)*(oc/realc)*gamz0**2*((vz0/realc)+dsqrt(dum1))
  104 continue
c     
      do 10 i1=1,101
      do 10 i2=1,101
      apfwd2(i1,i2)=-1.0d50
   10 continue
      pinmin=1.0d50
      pinmax=-1.0d50
ccc specify beam axial velocity spread tz (=delpz/pz0)
c   (always define itzmax, start itz from 1, and end itz at itzmax)
c   (to search for the saturation pin for the maximum saturated output 
c    power with a beam velocity spread, it is much faster to search for
c    the saturation pin for the cold beam case first.)
c   (this program will take the saturation pin obtained for case "itz"
c    as the starting pin for the search of the saturation pin for case
c    "itz+1".)
      itzmax=1
      do 200 itz=1,itzmax
      tz=0.05d0*(itz-1)
      atz(itz)=tz
c
ccc specify no of electrons
c
      itmax=7
      ipzmax=4*idint(tz*100.0d0)+1
      ircmax=1
      ipsimx=1
      iphimx=9
      if(im.eq.0) ipsimx=1
      if(tz.eq.0.0d0) ipzmax=1
      if(drc.eq.0.0d0) ircmax=1
      itotal=itmax*ipzmax*ircmax*ipsimx*iphimx
      itot(itz)=itotal
c
c   calculate maximum and minimum pp, pz, etc.
c
      delpz=tz*pz0
      dz=2.5
      pzmin=dmax1(0.01d0*delpz,pz0-dz*delpz)
      pzmax=dmin1(0.99d0*p0,pz0+dz*delpz)
      ppmax=dsqrt((p0/(me*realc))**2-(pzmin/(me*realc))**2)*me*realc
      ppmin=dsqrt((p0/(me*realc))**2-(pzmax/(me*realc))**2)*me*realc
      vpmin=ppmin/(gamma0*me)
      vpmax=ppmax/(gamma0*me)
      vzmin=pzmin/(gamma0*me)
      vzmax=pzmax/(gamma0*me)
      write(*,12) im,in,is,xmn,kcapmn
   12 format(/' mode: im=',i3,', in=',i3,', is=',i3,
     &', xmn=',1pd11.4,', kcapmn=',1pd11.4/)
      write(*,13)vbkv,ribamp,pbkw,alpha0,angle,gamma0,gamz0
   13 format(' vb=',1pd12.5,' kv, ib=',1pd10.3,' amp, pb=',
     &1pd10.3,' kw',/' alpha0=',1pd10.3,'(',1pd11.4,
     &' deg), gamma0=',1pd13.6,', gamz0=',1pd10.3/)
      write(*,14) rwref,nrcavg,ndrc,nrc1,nrc2,nrl0,
     &opnrc1,opnrc2,opnrc3,nrhomx,fcref,delta
   14 format(' rwref=reference wall radius=',1pd12.5,' cm'/
     &' rcavg/rwref=',1pd10.3,', drc/rwref=',1pd10.3/
     &' rc1/rwref=',1pd10.3,', rc2/rwref=',1pd10.3/
     &' rl0/rwref=',1pd10.3/
     &' nrc(1st optimum)=',1pd10.3/
     &' nrc(2nd optimum)=',1pd10.3/
     &' nrc(3rd optimum)=',1pd10.3//
     &' skin depth (for rho/rhocu=',1pd9.2,' and f=',1pd9.2,
     &' hz) =',1pd10.3,' cm'/)
      write(*,15) bz0bzg,bz0,bzg
   15 format(' bz0/bzg=',1pd12.5,', bz0=',1pd12.5,' gauss'/
     &' bzg (grazing magnetic field w.r.t. rwref)=',1pd12.5,' gauss')
      write(*,16) fcref,fcyc,f1,f2,kz1,kz2
   16 format(' cutoff freq. (w.r.t. rwref)=',1pd12.5,' hz,',
     &' cyc. freq.=',1pd12.5,' hz'/
     &' points of intersections in freq. vs kz diagram (w.r.t. rwref):'
     &/'   (f1, f2)=',1pd11.4,',',1pd11.4,' hz'/' (kz1, kz2)=',1pd11.4,
     &',',1pd11.4,' /cm'/)
      write(*,17) tz,p0,pp0,pz0,delpz,ppmin,ppmax,pzmin,pzmax,
     &v0,vp0,vz0,vpmin,vpmax,vzmin,vzmax
   17 format(' momentum and velocity are in gaussian units:'/
     &' tz=delpz/pz0=',1pd10.3/' p0=',1pd10.3,
     &' pp0=',1pd10.3,', pz0=',1pd10.3,', delpz=',1pd10.3/
     &' ppmin=',1pd10.3,', ppmax=',1pd10.3,/
     &' pzmin=',1pd10.3,', pzmax=',1pd10.3/
     &' v0=',1pd10.3,', vp0=',1pd10.3,', vz0=',1pd10.3/
     &' vpmin=',1pd10.3,', vpmax=',1pd10.3/
     &' vzmin=',1pd10.3,', vzmax=',1pd10.3/)
      write(*,19) itera,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
   19 format(' itera=',i3/' itmax=',i3,', ipzmax=',i3,', ircmax=',i3,
     &', ipsimx=',i3,', iphimx=',i3,', itotal=',i7/)

      do 30 i1=1,101
      do 30 i2=1,101
      apfwd1(i1,i2)=0.0
      agfwd1(i1,i2)=0.0
      aefwd1(i1,i2)=0.0
   30 continue

ccc specify incident wave frequency.  
c   (always define iwmax, start iw from 1, and end iw at iwmax)
      iwmax=13
      iwmax=1
      do 100 iw=1,iwmax
      fhzmin=33.2d9
      fhzmax=38.0d9
      if(iwmax.gt.1) fstep=(fhzmax-fhzmin)/dfloat(iwmax-1)
      if(iwmax.eq.1) fstep=0.0
      fhz=fhzmin+fstep*dfloat(iw-1)
      if(iw.eq.1) fhz=33.2d9
      if(iw.eq.2) fhz=33.4d9
      if(iw.eq.3) fhz=33.6d9
      if(iw.eq.4) fhz=33.8d9
      if(iw.eq.5) fhz=34.0d9
      if(iw.eq.6) fhz=34.1d9
      if(iw.eq.7) fhz=34.2d9
      if(iw.eq.8) fhz=34.3d9
      if(iw.eq.9) fhz=34.5d9
      if(iw.eq.10) fhz=35.0d9
      if(iw.eq.11) fhz=36.0d9
      if(iw.eq.12) fhz=37.0d9
      if(iw.eq.13) fhz=38.0d9
      if(iwmax.eq.1) fhz=34.2d9
      w=2.0*pi*fhz
      afreq(iw)=fhz
c   angular frequency of the incident wave (in radian/sec, denoted by w) 
c   is passed to subprograms through the common block.

      iskip=0
ccc specify input (forward) wave power (in watts) and phase (named phsin,
c   in radian) at left end. phase of output wave (named phsout) and phase
c   difference between output and input waves (named delphs) will be
c   returned to the main program by function cbc through the common block.
c      note:
c   1. always define ipinmx, start ipin from 1, and end ipin at ipinmx.
c   2. if there is difficulty in obtaining a solution with itera=1, try
c      to start pin in the small signal regime where it is numerically 
c      easier to obtain a solution by the iteration procedure without a
c      good guess for the reflection coefficient.
c   3. reflection coefficient obtained for the present pin will be supplied
c      to muller as the guessed reflection coefficient for the next pin.
c
      iphsmx=8
      iphsmx=1
      do 100 iphase=1,iphsmx
      phsin=0.0d0+(iphase-1)*2.0d0*pi/dfloat(iphsmx)
c
      ipinmx=31
      ipinmx=1
      do 100 ipin=1,ipinmx
      pin=8.0*1.30**(ipin-2)
      if(itz.gt.1) pin=apin2(iw,itz-1)*1.3**(ipin-4)
      if(ipin.eq.1) pin=1.0d-3
      if(ipinmx.eq.1.and.itz.eq.1) pin=20.0d0
      if(ipinmx.eq.1.and.itz.eq.2) pin=30.0d0

      apin1(ipin,iw)=pin
c   the following statement will terminate do loop ipin when the calculated
c   output power is lower than the maximum power by a certain amount (specified
c   near the end of do loop 100). cancellation of the following statement
c   will allow completion of do loop ipin.
      if(iskip.eq.1) go to 100
c
      if(pin.gt.pinmax) pinmax=pin
      if(pin.lt.pinmin) pinmin=pin
      if(ipinmx.eq.1) pinmax=2.0d0*pinmin

c   calculate no of cyclotron orbits (ncyc) and no of guide wavelengths 
c   (nwave) over the total length of the rf structure.
      ncyc=fcyc*(zmark(izmark)-zmark(1))/vz0  
      kz2=(w/realc)**2-xmn**2/rwref**2
      kz=dsqrt(kz2)
      lamdag=2.0d0*pi/kz
      nwave=(zmark(izmark)-zmark(1))/lamdag
      nlarge=dmax1(ncyc,nwave)

ccc specify no of steps for z-integration. 
      do 40 idiv=6,6
      izstep=idiv*idint(nlarge+1.0d0)
c
c   for a new problem, always check convergence with respect to izstep. 
c   do loop 40 is for such a test. "idiv" is the no of steps per cyclotron 
c   orbit or per wavelength (2*pi/kz), whichever makes izstep greater.

ccc guess complex reflection coefficient  
      cguess=dcmplx(0.5d0,0.5d0)
      if(ipin.eq.1) cg(1)=cguess
c   the complex reflection coefficient (denoted by variable cg(1)) has been
c   formulated to be the root of the equation cbc(cg)=0. cg(1) defined above 
c   is a guessed value to be supplied to subroutine muller as a starting 
c   value for the search of the correct reflection coefficient (referred to 
c   the z=0 plane).

c   calculate reflection coefficient (referred to the z=0 plane)
      ep1=1.0d-4
      ep2=ep1
      iroot=1
      imaxit=50
      icont=0
      idiag=0
      if(itera.eq.1) call muller(0,iroot,cg,imaxit,ep1,ep2,cbc,.false.)
      if(itera.ne.1) cg(1)=dcmplx(0.0d0,0.0d0)
      idiag=1
      value=cdabs(cbc(cg(1)))
c**************************************
c  see comments in subroutine pdiag for the the diagnostic functions of
c  this subroutine
c     call pdiag(21)
c**************************************
c   muller is called to solve equation cbc(cg)=0 for the complex root cg.
c   'idiag' is an instruction for diagnostics. if (and only if) idiag=1,
c   function cbc will diagnose and store the profiles of wall radius, resis-
c   tivity, rf field, wave powers, ohmic losses, gain, final beam coordinates,
c   and conservation/consistency factors, and pass them through the common 
c   block to the main program.
c   'icont' monitors the no of times that function cbc is called for root 
c   searching ('icont' goes up by one opon each call to cbc).
c   'value' monitors the accuracy of the root returned by muller.
      gr=cg(1)
      gi=-ci*cg(1)
      gamp=cdabs(cg(1))
      if(gr.ne.0.0d0.or.gi.ne.0.0d0) gphase=datan2(gi,gr)
      if(gr.eq.0.0d0.and.gi.eq.0.0d0) gphase=0.0
c   diagnose gain, efficiencies, and power conservation factor, etc.
      etab=0.0
      do 90 it=1,itmax
      do 90 ipz=1,ipzmax
      do 90 irc=1,ircmax
      do 90 ipsi=1,ipsimx
      do 90 iphi=1,iphimx
      gamma=agamma(it,ipz,irc,ipsi,iphi)
      etab=etab+weight(it,ipz,irc,ipsi,iphi)*(gamma0-gamma)
     &/(gamma0-1.0d0)
   90 continue
      etafwd=(apfwd(izstep+1)-pin)/(pbkw*1.0d3)
      etabwd=(apbwd(1)-apbwd(izstep+1))/(pbkw*1.0d3)
      etaohm=pohmsm/(pbkw*1.0d3)
      ratio1=pohmsm/apfwd(izstep+1)
      ratio2=apbwd(izstep+1)/apfwd(izstep+1)
      pconsv=(apnet(izstep+1)-apnet(1)+pohmsm-etab*pbkw*1.0d3)
     &/(etab*pbkw*1.0d3)
      gfwd=10.0d0*dlog10(apfwd(izstep+1)/pin)
      if(pohmsm.gt.0.0d0) pohmdb=10.0*dlog10(pohmsm/pin)
      delz=(zmark(izmark)-zmark(1))/dfloat(izstep)
      agfwd(1)=0.0
      agfwdz(1)=0.0
      if(apbwd(1).gt.0.0d0) agbwd(1)=10.0d0*dlog10(apbwd(1)/pin)
      if(apbwd(1).eq.0.0d0) agbwd(1)=-1.0d50
      if(apohmz(1).gt.0.0d0) agohm(1)=10.0d0*dlog10(apohmz(1)/pin)
      if(apohmz(1).eq.0.0d0) agohm(1)=-1.0d50
      do 80 iz=1,izstep
      agfwd(iz+1)=10.0d0*dlog10(apfwd(iz+1)/pin)
      agfwdz(iz+1)=(agfwd(iz+1)-agfwd(iz))/delz
      agbwd(iz+1)=10.0d0*dlog10(apbwd(iz+1)/pin)
      if(apohmz(iz+1).gt.0.0d0)
     &agohm(iz+1)=10.0d0*dlog10(apohmz(iz+1)/pin)
      if(apohmz(iz+1).eq.0.0d0) agohm(iz+1)=-1.0d50
   80 continue

      if(iplot.ne.1) go to 105
c   plot rf structure shape, wall resistivity, rf field and power profiles.
c
c   find rwmax, nrhomx, fampmx, pnetmx, pnetmn, pfwdmx, pfwdmn, gbwdmx, gbwdmn, 
c   gohmmx, gohmmn, gfwdmx, gfwdmn, gzmax, and gzmin
      rwmax=arw(1)
      nrhomx=anrho(1)
      fampmx=afamp(1)
      gbwdmx=agbwd(1)
      gbwdmn=agbwd(1)
      pnetmx=apnet(1)
      pnetmn=apnet(1)
      gohmmx=agohm(1)
      gohmmn=agohm(1)
      gfwdmx=agfwd(1)
      gfwdmn=agfwd(1)
      gzmax=agfwdz(1)
      gzmin=agfwdz(1)
      do 70 iz=1,izstep
      if(arw(iz+1).gt.rwmax) rwmax=arw(iz+1)
      if(anrho(iz+1).gt.nrhomx) nrhomx=anrho(iz+1)
      if(afamp(iz+1).gt.fampmx) fampmx=afamp(iz+1)
      if(agbwd(iz+1).gt.gbwdmx) gbwdmx=agbwd(iz+1)
      if(agbwd(iz+1).le.gbwdmn) gbwdmn=agbwd(iz+1)
      if(apnet(iz+1).gt.pnetmx) pnetmx=apnet(iz+1)
      if(apnet(iz+1).le.pnetmn) pnetmn=apnet(iz+1)
      if(agohm(iz+1).gt.gohmmx) gohmmx=agohm(iz+1)
      if(agohm(iz+1).le.gohmmn) gohmmn=agohm(iz+1)
      if(agfwd(iz+1).gt.gfwdmx) gfwdmx=agfwd(iz+1)
      if(agfwd(iz+1).le.gfwdmn) gfwdmn=agfwd(iz+1)
      if(dabs(agfwdz(iz+1)).gt.10.0d0) go to 70
      if(agfwdz(iz+1).gt.gzmax) gzmax=agfwdz(iz+1)
      if(agfwdz(iz+1).le.gzmin) gzmin=agfwdz(iz+1)
   70 continue
      if(pnetmn.ge.0.0d0) pnetmn=0.0d0
c     if(gfwdmn.ge.0.0d0) gfwdmn=0.0d0
c     if(gzmin.ge.0.0d0) gzmin=0.0d0
c     if(gbwdmn.ge.0.0d0) gbwdmn=0.0d0
c     if(gohmmn.ge.0.0d0) gohmmn=0.0d0
      if(gbwdmn.lt.(gbwdmx-30.0d0)) gbwdmn=gbwdmx-30.0d0
      if(gohmmn.lt.(gohmmx-30.0d0)) gohmmn=gohmmx-30.0d0

      call sscale(zmark(1),zmark(izmark),0.0d0,rwmax)
      call splot(az,arw,izstep+1,'rr  rw(cm) vs z(cm)$')
      if(nrhomx.eq.0.0d0) go to 106
      call sscale(zmark(1),zmark(izmark),0.0d0,nrhomx)
      call splot(az,anrho,izstep+1,'**  rho(ohm-m)/1.72e-8 vs z(cm)$')
  106 continue
      call sscale(zmark(1),zmark(izmark),0.0d0,fampmx)
      call splot(az,afamp,izstep+1,'ff  amplitude of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),-pi,pi)
      call splot(az,afphse,izstep+1,'pp  phase(radian) of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),pnetmn,pnetmx)
      call splot(az,apnet,izstep+1,'nn  net wave power(w) vs z(cm)$') 
      write(*,21)
   21 format(/ ' note: forward wave power (hence the plot below) may be
     &inaccurately diagnosed'/'       in regions where w is slightly abo
     &ve cutoff or where the wall is tapered.'/
     &'       it was artificially set to 1.0d-10 watt where w < wcmn.')
      call sscale(zmark(1),zmark(izmark),gfwdmn,gfwdmx)
      call splot(az,agfwd,izstep+1,
     &'gg  gain=10*log(forward wave power/input power) vs z(cm)$') 
      write(*,21)
      call sscale(zmark(1),zmark(izmark),gzmin,gzmax)
      call splot(az,agfwdz,izstep+1,
     &'gg  forward wave growth rate(db/cm) vs z(cm)$')
      if(pohmsm.eq.0.0d0) go to 107
      call sscale(zmark(1),zmark(izmark),gohmmn,gohmmx)
      call splot(az,agohm,izstep+1,
     &'oo  10*log(ohmic loss per cm/input power) vs z(cm)$')
  107 continue
      write(*,22)
   22 format(/ ' note: backward wave power (hence the plot below) may be
     & inaccurately diagnosed'/'       in regions where w is slightly ab
     &ove cutoff or where the wall is tapered.'/
     &'       it was artificially set to 1.0d-10 watt where w < wcmn.')
      call sscale(zmark(1),zmark(izmark),gbwdmn,gbwdmx)
      call splot(az,agbwd,izstep+1, 
     &'bb  10*log(backward wave power/input power) vs z(cm)$') 
      write(*,23)
   23 format(/)
  105 continue

      if(iwrite.ne.1) go to 108
c   write results
      write(*,26) fhz,pin,itotal,itera
   26 format(' freq=',1pd12.5,' hz, pin(incident power)=',1pd10.3,
     &' w, itotal=',i7,', itera=',i1)
      write(*,25) ncyc,nwave,lamdag,idiv,izstep
   25 format(' ncyc=',1pd9.2,', nwave=',1pd9.2,', lamdag=',1pd9.2,
     &' cm, idiv=',i4,', izstep=',i7)
      write(*,24) gamp,gphase,icont,value,ep1
   24 format(' gamp=',1pd11.4,', gphase(radian)=',1pd10.3,
     &'(',i3,',',1pd8.1,'), ep1=',1pd9.2)
      write(*,97) phsin,phsout,delphs
   97 format(' phsin(radian)=',1pd10.3,', phsout(radian)=',1pd10.3,
     &', delphs(radian)=',1pd10.3)
      write(*,27) apnet(1),apnet(izstep+1),apfwd(1),
     &apfwd(izstep+1),apbwd(1),apbwd(izstep+1),pohmsm,pconsv
   27 format(' pnet (net wave power):     at left end=',1pd10.3,
     &' w, at right end=',1pd10.3,' w'/ 
     &' pfwd (forward wave power): at left end=',1pd10.3,
     &' w, at right end=',1pd10.3,' w'/
     &' pbwd (backward wave power):at left end=',1pd10.3,
     &' w, at right end=',1pd10.3,' w!'/
     &' pohmsm (sum of wall losses)=',1pd10.3,
     &' w,        (power conservation=',1pd9.2,')')
      write(*,28) gfwd,gmcsv
   28 format(' gfwd=10*log(pfwd at right end/pin)= ',1pd10.3,' db,',
     &'(gamma consistency=',1pd9.2,')')
      if(im.eq.0) write(*,29) ptcsv
   29 format(49x,'(ptheta conservation=',1pd9.2,')')
      write(*,34) agbwd(1)
   34 format(' gbwd=10*log(pbwd at left end/pin)=',1pd10.3,' db')
      write(*,35) ratio1,ratio2
   35 format(' pohmsm/pfwd(right end)=',1pd10.3,
     &', pbwd(right end)/pfwd(right end)=',1pd10.3,'!')
      write(*,36) etafwd,etab,etabwd,etaohm
   36 format(' etafwd=',1pd11.4,', etab=',1pd11.4,
     &'  (etabwd=',1pd11.4,', etaohm=',1pd11.4,')')
      write(*,37)
   37 format('**********************************************************
     &**********************')
  108 continue
   40 continue
c   for a fixed tz (beam velocity spread), store pfwd (forward wave power at
c   right end), gfwd (=10*log10(pfwd/pin)), etafwd (=(pfwd-pin)/beam power)
c   as functions of the specified pin (input power) for different frequencies
c   as specified. 
c   store gfwd at the minimum specified pin as a function of frequency for
c   different tz as specified.
c   find and store the maximized pfwd (within the specified pin range) at the
c   output end and the associated pin, gfwd, etafwd, etabwd (=(pbwd at left end
c   -pbwd at right end)/beam power), and etaohm (total ohmic power/beam power)
c   as functions of the specified frequency for different tz as specified.
      apfwd1(ipin,iw)=apfwd(izstep+1)
      agfwd1(ipin,iw)=agfwd(izstep+1)
      aefwd1(ipin,iw)=etafwd
      if(ipin.eq.1)
     &agfwd0(iw,itz)=10.0*dlog10(apfwd(izstep+1)/pin)
      if(ipin.eq.1)
     &agbwd0(iw,itz)=agbwd(1)
      if(apfwd1(ipin,iw).lt.(0.9d0*apfwd2(iw,itz))) iskip=1
      if(apfwd1(ipin,iw).lt.apfwd2(iw,itz))  go to 100
      apin2(iw,itz)=pin
      agfwd2(iw,itz)=agfwd1(ipin,iw)
      aefwd2(iw,itz)=aefwd1(ipin,iw)
      aebwd2(iw,itz)=etabwd
      aeohm2(iw,itz)=etaohm
      agbwd2(iw,itz)=agbwd(1)
      icont2(iw,itz)=icont
      avalu2(iw,itz)=value
      apcnsv(iw,itz)=pconsv
      apfwd2(iw,itz)=apfwd1(ipin,iw)
  100 continue
      if(ipinmx.eq.1) go to 110
c   plot pfwd, gfwd, and efwd at the output end as functions of pin for
c   fixed tz and different frequencies as specified.
c
c   find pfwdmx, gfwdmx, gfwdmn, efwdmx, and efwdmn
      pfwdmx=-1.0d50
      pfwdmn=1.0d50
      gfwdmx=-1.0d50
      gfwdmn=1.0d50
      efwdmx=-1.0d50
      efwdmn=1.0d50
      do 57 iw=1,iwmax
      do 57 ipin=1,ipinmx
      if(apfwd1(ipin,iw).gt.pfwdmx) pfwdmx=apfwd1(ipin,iw)
      if(apfwd1(ipin,iw).le.pfwdmn) pfwdmn=apfwd1(ipin,iw)
      if(agfwd1(ipin,iw).gt.gfwdmx) gfwdmx=agfwd1(ipin,iw)
      if(agfwd1(ipin,iw).le.gfwdmn) gfwdmn=agfwd1(ipin,iw)
      if(aefwd1(ipin,iw).gt.efwdmx) efwdmx=aefwd1(ipin,iw)
      if(aefwd1(ipin,iw).le.efwdmn) efwdmn=aefwd1(ipin,iw)
   57 continue
      if(pfwdmn.eq.pfwdmx) pfwdmx=2.0d0*dabs(pfwdmn)
      if(gfwdmn.eq.gfwdmx) gfwdmx=2.0d0*dabs(gfwdmn)
      if(efwdmn.eq.efwdmx) efwdmx=2.0d0*dabs(efwdmn)
      call sscale(pinmin,pinmax,pfwdmn,pfwdmx)
      imax=iwmax-1
      if(iwmax.eq.1) go to 52
      do 51 i=1,imax
   51 call splot(apin1(1,i),apfwd1(1,i),ipinmx,'*   $')
   52 call splot(apin1(1,iwmax),apfwd1(1,iwmax),ipinmx,
     &'pp  pfwd(w) vs pin(w)$')
      call sscale(pinmin,pinmax,gfwdmn,gfwdmx)
      if(iwmax.eq.1) go to 54
      do 53 i=1,imax
   53 call splot(apin1(1,i),agfwd1(1,i),ipinmx,'*   $')
   54 call splot(apin1(1,iwmax),agfwd1(1,iwmax),ipinmx,
     &'gg  gfwd(db) vs pin(w)$')
      call sscale(pinmin,pinmax,efwdmn,efwdmx)
      if(iwmax.eq.1) go to 56
      do 55 i=1,imax
   55 call splot(apin1(1,i),aefwd1(1,i),ipinmx,'*   $')
   56 call splot(apin1(1,iwmax),aefwd1(1,iwmax),ipinmx,
     &'ee  etafwd vs pin(w)$')
      write(*,37)
      if(itz.eq.itzmax) go to 109
      write(*,38)
   38 format(/'       -- intermediate summary for do loop 200 --')
c   write pfwd at right end (maximized within the range of specified pin)
c   and the associated pin, gfwd, etafwd, etabwd, etaohm as functions of
c   the specified frequency for a fixed tz. numbers in parentheses are
c   gfwd for the minimum specified pin (which is the small signal gain if
c   pinmin is sufficiently small).
      write(*,39) atz(itz),itot(itz),itera,idiv,ep1
   39 format(/'tz=',1pd9.2,', itotal=',i7,', itera=',i2,
     &', idiv=',i4,', ep1=',1pd8.1/)
      write(*,41)
   41 format(' freq(hz):','  pfwd(w)','    etafwd ','  gfwd(db)',
     &' (at pinmin)  ',' gbwd(db)',' (at pinmin)')
      do 91 iw=1,iwmax
      write(*,42) afreq(iw),apfwd2(iw,itz),aefwd2(iw,itz),        
     &agfwd2(iw,itz),agfwd0(iw,itz),agbwd2(iw,itz),agbwd0(iw,itz)
   42 format(1pd10.4,':',1pd9.2,' ',1pd9.2,' ',1pd9.2,' (',1pd9.2,
     &')  ',1pd9.2,' (',1pd9.2,')')
   91 continue
      write(*,43)
   43 format(/' freq(hz):','   pin(w) ','   etabwd ','   etaohm',
     &'   icont','   value','     pconsv')
      do 92 iw=1,iwmax
      write(*,44) afreq(iw),apin2(iw,itz),aebwd2(iw,itz),aeohm2(iw,itz),
     &icont2(iw,itz),avalu2(iw,itz),apcnsv(iw,itz)
   44 format(1pd10.4,':',1pd9.2,' ',1pd9.2,' ',1pd9.2,
     &'  ',i4,'  ',1pd9.2,' ',1pd9.2)
   92 continue
      write(*,37)
  109 continue
  110 continue
  200 continue
      if(iwmax.eq.1) go to 111
c   plot gfwd for the minimum specified pin (which is the small signal
c   gain if pinmin is sufficiently amall) as a function of frequency for 
c   different tz as specified.
c   plot pfwd at the output end (maximized w.r.t. pin in the specified pin
c   range) and the associated pin, gfwd, etafwd, etabwd, and etaohm as
c   functions of frequency for different tz as specified.
c
c   find pfwdmx, gfwdmx, gfwdmn, efwdmx, efwdmn, ebwdmx, ebwdmn, gbwdmx,
c   gbwdmn, and eohmmx
      pfwdmx=-1.0d50
      pfwdmn=1.0d50
      gfwdmx=-1.0d50
      gfwdmn=1.0d50
      efwdmx=-1.0d50
      efwdmn=1.0d50
      ebwdmx=-1.0d50
      ebwdmn=1.0d50
      gbwdmx=-1.0d50
      gbwdmn=1.0d50
      eohmmx=-1.0d50
      eohmmn=1.0d50
      do 59 itz=1,itzmax
      do 59 iw=1,iwmax
      if(apfwd2(iw,itz).gt.pfwdmx) pfwdmx=apfwd2(iw,itz)
      if(apfwd2(iw,itz).le.pfwdmn) pfwdmn=apfwd2(iw,itz)
      if(agfwd0(iw,itz).gt.gfwdmx) gfwdmx=agfwd0(iw,itz)
      if(agfwd2(iw,itz).gt.gfwdmx) gfwdmx=agfwd2(iw,itz)
      if(agfwd0(iw,itz).le.gfwdmn) gfwdmn=agfwd0(iw,itz)
      if(agfwd2(iw,itz).le.gfwdmn) gfwdmn=agfwd2(iw,itz)
      if(aefwd2(iw,itz).gt.efwdmx) efwdmx=aefwd2(iw,itz)
      if(aefwd2(iw,itz).le.efwdmn) efwdmn=aefwd2(iw,itz)
      if(aebwd2(iw,itz).gt.ebwdmx) ebwdmx=aebwd2(iw,itz)
      if(aebwd2(iw,itz).le.ebwdmn) ebwdmn=aebwd2(iw,itz)
      if(agbwd2(iw,itz).gt.gbwdmx) gbwdmx=agbwd2(iw,itz)
      if(agbwd0(iw,itz).gt.gbwdmx) gbwdmx=agbwd0(iw,itz)
      if(agbwd2(iw,itz).le.gbwdmn) gbwdmn=agbwd2(iw,itz)
      if(agbwd0(iw,itz).le.gbwdmn) gbwdmn=agbwd0(iw,itz)
      if(aeohm2(iw,itz).gt.eohmmx) eohmmx=aeohm2(iw,itz)
      if(aeohm2(iw,itz).le.eohmmn) eohmmn=aeohm2(iw,itz)
   59 continue
      if(pfwdmn.ge.0.0) pfwdmn=0.0
c     if(gfwdmn.ge.0.0d0) gfwdmn=0.0d0
c     if(gbwdmn.ge.0.0d0) gbwdmn=0.0d0
c     if(efwdmn.ge.0.0d0) efwdmn=0.0d0
c     if(ebwdmn.ge.0.0d0) ebwdmn=0.0d0
c     if(eohmmn.ge.0.0d0) eohmmn=0.0d0
      imax=itzmax-1
      call sscale(afreq(1),afreq(iwmax),gfwdmn,gfwdmx)
      if(itzmax.eq.1) go to 62
      do 61 i=1,imax
   61 call splot(afreq,agfwd0(1,i),iwmax,'*   $')
   62 call splot(afreq,agfwd0(1,itzmax),iwmax,
     &'ss  gfwd(in db, at pin=pinmin) vs freq(hz)$')
      call sscale(afreq(1),afreq(iwmax),pfwdmn,pfwdmx)
      if(itzmax.eq.1) go to 64
      do 63 i=1,imax
   63 call splot(afreq,apfwd2(1,i),iwmax,'*   $')
   64 call splot(afreq,apfwd2(1,itzmax),iwmax,
     &'pp  pfwd(in w, maximized in specified range of pin) vs freq(hz)
     &$')
      call sscale(afreq(1),afreq(iwmax),pinmin,pinmax)
      if(itzmax.eq.1) go to 66
      do 65 i=1,imax
   65 call splot(afreq,apin2(1,i),iwmax,'*   $')
   66 call splot(afreq,apin2(1,itzmax),iwmax,
     &'dd  pin(in w, for above pfwd) vs freq(hz)$')
      call sscale(afreq(1),afreq(iwmax),gfwdmn,gfwdmx)
      if(itzmax.eq.1) go to 68
      do 67 i=1,imax
   67 call splot(afreq,agfwd2(1,i),iwmax,'*   $')
   68 call splot(afreq,agfwd2(1,itzmax),iwmax,
     &'gg  gfwd(in db, for above pfwd) vs freq(hz)$')
      call sscale(afreq(1),afreq(iwmax),efwdmn,efwdmx)
      if(itzmax.eq.1) go to 72
      do 71 i=1,imax
   71 call splot(afreq,aefwd2(1,i),iwmax,'*   $')
   72 call splot(afreq,aefwd2(1,itzmax),iwmax,
     &'ee  etafwd(for above pfwd) vs freq(hz)$')
      if(eohmmx.le.0.0d0) go to 112
      call sscale(afreq(1),afreq(iwmax),eohmmn,eohmmx)
      if(itzmax.eq.1) go to 74
      do 73 i=1,imax
   73 call splot(afreq,aeohm2(1,i),iwmax,'*   $')
   74 call splot(afreq,aeohm2(1,itzmax),iwmax,
     &'oo  etaohm(for above pfwd) vs freq(hz)$')
  112 continue
      call sscale(afreq(1),afreq(iwmax),ebwdmn,ebwdmx)
      if(itzmax.eq.1) go to 76
      do 75 i=1,imax
   75 call splot(afreq,aebwd2(1,i),iwmax,'*   $')
   76 call splot(afreq,aebwd2(1,itzmax),iwmax,
     &'bb  etabwd(for above pfwd) vs freq(hz)$')
      call sscale(afreq(1),afreq(iwmax),gbwdmn,gbwdmx)
      if(itzmax.eq.1) go to 78
      do 77 i=1,imax
   77 call splot(afreq,agbwd2(1,i),iwmax,'*   $')
   78 call splot(afreq,agbwd2(1,itzmax),iwmax,
     &'bb  10*log(pbwd at left end(for above pfwd)/pin) vs freq(hz)$')
      call sscale(afreq(1),afreq(iwmax),gbwdmn,gbwdmx)
      if(itzmax.eq.1) go to 82
      do 81 i=1,imax
   81 call splot(afreq,agbwd0(1,i),iwmax,'*   $')
   82 call splot(afreq,agbwd0(1,itzmax),iwmax,
     &'ss  10*log(pbwd at left end(at pin=pinmin)/pin) vs freq(hz)$')
  111 continue
c              -- final summary for do loop 200 --
      write(*,45)
   45 format(/'          -- final summary for do loop 200 --')
c   write pfwd at right end (maximized within the range of specified pin)
c   and the associated pin, gfwd, etafwd, etabwd, etaohm as functions of
c   the specified frequency for different tz. numbers in parentheses are
c   gfwd for the minimum specified pin (which is the small signal gain if
c   pinmin is sufficiently small).
      do 93 itz=1,itzmax
      write(*,39) atz(itz),itot(itz),itera,idiv,ep1
      write(*,41)
      do 94 iw=1,iwmax
      write(*,42) afreq(iw),apfwd2(iw,itz),aefwd2(iw,itz),        
     &agfwd2(iw,itz),agfwd0(iw,itz),agbwd2(iw,itz),agbwd0(iw,itz)
   94 continue
      write(*,43)
      do 95 iw=1,iwmax
      write(*,44) afreq(iw),apin2(iw,itz),aebwd2(iw,itz),aeohm2(iw,itz),
     &icont2(iw,itz),avalu2(iw,itz),apcnsv(iw,itz)
   95 continue
      write(*,37)
   93 continue
 1000 continue
      stop
      end
c ******************************************************************
      function cbc(cg)
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(500000),dy(500000),q(500000)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(20001),arw(20001),anrho(20001),afamp(20001),
     &afphse(20001),apnet(20001),apfwd(20001),apbwd(20001),
     &apohmz(20001)
      dimension arc(11,41,1,11,11),apsi(11,41,1,11,11),
     &app(11,41,1,11,11),aphi(11,41,1,11,11),apz(11,41,1,11,11),
     &atime(11,41,1,11,11),agamma(11,41,1,11,11),
     &weight(11,41,1,11,11)
      dimension apti(11,41,1,11,11),aptf(11,41,1,11,11),
     &aptcsv(11,41,1,11,11),agmcsv(11,41,1,11,11)
      dimension adphi(20001,301),apx(301,20001),
     &apy(301,20001),azz(20001)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,is,im
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      common/diagf/az,arw,anrho,afamp,afphse,apnet,apfwd,apbwd,apohmz,
     &pohmsm,idiag,icont
      common/diagb/arc,apsi,app,aphi,apz,atime,agamma
      common/phase/phsin,phsout,delphs
      common/check/agmcsv,apti,aptf,aptcsv,gmcsv,ptcsv
      common/diagp/adphi,apx,apy,azz,izmax
      external difeq
      z1=zmark(1)
      rw=radius(z1)
      if((w/realc)**2.le.xmn**2/rw**2) write(*,1)
      if((w/realc)**2.le.xmn**2/rw**2) stop 100
    1 format(' *** w below cutoff in first section ***')
c
      ie=7*itotal
      ieqlst=ie+6
c initialize coordinate part of coordinate/field array y at left end 
c according to the initial electron beam distribution in real and momentum 
c spaces.
      call coord(arc,apsi,app,aphi,apz,atime,agamma,weight)
      i=0
      do 40 it=1,itmax
      do 40 ipz=1,ipzmax
      do 40 irc=1,ircmax
      do 40 ipsi=1,ipsimx
      do 40 iphi=1,iphimx
      y(i+1)=arc(it,ipz,irc,ipsi,iphi)
      y(i+2)=apsi(it,ipz,irc,ipsi,iphi)
      y(i+3)=app(it,ipz,irc,ipsi,iphi)
      y(i+4)=aphi(it,ipz,irc,ipsi,iphi)
      y(i+5)=apz(it,ipz,irc,ipsi,iphi)
      y(i+6)=atime(it,ipz,irc,ipsi,iphi)
      y(i+7)=agamma(it,ipz,irc,ipsi,iphi)
      i=i+7
c *****************************************************************
      if(idiag.ne.1) go to 399
c  diagnose electron phase space distribution (for subroutine pdiag)
      azz(1)=z1
      if(it.ne.1) go to 399 
      if(ipz.ne.1) go to 399 
      if(irc.ne.1) go to 399 
      if(ipsi.ne.1) go to 399 
      phase=aphi(it,ipz,irc,ipsi,iphi)
      apx(iphi,1)=app(it,ipz,irc,ipsi,iphi)*dcos(phase)
      apy(iphi,1)=app(it,ipz,irc,ipsi,iphi)*dsin(phase)
  399 continue
c *****************************************************************
   40 continue

c initialize field part of coordinate/field array y at left end according
c to the input forward wave power (pin) specified in the main program, and
c store radius, resistivity, rf field, and powers at left end.
c (wave power in w, wall loss per unit length in w/cm)
      ckz2=(w/realc)**2-xmn**2*closs(z1)/rw**2
      ckz=cdsqrt(ckz2)
      cdum1=cdexp(ci*ckz*z1)
      cdum2=cdexp(-ci*ckz*z1)
      cf0=dcmplx(1.0,0.0)*cdexp(ci*phsin)
      do 20 i=1,2
      cdum3=cf0*(cdum1+cg*cdum2)
      cdum4=cf0*(cdum1-cg*cdum2)
      y(ie+1)=cdum3
      y(ie+2)=-ci*cdum3
      y(ie+3)=ci*ckz*cdum4
      y(ie+4)=ckz*cdum4
      y(ie+5)=0.0
      y(ie+6)=z1
      cf=dcmplx(y(ie+1),y(ie+2))
      cfp=dcmplx(y(ie+3),y(ie+4))
      kzr=ckz
      cdum10=cdexp(ci*kzr*z1)
      cdum20=cdexp(-ci*kzr*z1)
      cdum5=ci*kzr*cf
      cdum6=2.0*ci*kzr
      dum17=w*xmn**2*kcapmn/8.0d7
      delta=dsqrt(realc**2*rho(z1)/(2.0d0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      cffwd=(cdum5+cfp)*cdum20/cdum6
      cfbwd=(cdum5-cfp)*cdum10/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z1)
c       the backward component of cf is cfbwd*cdexp(-ci*kzr*z1)
      cfstar=dcmplx(y(ie+1),-y(ie+2))
      cfpstr=dcmplx(y(ie+3),-y(ie+4))
      cdum30=cf*cfpstr-cfstar*cfp
      az(1)=z1
      arw(1)=radius(z1)
      anrho(1)=rho(z1)/1.72e-8
      afamp(1)=dsqrt(y(ie+1)**2+y(ie+2)**2)
      afphse(1)=datan2(y(ie+2),y(ie+1))
      apfwd(1)=kzr*dum17*cdabs(cffwd)**2
      apbwd(1)=kzr*dum17*cdabs(cfbwd)**2
      apnet(1)=0.5*ci*dum17*cdum30
      apohmz(1)=dum17*dum18*dum19*(y(ie+1)**2+y(ie+2)**2)
      cf0=cf0*dsqrt(pin/apfwd(1))
   20 continue

      if(im.ne.0) go to 51
c store ptheta of all electrons at left end
      i=0
      do 50 it=1,itmax
      do 50 ipz=1,ipzmax
      do 50 irc=1,ircmax
      do 50 ipsi=1,ipsimx
      do 50 iphi=1,iphimx
      rc=y(i+1)
      psi=y(i+2)
      pp=y(i+3)
      phi=y(i+4)
      time=y(i+6)
      gamma=y(i+7)
      cf=dcmplx(y(ie+1),y(ie+2))
      call bfield(z1,bz,bzp)
      oe=e*bz/(me*realc)
      vp=pp/(gamma*me)
      oc=oe/gamma
      rl=vp/oc
      r=dsqrt(rc**2+rl**2+2.0*rc*rl*dsin(phi-psi))
      kmn=xmn/radius(z1)
      ptheta=0.5*me*oe*(rl**2-rc**2)
     &-(e/realc)*r*kmn*cf*bj(1,kmn*r)*cdexp(-ci*w*time)
      apti(it,ipz,irc,ipsi,iphi)=ptheta
      i=i+7
   50 continue
   51 continue

c integrate differential equations for coordinate/field array y from left 
c end to right end.
      do 10 i=1,ieqlst
   10 q(i)=0.0
      l=zmark(izmark)-z1
      delz=l/dfloat(izstep)
      do 100 iz=1,izstep
c advance coordinate/field array y by one step.
      call rkint(difeq,y,dy,q,1,ieqlst,delz)
      if(idiag.ne.1) go to 99
c store radius, resistivity, rf field, and powers as functions of z
c (wave power in w, wall loss per unit length in w/cm). 
c forward wave power (apfwd) and backward wave power (apbwd) are 
c artificially set to 1.0d-10 where w < wcmn.
c arrays az, arw, afamp, afphse, apfwd, apbwd, apnet, and apohmz are 
c passed to the main program through the common block.
      z=y(ie+6)
      rw=radius(z)
      cf=dcmplx(y(ie+1),y(ie+2))
      cfp=dcmplx(y(ie+3),y(ie+4))
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      az(iz+1)=z
      arw(iz+1)=radius(z)
      anrho(iz+1)=rho(z)/1.72e-8
      afamp(iz+1)=dsqrt(y(ie+1)**2+y(ie+2)**2)
      afphse(iz+1)=datan2(y(ie+2),y(ie+1))
      cfstar=dcmplx(y(ie+1),-y(ie+2))
      cfpstr=dcmplx(y(ie+3),-y(ie+4))
      cdum30=cf*cfpstr-cfstar*cfp
      apnet(iz+1)=0.5*ci*dum17*cdum30
      apohmz(iz+1)=dum17*dum18*dum19*(y(ie+1)**2+y(ie+2)**2)
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
      apfwd(iz+1)=1.0d-10
      apbwd(iz+1)=1.0d-10
      if(iz.ne.izstep) go to 39
      phsout=1.0d10
      delphs=1.0d10
   39 continue
*******************************************************************
c  diagnose electron phase space distribution (for subroutine pdiag)
      izmax=izstep+1
      azz(iz+1)=z
      i=0
      do 398 it=1,itmax
      do 398 ipz=1,ipzmax
      do 398 irc=1,ircmax
      do 398 ipsi=1,ipsimx
      do 398 iphi=1,iphimx
      if(it.ne.1) go to 397 
      if(ipz.ne.1) go to 397 
      if(irc.ne.1) go to 397 
      if(ipsi.ne.1) go to 397 
      app(it,ipz,irc,ipsi,iphi)=y(i+3)
      aphi(it,ipz,irc,ipsi,iphi)=y(i+4)
      phase=aphi(it,ipz,irc,ipsi,iphi)
      apx(iphi,iz+1)
     &=app(it,ipz,irc,ipsi,iphi)*dcos(phase)
      apy(iphi,iz+1)
     &=app(it,ipz,irc,ipsi,iphi)*dsin(phase)
      call bfield(z1,bz,bzp)
      oe=e*bz/(me*realc)
      oc0=oe/gamma0
      vz0=pz0/(gamma0*me)
      adphi(iz+1,iphi)=phase-oc0*(z-z1)/vz0
  397 continue
      i=i+7
  398 continue
c *****************************************************************
   99 continue
  100 continue
c store sum of wall losses (in w). pohmsm is passed to the main program 
c through the common block.
      pohmsm=y(ie+5)

      if(idiag.ne.1) go to 58
c store final coordinate arrays (arc,..agamma) at right end. these arrays
c are passed to the main program through the common block for diagnosis.
      gmcsv=0.0
      ptcsv=0.0
      i=0
      do 60 it=1,itmax
      do 60 ipz=1,ipzmax
      do 60 irc=1,ircmax
      do 60 ipsi=1,ipsimx
      do 60 iphi=1,iphimx
      rc=y(i+1)
      psi=y(i+2)
      pp=y(i+3)
      phi=y(i+4)
      pz=y(i+5)
      time=y(i+6)
      gamma=y(i+7)
      arc(it,ipz,irc,ipsi,iphi)=y(i+1)
      apsi(it,ipz,irc,ipsi,iphi)=y(i+2)
      app(it,ipz,irc,ipsi,iphi)=y(i+3)
      aphi(it,ipz,irc,ipsi,iphi)=y(i+4)
      apz(it,ipz,irc,ipsi,iphi)=y(i+5)
      atime(it,ipz,irc,ipsi,iphi)=y(i+6)
      agamma(it,ipz,irc,ipsi,iphi)=y(i+7)
c check constancy of ptheta(for m=0 mode only) and consistency of gamma.
      npp=pp/(me*realc)
      npz=pz/(me*realc)
      gamma2=dsqrt(1.0d0+npp**2+npz**2)
      agmcsv(it,ipz,irc,ipsi,iphi)=(gamma-gamma2)/gamma0
      gmcsv=gmcsv+dabs(agmcsv(it,ipz,irc,ipsi,iphi))/dfloat(itotal)
      if(im.ne.0) go to 59
      cf=dcmplx(y(ie+1),y(ie+2))
      z=y(ie+6)
      call bfield(z,bz,bzp)
      oe=e*bz/(me*realc)
      vp=pp/(gamma*me)
      oc=oe/gamma
      rl=vp/oc
      r=dsqrt(rc**2+rl**2+2.0*rc*rl*dsin(phi-psi))
      kmn=xmn/radius(z)
      ptheta=0.5*me*oe*(rl**2-rc**2)
     &-(e/realc)*r*kmn*cf*bj(1,kmn*r)*cdexp(-ci*w*time)
      aptf(it,ipz,irc,ipsi,iphi)=ptheta
      aptcsv(it,ipz,irc,ipsi,iphi)=(aptf(it,ipz,irc,ipsi,iphi)
     &-apti(it,ipz,irc,ipsi,iphi))/apti(it,ipz,irc,ipsi,iphi)
      ptcsv=ptcsv+dabs(aptcsv(it,ipz,irc,ipsi,iphi))/dfloat(itotal)
   59 continue
      i=i+7
   60 continue
   58 continue

c calculate cbc(cg).
      fr=y(ie+1)
      fi=y(ie+2)
      fpr=y(ie+3)
      fpi=y(ie+4)
      cf=dcmplx(fr,fi)
      cfp=dcmplx(fpr,fpi)
      zf=y(ie+6)
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
c **********************************************************************
      subroutine difeq(y,dy,ieqfst,ieqlst)
c this subroutine is called by subroutine rkint to calculate the 
c derivatives of the coordinate/field array y (defined as array dy).
c subroutine rkint is called by function cbc to integrate the differential
c equations for the coordinate/field array.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(500000),dy(500000)
      dimension weight(11,41,1,11,11)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension t(1000)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,is,im
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      ie=7*itotal
      z=y(ie+6)
      rw=radius(z)
      cf=dcmplx(y(ie+1),y(ie+2))
      cfp=dcmplx(y(ie+3),y(ie+4))
      call bfield(z,bz,bzp)
      oe=e*bz/(me*realc)
      oep=e*bzp/(me*realc)
      oepoe=oep/oe
      kmn=xmn/rw
      wcmn=kmn*realc
      cx=dcmplx(0.0d0,0.0d0)
      i=0
      do 60 it=1,itmax
      do 60 ipz=1,ipzmax
      do 60 irc=1,ircmax
      do 60 ipsi=1,ipsimx
      do 60 iphi=1,iphimx
      rc=y(i+1)
      psi=y(i+2)
      pp=y(i+3)
      phi=y(i+4)
      pz=y(i+5)
      time=y(i+6)
      gamma=y(i+7)
      vp=pp/(gamma*me)
      vz=pz/(gamma*me)
      oc=oe/gamma
      rl=vp/oc
      kmnrc=kmn*rc
      kmnrl=kmn*rl
      ism=is-im
      as=w*time-dfloat(is)*phi+dfloat(is-im)*psi
     &+(dfloat(im)-dfloat(is)/2.0d0)*pi
      cdum18=cdexp(-ci*as)
      cdum1=w*cf+ci*vz*cfp
      duma=kmn*e/(pz*realc)
ccc the following two statements model the converter section by tapering 
c   the electron charge in the first and last section
      if(z.gt.zmark(izmark-1)) 
     &duma=duma*(zmark(izmark)-z)/(zmark(izmark)-zmark(izmark-1))
      if(z.lt.zmark(2)) duma=duma*(z-zmark(1))/(zmark(2)-zmark(1))
      dumb=duma/rc
      dumc=duma*gamma*me
      dumd=dumc/pp
      dume=duma*pp
      dumf=duma*pp*w/(me*realc**2)
      call bes(iabs(ism)+2,kmnrc,0,result,t)
      if(ism.ge.0) dum11=t(ism+1)
      if(ism.lt.0) dum11=(-1)**ism*t(-ism+1)
      if((ism-1).ge.0) dum4=t(ism)
      if((ism-1).lt.0) dum4=(-1)**(ism-1)*t(2-ism)
      if((ism+1).ge.0) dum6=t(ism+2)
      if((ism+1).lt.0) dum6=(-1)**(ism+1)*t(-ism)
      call bes(iabs(is)+2,kmnrl,0,result,t)
      if(is.ge.0) dum3=t(is+1)
      if(is.lt.0) dum3=(-1)**is*t(-is+1)
      if((is-1).ge.0) dum5=t(is)
      if((is-1).lt.0) dum5=(-1)**(is-1)*t(2-is)
      if((is+1).ge.0) dum7=t(is+2)
      if((is+1).lt.0) dum7=(-1)**(is+1)*t(-is)
      dum2=0.5d0*(dum4-dum6)
      dum12=0.5d0*(dum5-dum7)
      dum23=dum2*dum3
      dum45=dum4*dum5
      dum67=dum6*dum7
      dum113=dum11*dum3
      dum112=dum11*dum12
      dy(i+1)=duma*ci*(-cdum1*dum23/oc
     &+0.5d0*kmnrl*cf*(dum45-dum67))*cdum18
     &-0.5d0*oepoe*(rc-rl*dsin(phi-psi))
      dy(i+2)=dumb*(-cdum1*dfloat(is-im)*dum113/(oc*kmnrc)
     &+0.5d0*kmnrl*cf*(dum45+dum67))*cdum18
     &-0.5d0*oepoe*(rl/rc)*dcos(phi-psi)
      dy(i+3)=dumc*ci*cdum1*dum112*cdum18
     &+0.5d0*me*oep*(rl+rc*dsin(phi-psi))
      dy(i+4)=dumd*(-dfloat(is)*cdum1/kmnrl+kmn*vp*cf)*dum113*cdum18
     &+me*oe/pz+0.5d0*me*(oep/pp)*rc*dcos(phi-psi)
      dy(i+5)=dume*cfp*dum112*cdum18
     &-0.5d0*me*oep*(pp/pz)*(rl+rc*dsin(phi-psi))
      dy(i+6)=1.0d0/vz
      dy(i+7)=dumf*ci*cf*dum112*cdum18
      dumg=8.0d0*kmn*rib/(xmn**2*kcapmn*realc)
      cx=cx+dumg*weight(it,ipz,irc,ipsi,iphi)*(pp/pz)
     &*dum112*cdexp(ci*as)
      i=i+7
   60 continue
ccc the following two statements model the converter section by tapering 
c   the electron charge in the first and last section
      if(z.gt.zmark(izmark-1)) 
     &cx=cx*(zmark(izmark)-z)/(zmark(izmark)-zmark(izmark-1))
      if(z.lt.zmark(2)) cx=cx*(z-zmark(1))/(zmark(2)-zmark(1))
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      kz2r=ckz2
      kz2i=-ci*ckz2
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      dum17=w*xmn**2*kcapmn/8.0d7
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0d0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      xr=cx
      xi=-ci*cx
      dy(ie+1)=y(ie+3)
      dy(ie+2)=y(ie+4)
      dy(ie+3)=-kz2r*y(ie+1)+kz2i*y(ie+2)-xr
      dy(ie+4)=-kz2i*y(ie+1)-kz2r*y(ie+2)-xi
      dy(ie+5)=dum17*dum18*dum19*(y(ie+1)**2+y(ie+2)**2)
      dy(ie+6)=1.0d0
      return
      end
c    ***********************************************************************
      subroutine cldtst(fhzmin,fhzmax,iwmax)
c    ***********************************************************************
c    this subroutine calculates and plots the return and transmission losses
c    of the interaction structure as a function of the frequency, in response
c    to an incident wave of constant power continuously injected from the
c    left end into the interaction structure.
c
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension cg(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension afreq(20001),atldb(20001),arldb(20001)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,is,im
      common/local1/pbwdi,pfwdf
      external cbcckt
c
      if(iwmax.eq.1) return
      z1=zmark(1)
      rw=radius(z1)
      wmin=fhzmin*2.0d0*pi
      if((wmin/realc)**2.le.xmn**2/rw**2) write(*,1)
      if((wmin/realc)**2.le.xmn**2/rw**2) return
    1 format(' * fhzmin below cutoff in first section, no cold test *')
c
c    specify no of steps for z-integration.
      izstep=2001
c    specify incident (forward) wave power in w at left end
      pin=1.0
      tlmin=0.0
      rlmin=0.0
      do 100 iw=1,iwmax
      fstep=(fhzmax-fhzmin)/dfloat(iwmax-1)
      fhz=fhzmin+fstep*(iw-1)
      afreq(iw)=fhz
      w=2.0*pi*fhz
c    angular frequency of the incident wave (in radian/sec, denoted by w)
c    is passed to subprograms through the common block.
c
c    guess complex reflection coefficient
      cguess=dcmplx(0.5,0.5)
      if(iw.eq.1) cg(1)=cguess

c    calculate reflection coefficient (referred to the z=0 plane)
      ep1=1.0d-5
      ep2=ep1
      iroot=1
      imaxit=50
      call muller(0,iroot,cg,imaxit,ep1,ep2,cbcckt,.false.)
      value=cdabs(cbcckt(cg(1)))
      tldb=10.0*dlog10(pfwdf/pin)
      rldb=10.0*dlog10(pbwdi/pin)
      atldb(iw)=tldb
      arldb(iw)=rldb
      if(atldb(iw).lt.tlmin) tlmin=atldb(iw)
      if(arldb(iw).lt.rlmin) rlmin=arldb(iw)
  100 continue
      call sscale(afreq(1),afreq(iwmax),tlmin,0.0d0)
      call splot(afreq,atldb,iwmax,
     &'tt  transmission loss(db) vs freq.(hz)$')
      call sscale(afreq(1),afreq(iwmax),rlmin,0.0d0)
      call splot(afreq,arldb,iwmax,
     &'rr  return loss(db) vs freq.(hz)$')
      return
      end
c ******************************************************************
      function cbcckt(cg)
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(20),dy(20),q(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode/w,pin,xmn,kcapmn,is,im
      common/local1/pbwdi,pfwdf
      external difckt
c
c initialize rf field array y at left end according to the incident
c forward wave power (pin) specified in the main program and store
c backward wave power at left end
      z1=zmark(1)
      rw=radius(z1)
      ckz2=(w/realc)**2-xmn**2*closs(z1)/rw**2
      ckz=cdsqrt(ckz2)
      cdum1=cdexp(ci*ckz*z1)
      cdum2=cdexp(-ci*ckz*z1)
      cf0=dcmplx(1.0,0.0)
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
      cffwd=(cdum5+cfp)*cdum20/cdum6
      cfbwd=(cdum5-cfp)*cdum10/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z1)
c       the backward component of cf is cfbwd*cdexp(-ci*kzr*z1)
      pfwdi=kzr*dum17*cdabs(cffwd)**2
      pbwdi=kzr*dum17*cdabs(cfbwd)**2
      cf0=cf0*dsqrt(pin/pfwdi)
   20 continue
c integrate differential equations for rf field array y from left end
c to right end.
      do 10 i=1,6
   10 q(i)=0.0
      l=zmark(izmark)-z1
      delz=l/dfloat(izstep)
      do 100 iz=1,izstep
c advance rf field array y by one step.
      call rkint(difckt,y,dy,q,1,6,delz)
  100 continue
c store forward power at final step 
      z=y(6)
      rw=radius(z)
      cf=dcmplx(y(1),y(2))
      cfp=dcmplx(y(3),y(4))
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      ckz=cdsqrt(ckz2)
      kzr=ckz
      cdum10=cdexp(ci*kzr*z)
      cdum20=cdexp(-ci*kzr*z)
      cdum5=ci*kzr*cf
      cdum6=2.0*ci*kzr
      cffwd=(cdum5+cfp)*cdum20/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z)
      pfwdf=kzr*dum17*cdabs(cffwd)**2
c calculate cbcckt(cg).
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
      cbcckt=cfp-ci*ckz*cf
      go to 201
  202 continue
      ckapa2=xmn**2*closs(zf)/rwf**2-(w/realc)**2
      ckapa=cdsqrt(ckapa2)
      cbcckt=cfp+ckapa*cf
  201 continue
      return
      end
c ******************************************************************
      subroutine difckt(y,dy,ieqfst,ieqlst)
c this subroutine is called by subroutine rkint to calculate the
c derivatives of rf field array y (defined as array dy). subroutine
c rkint is called by function cbcckt to integrate the differential
c equations for the rf field array.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(20),dy(20)
      common/const/ci,pi,realc,me,e
      common/mode/w,pin,xmn,kcapmn,is,im
      z=y(6)
      rw=radius(z)
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      kz2r=ckz2
      kz2i=-ci*ckz2
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
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
c **********************************************************************
c               ---  double precision version  ---
c **********************************************************************
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
c
c default values
      if(z.lt.zmark(1)) radius=rwl(1)
      if(z.gt.zmark(izmark)) radius=rwr(izmark)
c
      return
      end
c **********************************************************************
      function rho(z)
c **********************************************************************
c               ---  double precision version  ---
c **********************************************************************
c rho(z) is the wall resistivity at position z.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      imax=izmark-1
      do 100 i=1,imax
      if(z.ge.zmark(i).and.z.le.zmark(i+1))
     &rho=rhor(i)+(rhol(i+1)-rhor(i))*(z-zmark(i))/(zmark(i+1)-zmark(i))
  100 continue
c
c default values
      if(z.lt.zmark(1)) rho=rhol(1)
      if(z.gt.zmark(izmark)) rho=rhor(izmark)
c
      return
      end
c **********************************************************************
      function closs(z)
c **********************************************************************
c               ---  double precision version  ---
c **********************************************************************
c closs(z) is a factor to account for wall losses of the te(m,n) mode
c at position z. closs(z)=1.0 for zero wall resistivity.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      common/const/ci,pi,realc,me,e
      common/mode/w,pin,xmn,kcapmn,is,im
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      rw=radius(z)
      wcmn=xmn*realc/rw
      cdum1=(1.0+ci)*delta/rw
      cdum2=(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w**2/wcmn**2)
      closs=1.0-cdum1*(1.0+cdum2)
      return
      end
c **********************************************************************
      subroutine bfield(z,bz,bzp)
c this subroutine calculates the external b-field (on the axis) and its
c derivative at position z.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      common/bext/bz0
      bz=bz0
      bzp=0.0
      return
      end
c ***********************************************************************
      subroutine coord(arc,apsi,app,aphi,apz,atime,agamma,weight)
      implicit  real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension arc(11,41,1,11,11),apsi(11,41,1,11,11),
     &app(11,41,1,11,11),aphi(11,41,1,11,11),apz(11,41,1,11,11),
     &atime(11,41,1,11,11),agamma(11,41,1,11,11),
     &weight(11,41,1,11,11)
      common/const/ci,pi,realc,me,e
      common/mode/w,pin,xmn,kcapmn,is,im
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
c
c this subroutine sets initial electron coordinates (in real and momentum  
c spaces) and weighing factors
      suma=0.0
      tstep=2.0*pi/(w*dfloat(itmax))
      phistp=2.0d0*pi/dfloat(iphimx)
      if(im.ne.0) psistp=2.0d0*pi/dfloat(ipsimx*im)
      pzstep=(pzmax-pzmin)/dfloat(ipzmax)
      rcstep=(rc2-rc1)/dfloat(ircmax)
      do 120 it=1,itmax
      do 120 ipz=1,ipzmax
      pz=pzmin+pzstep*(dfloat(ipz)-0.5)
      pp=dsqrt(p0-pz)*dsqrt(p0+pz)
      do 120 irc=1,ircmax
      rc=rc1+rcstep*(dfloat(irc)-0.5)
      if(rc.eq.0.0d0) rc=1.0d-7
      if(delpz.eq.0.0) fac2=rc
      if(delpz.ne.0.0)
     &fac2=dexp(-0.5*((pz-pz0)/delpz)**2)*rc
      do 120 ipsi=1,ipsimx
      do 120 iphi=1,iphimx
      arc(it,ipz,irc,ipsi,iphi)=rc
      if(ipsimx.eq.1) apsi(it,ipz,irc,ipsi,iphi)=0.0d0
      if(ipsimx.ne.1)
     &apsi(it,ipz,irc,ipsi,iphi)=psistp*(dfloat(ipsi)-0.5)
      app(it,ipz,irc,ipsi,iphi)=pp
      aphi(it,ipz,irc,ipsi,iphi)=phistp*(dfloat(iphi)-0.5d0)
     &+apsi(it,ipz,irc,ipsi,iphi)
c    &+2.0d0*pi*dfloat(it-1)/dfloat(itmax)
c    &+2.0d0*pi*dfloat(ipz-1)/dfloat(ipzmax)
c    &+2.0d0*pi*dfloat(irc-1)/dfloat(ircmax)
      apz(it,ipz,irc,ipsi,iphi)=pz
      atime(it,ipz,irc,ipsi,iphi)=tstep*dfloat(it-1)
      agamma(it,ipz,irc,ipsi,iphi)=gamma0
      weight(it,ipz,irc,ipsi,iphi)=fac2
      suma=suma+fac2
  120 continue
      do 130 it=1,itmax
      do 130 ipz=1,ipzmax
      do 130 irc=1,ircmax
      do 130 ipsi=1,ipsimx
      do 130 iphi=1,iphimx
      weight(it,ipz,irc,ipsi,iphi)=weight(it,ipz,irc,ipsi,iphi)/suma
  130 continue
      return
      end
c *********************************************************************
      function oprc(is,im,in,ipeak)
      implicit  real*8  (a-h,k-z)
      dimension axmn(9,8)
      data axmn/
     & 3.832d0, 1.841d0, 3.054d0, 4.2 1d0, 5.318d0, 6.416d0,
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
      xmn=axmn(iabs(im)+1,in)
      oprc=axmn(iabs(is-im)+1,ipeak)/xmn
      if(is.eq.im.and.ipeak.eq.1) oprc=0.0d0
      if(is.eq.im.and.ipeak.ne.1) oprc=axmn(1,ipeak-1)/xmn
      return
      end
c *********************************************************************
      subroutine couple(is,im,x,y,hsm,tsm,usm)
c ******************************************************************
c               ---  double precision version  ---
c ******************************************************************
      implicit  real*8 (a-h,j-z)
      xx=x
      if(x.eq.0.0e0) xx=1.0d-6
      d1=bjp(is,y)
      d2=bj(is-im,xx)
      d3=bjp(is-im,xx)
      d4=bj(is-im-1,xx)
      d5=bj(is-im+1,xx)
      hsm=d2**2*d1**2
      tsm=2.0d0*hsm+y*d1*(2.0d0*bjpp(is,y)*d2**2
     &-bj(is,y)*(d2*d3/xx+d3**2+d2*bjpp(is-im,xx)))
      usm=-0.5d0*y*d1*(bj(is-1,y)*(d4*d2/xx+bjp(is-im-1,xx)*d2+d4*d3)
     &-bj(is+1,y)*(d5*d2/xx+bjp(is-im+1,xx)*d2+d5*d3))
      return
      end
c **********************************************************************
      subroutine rkint (derivy, y, dy, q, neqfst, neqlst, dx)
c **********************************************************************
c               ---  double precision version  ---
c **********************************************************************
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
c ********************************************************************
c               ---  double precision version  ---
c ********************************************************************
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
c *********************************************************************
c               ---  double precision version  ---
c *********************************************************************
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
    9 format(1pd14.4,1x,101a1)
   10 format(15x,101a1)
      write(*,12)  (xx(n),n=1,101,20)
   12 format(1p6e20.4)
      return
      end
c *********************************************************************
      subroutine splot(x, y, npt, tit)
c *********************************************************************
c               ---  double precision version  ---
c *********************************************************************
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
    9 format(1x,1pd11.4,101a1)
   10 format(12x,101a1)
      write(*,12)  (xx(n),n=1,nxap,nxp2)
   12 format(1x,1pd17.4,1p5e12.4)
      return
      end
c ******************************************************************
      subroutine mxmn (y,j ,nmax, nmin)
c ******************************************************************
c               ---  double precision version  ---
c ******************************************************************
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
      if(dx.eq.0.) dx =dabs(x1)
      if(dx.eq.0.) dx = 1
    2 ix =dlog10(dx)
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
c ***********************************************************************
c               ---  double precision version  ---
c ***********************************************************************
c
c bj: bessel function j(x) of integer order.
c in: positive or negative order of bj (-40< in <40 ).
c x : positive or negative argument of bj (x=0 or 1.0e-10<x <300 ).
c
      implicit  real*8  (a,b,d-h,j-z)
      dimension t(1000)
      isign1=1
      isign2=1
c     sx=sngl(x)
c     modified by c.s.hsue
      xabs=dabs(x)
      inabs=iabs(in)
      call bes(inabs,xabs,0,sbj,t)
      if(in.lt.0) isign1=(-1)**inabs
      if(x.lt.0.0e0) isign2=(-1)**inabs
      bj=sbj*dfloat(isign1*isign2)
      return
      end
c **********************************************************************
      function bjp(in,x)
c **********************************************************************
c               ---  double precision version  ---
c **********************************************************************
c bjp=first derivative of bj(x)
      implicit  real*8 (a,b,d-h,j-z)
      bjp=0.5d0*(bj(in-1,x)-bj(in+1,x))
      return
      end
c **********************************************************************
      function bjpp(in,x)
c **********************************************************************
c               ---  double precision version  ---
c **********************************************************************
c bjpp=second derivative of bj(x)
      implicit  real*8 (a,b,d-h,j-z)
      bjpp=0.25d0*(bj(in-2,x)-2.0d0*bj(in,x)+bj(in+2,x))
      return
      end
c *********************************************************************
      subroutine bes(no,x,kode,result,t)
c *********************************************************************
c               ---  double precision version  ---
c *********************************************************************
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
      t(lub) = 1.d-75
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
************************************************************
      subroutine pdiag(ipypx)

      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension adphi(20001,301),apx(301,20001),
     &apy(301,20001),az(20001),azz(20001)
      common/diagp/adphi,apx,apy,azz,izmax
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
c
c this subprogram plots delta phi (deviation of electron phase
c angle with respect to the zero order orbit phase angle) of 
c representative electrons vs z, and py vs px at selected z
c coordinates.
c
c calculations are done in subroutine cbc and results are
c passed to here through common statement denoted "diagp".
c magnetic field is assumed to be constant along z for the 
c calculation of delta phi.
c
c only those electrons with it=ipz=irc=ipsi=1 are diagnosed.
c 
c ipypx is a user supplied integer specifying the no of axial
c positions at which py vs px will be plotted. py vs px at the
c 2 end positions will be automatically plotted.
c
      apmax=0.0
      aphsmx=0.0
      aphsmn=0.0
      do 396 iz=1,izmax
      az(iz)=azz(iz)
      do 396 iphi=1,iphimx
      if(adphi(iz,iphi).gt.aphsmx) 
     &aphsmx=adphi(iz,iphi)
      if(adphi(iz,iphi).lt.aphsmn) 
     &aphsmn=adphi(iz,iphi)
      if(dabs(apx(iphi,iz)).gt.apmax) 
     &apmax=dabs(apx(iphi,iz))
      if(dabs(apy(iz,iphi)).gt.apmax) 
     &apmax=dabs(apy(iz,iphi))
  396 continue
      call sscale(az(1),az(izmax),aphsmn,aphsmx)
      do 393 iphi=1,iphimx
      call splot(az,adphi(1,iphi),izmax,'p   $')
  393 continue
      call splot(az,adphi(1,iphimx),izmax,
     &'pp  linear phase diagram$')
      call sscale(az(1),az(izmax),aphsmn,aphsmx)
      do 392 iphi=1,iphimx,5
      call splot(az,adphi(1,iphi),izmax,'p   $')
  392 continue
      call splot(az,adphi(1,iphimx),izmax,
     &'pp  linear phase diagram$')
      call sscale(az(1),az(izmax),aphsmn,aphsmx)
      do 391 iphi=1,iphimx,10
      call splot(az,adphi(1,iphi),izmax,'p   $')
  391 continue
      call splot(az,adphi(1,iphimx),izmax,
     &'pp  linear phase diagram$')
      iskip=idint(dfloat(izmax-1)/dfloat(ipypx))
      if(iskip.eq.0) iskip=1
      do 395 iz=1,izmax,iskip
      call sscale(-apmax,apmax,-apmax,apmax)
      call splot(apx(1,iz),apy(1,iz),iphimx,
     &'pp  polar phase diagram$')
      z=az(iz)
      write(6,394) z
  394 format(' z=',1pd10.3,' cm')
  395 continue
      if((izmax-1).eq.(iskip*ipypx)) go to 401
      call sscale(-apmax,apmax,-apmax,apmax)
      call splot(apx(1,izmax),apy(1,izmax),iphimx,
     &'pp  polar phase diagram$')
      z=az(izmax)
      write(6,394) z
  401 continue
      return
      end

