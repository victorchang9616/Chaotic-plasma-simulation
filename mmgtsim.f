      program mmgtsim
c   
c   (all in double precision)
c   ***********************************************************************
c
c   this program calculates the gain, growth rate, output power, efficiency,
c   reflection coefficient, etc. of a multi-mode gyro-twt or carm amplifier  
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
c   the program can handle an arbitrary number of te(m,n) modes. since it
c   looks for a quasi steady-state solution(i.e. everything repeats itself 
c   over a finite period), the use of multiple modes is rigorously justified
c   only when these modes have harmonic frequencies(i.e. higher frequencies
c   are harmonics of the lowest frequency). in such cases, the period over
c   which the solution repeats itself is the period of the lowest frequency 
c   mode, and the quasi steady-state solution is obtained by averaging the 
c   electron performance over this period.
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
c   date of first version: july 31, 1995
c   date of this version: september 25, 2000
c   ************************************************************************
c
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
c   all variables beginning with a,b,d-h,j-z are real numbers. all
c   variables beginning with i are integers. all variables beginning 
c   with c are complex numbers.
      dimension cg(20)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100) 
      dimension az(20001),arw(20001),anrho(20001)
      dimension afamp(20001,3),afphse(20001,3),apnet(20001,3),
     &apfwd(20001,3),apbwd(20001,3),apohmz(20001,3)
      dimension agfwd(20001,3),agfwdz(20001,3),agbwd(20001,3),
     &agohm(20001,3)
      dimension atz(101),afreq(101,3),apcnsv(101,101),aetab(101,101),
     &apin1(101,101),apin2(101,101,3),apfwd1(101,101,3),
     &apfwd2(101,101,3),aefwd1(101,101,3),aefwd2(101,101,3),
     &aebwd2(101,101,3),aeohm2(101,101,3),agbwd0(101,101,3),
     &agbwd2(101,101,3),agfwd0(101,101,3),agfwd1(101,101,3),
     &agfwd2(101,101,3),icont2(101,101,3),avalu2(101,101,3)
c   arc(it,ipz,irc,ipsi,iphi)
      dimension arc(21,41,1,1,23),apsi(21,41,1,1,23),
     &app(21,41,1,1,23),aphi(21,41,1,1,23),apz(21,41,1,1,23),
     &atime(21,41,1,1,23),agamma(21,41,1,1,23),
     &weight(21,41,1,1,23)
      dimension apti(21,41,1,1,23),aptf(21,41,1,1,23),
     &aptcsv(21,41,1,1,23),agmcsv(21,41,1,1,23)
      dimension itot(100)
      dimension axmn(9,8),aw(11),apin(11),iis(11),iim(11),iin(11)        
      dimension arwref(11),pohmsm(11),pohmdb(11),etaohm(11),
     &etafwd(11),etabwd(11),ratio1(11),ratio2(11),gfwd(11),
     &fcldmn(11),fcldmx(11),nwave(11),lamdag(11),
     &fcref(11),bzg(11),fg(11)
      dimension phsin(11),phsout(11),delphs(11)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode1/axmn,aw,apin,iis,iim,iin,imdmax,imdbcl,imdbcr
      common/mode2/w,pin,xmn,kcapmn,is,im
      common/phase/phsin,phsout,delphs
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      common/bext/bz0
      common/cutoff/icut
      common/diagf/az,arw,anrho,afamp,afphse,apnet,apfwd,apbwd,apohmz,
     &pohmsm,idiag,icont
      common/diagb/arc,apsi,app,aphi,apz,atime,agamma
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
      iplot=0
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
c
ccc specify method of the rf field calculation in cutoff regions
      icut=1
c   icut=1: artificially set rf field to 0 in cutoff regions (this procedure 
c           may be necessary in long cutoff regions, because the field
c           equation with wrong initial conditions at the left end of the
c           cutoff region will yield exponentially growing solutions (which 
c           is the evanescent wave solution for a wave incident from the
c           right end of the cutoff region).
c   icut=0: continue rf field calculations in cutoff regions
c
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
ccc specify rf structure dimension arrays zmark, rwl, and rwr in cm. 
      r1=2.044d0
      r2=1.45d0
      rout=1.75d0
      tdeg=1.5d0
      la=9.3d0
      lb=5.0d0
      do 1000 ilc=1,1
      lc=31.0d0
      ld=7.0d0
      le=8.0d0
      lf=(rout-r2)/dtan(tdeg*2.0d0*pi/360.0d0)
      lg=1.0d0
      write(6,113) r1,r2,rout,tdeg,la,lb,lc,ld,le,lf,lg
  113 format(' r1=',1pd10.3,' cm, r2=',1pd10.3,
     &' cm, rout=',1pd10.3,' cm, tdeg=',1pd10.3,' degree'/
     &' la=',1pd10.3,' cm, lb=',1pd10.3,' cm, lc=',1pd10.3,
     &' cm, ld=',1pd10.3,' cm'/' le=',1pd10.3,' cm, lf=',1pd10.3,
     &' cm, lg=',1pd10.3,' cm')
      zmark(1)=0.0d0
      zmark(2)=zmark(1)+la
      zmark(3)=zmark(2)+lb
      zmark(4)=zmark(3)+lc
      zmark(5)=zmark(4)+ld
      zmark(6)=zmark(5)+le
      zmark(7)=zmark(6)+lf
      zmark(8)=zmark(7)+lg
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
      rwl(7)=rout
      rwr(7)=rwl(7)
      rwl(8)=rout
      rwr(8)=rwl(8)

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
      rhocu=1.72d-8
      rho1=2.0d5*rhocu
      rho2=5.7d4*rhocu
c
      rhol(1)=rho1
      rhor(1)=rho1
      rhol(2)=rho1
      rhor(2)=rho1
      rhol(3)=rho2
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
      write(6,3) izmark
    3 format(/' rf structure dimension arrays: (izmark=',i4,')')
      do 50 i=1,izmark,6
      if(izmark.ge.(i+5)) imax=i+5
      if(izmark.lt.(i+5)) imax=izmark
      write(6,4) (zmark(ii),ii=i,imax)
    4 format(/'   zmark(cm)=',6(1pd10.3,','))
      write(6,5) (rwl(ii),ii=i,imax)
    5 format('     rwl(cm)=',6(1pd10.3,','))
      write(6,6) (rwr(ii),ii=i,imax)
    6 format('     rwr(cm)=',6(1pd10.3,','))
   50 continue
      write(6,7) izmark
    7 format(/' wall resistivity arrays: (izmark=',i4,')')
      do 60 i=1,izmark,6
      if(izmark.ge.(i+5)) imax=i+5
      if(izmark.lt.(i+5)) imax=izmark
      write(6,8) (zmark(ii),ii=i,imax)
    8 format(/'   zmark(cm)=',6(1pd10.3,','))
      write(6,9) (rhol(ii),ii=i,imax)
    9 format(' rhol(ohm-m)=',6(1pd10.3,','))
      write(6,11) (rhor(ii),ii=i,imax)
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

ccc specify the total number of modes (each mode is characterized by
c   mode numbers "im", "in", and cyclotron harmonic number "is" at which
c   the te(im,in) mode interacts with the electron beam)
      imdmax=2
ccc specify the mode with reference to which the normalized beam guiding
c   center positions nrcavg, ndrc, nrc1, nrc2, etc. are to be converted
c   into (unnormalized) physical quantities
      imdrc=1
ccc specify the mode number with reference to which the normalized magnetic
c   field bz0bzg is to be converted into a (unnormalized) physical quantity
      imdbz=1
ccc specify the mode number whose input power (denoted by pref) will be 
c   used for gain calculations of all modes (note: pref must be > 0.0)
      imdg=1
ccc specify the mode number whose output power will be maximized (if needed)
      imdp=2
ccc specify the mode number whose reflection coefficient cg (when pin>0)
c   or backward wave amplitude cg (when pin=0) at the left end will be
c   iteratively evaluated (if itera=1) by muller to satisfy the boundary 
c   condition cbc(cg)=0 at the right end for mode imdbcr 
      imdbcl=2
ccc specify the mode number whose cbc value will be evaluated and, if itera=1,
c   minimized by muller to satisfy the boundary condition at the right end
      imdbcr=2
c 
ccc specify im and in of each te(im,in) mode and its associated cyclotron
c   harmonic number is
c
      iim(1)=0
      iin(1)=2
      iis(1)=1

      iim(2)=0
      iin(2)=3
      iis(2)=2
c 
ccc specify the reference wall radius of each mode for the calculation of
c   cut-off frequency, grazing magnetic field, intersection points between
c   waveguide mode and beam-wave resonance line, etc.
c
c   note: the values of arwref(imdrc) and arwref(imdbz) will affect the values
c         of rcavg and bz0 etc. and hence the beam dynamics. other values of
c         arwref will not affect the beam dynamics
c
      arwref(1)=r1
      arwref(2)=r2
c
ccc instruction for cold test of the circuit
c   (icold=1: yes, icold=0: no)
      icold=0
      iwmax=41
      fcldmn(1)=16.0d9
      fcldmx(1)=21.0d9
      
      fcldmn(2)=34.0d9
      fcldmx(2)=39.0d9
      do 156 imode=2,imdmax
      if(icold.eq.1) 
     &call cldtst(imode,fcldmn(imode),fcldmx(imode),iwmax)
c   note: for this sample program, section lc is cutoff to mode 1.
  156 continue
c
ccc specify normalized beam electron guiding center position (nrcavg) and
c   normalized guiding center spread (ndrc). both are normalized to rwref
c   of mode imdrc.
c   (rc can not be set to 0 because that would make psi meaningless)
      do 1000 ircavg=1,1
      nrcavg=0.2624d0
      ndrc=0.0d0
c     rcavg=0.54d0
c     drc=0.0d0
c     nrcavg=rcavg/arwref(imdrc)
c     ndrc=drc/arwref(imdrc)
c
c   determine beam guiding center positions rcavg, drc, rc1, rc2, etc. in cm.
      nrc1=nrcavg-0.5d0*ndrc
      nrc2=nrcavg+0.5d0*ndrc
      rwref=arwref(imdrc)
      rcavg=nrcavg*rwref
      drc=ndrc*rwref
      nrc1=nrcavg-0.5d0*ndrc
      nrc2=nrcavg+0.5d0*ndrc
      rc1=nrc1*rwref
      rc2=nrc2*rwref
c
ccc specify beam voltage in kv
c   (all electrons are assumed to have the same initial energy)
      do 1000 ivbkv=1,1
      vbkv=60.0d0
c
ccc specify beam current in amp
      do 1000 irib=1,1
      ribamp=13.0d0
c
ccc specify beam alpha0
c   (alpha0 can not be set to 0 because that would make phi meaningless)
      do 1000 ialpha=1,1
      alpha0=1.5d0
c
c   calculate beam parameters
      gamma0=1.0d0+vbkv/511.0d0
      theta0=datan2(alpha0,1.0d0)
      angle=datan2(alpha0,1.0d0)*180.0d0/pi
      p0=dsqrt(gamma0**2-1.0d0)*me*realc
      pp0=p0*dsin(theta0)
      pz0=p0*dcos(theta0)
      v0=p0/(gamma0*me)
      vp0=v0*dsin(theta0)
      vz0=v0*dcos(theta0)
      gamz0=1.0d0/dsqrt(1.0d0-(vz0/realc)**2)
      rib=ribamp*3d9
      pbkw=vbkv*ribamp
c
c   for each mode, determine its cutoff frequency (fcref), grazing magnetic
c   field (bzg) and grazing frequency (fg), etc. w.r.t. the reference wall
c   radius of this mode.
      do 251 imode=1,imdmax
      im=iim(imode)
      in=iin(imode)
      is=iis(imode)
      xmn=axmn(iabs(im)+1,in)
      rwref=arwref(imode)
      wcref=xmn*realc/rwref
      fcref(imode)=wcref/(2.0d0*pi)
      oeg=gamma0*wcref/(gamz0*dfloat(is))
      bzg(imode)=oeg*me*realc/e
      fg(imode)=gamz0*fcref(imode)
  251 continue
c
ccc specify external magnetic field (bz0) in gauss or in ratio of bz0 to bzg
c   of mode imdbz.
      do 1000 ib=1,1
      bz0bzg=1.012d0+0.002d0*(ib-1)
      bz0=bz0bzg*bzg(imdbz)
c     bz0=1.232d4
c     bz0bzg=bz0/bzg(imdbz)
c
      oe=e*bz0/(me*realc)
      oc=oe/gamma0
      fcyc=oc/(2.0d0*pi)
      rl0=vp0/oc
      write(6,13)vbkv,ribamp,pbkw,alpha0,angle,gamma0,gamz0
   13 format(/' vb=',1pd12.5,' kv, ib=',1pd10.3,' amp, pb=',
     &1pd10.3,' kw',/' alpha0=',1pd10.3,'(',1pd11.4,
     &' deg), gamma0=',1pd13.6,', gamz0=',1pd10.3/)
      write(6,257) bz0,fcyc,rcavg,drc,rl0
  257 format(' bz0=',1pd12.5,' gauss, fcyc=',1pd9.2,
     &' hz, rcavg=',1pd9.2,' cm, drc=',1pd9.2,' cm'/
     &' rl0=',1pd10.3,' cm')
c
c   calculate and print key information relevant to each mode
      write(6,254)
  254 format(' mode im in is    rwref(cm)  rcavg/rwref',
     &'    fcref(hz)     fg(hz)      bz0/bzg')
      do 256 imode=1,imdmax
      nrcavg=rcavg/arwref(imode)
      bz0bzg=bz0/bzg(imode)
      if(imode.eq.imdrc.and.imode.ne.imdbz)
     &write(6,258) imode,iim(imode),iin(imode),iis(imode),
     &arwref(imode),nrcavg,fcref(imode),fg(imode),bz0bzg
      if(imode.ne.imdrc.and.imode.eq.imdbz)
     &write(6,259) imode,iim(imode),iin(imode),iis(imode),
     &arwref(imode),nrcavg,fcref(imode),fg(imode),bz0bzg
      if(imode.eq.imdrc.and.imode.eq.imdbz)
     &write(6,260) imode,iim(imode),iin(imode),iis(imode),
     &arwref(imode),nrcavg,fcref(imode),fg(imode),bz0bzg
      if(imode.ne.imdrc.and.imode.ne.imdbz)
     &write(6,261) imode,iim(imode),iin(imode),iis(imode),
     &arwref(imode),nrcavg,fcref(imode),fg(imode),bz0bzg
  258 format(1x,i3,':',i3,i3,i3,2x,1pd11.4,' <',1pd11.4,'> ',
     &1pd11.4,2x,1pd11.4,2x,1pd11.4)
  259 format(1x,i3,':',i3,i3,i3,2x,1pd11.4,2x,1pd11.4,2x,
     &1pd11.4,2x,1pd11.4,' <',1pd11.4,'>')
  260 format(1x,i3,':',i3,i3,i3,2x,1pd11.4,' <',1pd11.4,'> ',
     &1pd11.4,2x,1pd11.4,' <',1pd11.4,'>')
  261 format(1x,i3,':',i3,i3,i3,2x,1pd11.4,2x,1pd11.4,2x,
     &1pd11.4,2x,1pd11.4,2x,1pd11.4)
  256 continue
      write(6,252)
  252 format(' ')
c   calculate and print detailed information relevant to each mode
      do 163 imode=1,imdmax
      im=iim(imode)
      in=iin(imode)
      is=iis(imode)
      xmn=axmn(iabs(im)+1,in)
      rwref=arwref(imode)
      nrcavg=rcavg/rwref
      ndrc=drc/rwref
      nrc1=rc1/rwref
      nrc2=rc2/rwref
      opnrc1=oprc(is,im,in,1)
      opnrc2=oprc(is,im,in,2)
      opnrc3=oprc(is,im,in,3)
      bz0bzg=bz0/bzg(imode)
c   calculate points of intersections between waveguide mode (radius=rwref)
c   and beam-wave resonance line.
      wcref=2.0d0*pi*fcref(imode)
      if(bz0bzg.gt.1.0) go to 103
      kzg=gamz0*wcref*vz0/realc**2
      f1=fg(imode)
      f2=fg(imode)
      kz1=kzg
      kz2=kzg
      go to 104
  103 continue
      dum1=1.0d0-(wcref/(dfloat(is)*oc*gamz0))**2
      f1=dfloat(is)*oc*gamz0**2*(1.0-(vz0/realc)*dsqrt(dum1))/(2.0*pi)
      f2=dfloat(is)*oc*gamz0**2*(1.0+(vz0/realc)*dsqrt(dum1))/(2.0*pi)
      kz1=dfloat(is)*(oc/realc)*gamz0**2*((vz0/realc)-dsqrt(dum1))
      kz2=dfloat(is)*(oc/realc)*gamz0**2*((vz0/realc)+dsqrt(dum1))
  104 continue
c   calculate skin depth for rhomax at cutoff frequency (w.r.t. rwref)
      rhomax=nrhomx*1.72d-8
      delta=dsqrt(realc**2*rhomax/(2.0d0*pi*wcref*9.0d9))
c
      write(6,164) imode,im,in,is,xmn,imode,rwref
  164 format(' imode=',i2,' (im=',i2,', in=',i2,' is=',i2,
     &'), xmn=',1pd10.3,', rwref(',i2,')=',1pd11.4,' cm')
      write(6,15) imode,bz0bzg,imode,imode,imode,bzg(imode)
   15 format('   bz0/bzg(',i2,')=',1pd12.5/
     &'   bzg(',i2,')=grazing mag. field of mode',i2,' w.r.t. rwref(',
     &i2,')=',1pd12.5,' gauss')
      write(6,16) imode,imode,imode,fcref(imode),imode,imode,
     &f1,f2,kz1,kz2
   16 format('   fcref(',i2,')=cut-off freq. of mode',i2,
     &' w.r.t. rwref(',i2,')=',1pd12.5,' hz'/3x,
     &'points of intersections in freq. vs kz diagram of mode',
     &i2,' w.r.t. rwref(',i2,'):'
     &/'      (f1, f2)=',1pd11.4,',',1pd11.4,' hz'/'    (kz1, kz2)=',
     &1pd11.4,',',1pd11.4,' /cm')
      if(bz0bzg.lt.1.0d0) write(6,166) imode
  166 format('   ** note: bz0/bzg(',i2,') less than 1.0d0 **')
      write(6,14) imode,nrcavg,imode,ndrc,imode,nrc1,imode,nrc2
   14 format('   rcavg/rwref(',i2,')=',1pd10.3,', drc/rwref(',i2,
     &')=',1pd10.3/'   rc1/rwref(',i2,')=',1pd10.3,', rc2/rwref(',
     &i2,')=',1pd10.3)
      write(6,169) opnrc1,opnrc2,opnrc3
  169 format('   nrc(1st optimum)=',1pd10.3/
     &'   nrc(2nd optimum)=',1pd10.3/
     &'   nrc(3rd optimum)=',1pd10.3)
      write(6,168) nrhomx,fcref(imode),delta
  168 format('   skin depth (for rho/rhocu=',1pd9.2,' and f=',1pd9.2,
     &' hz) =',1pd10.3,' cm'/)
  163 continue
c     
      do 10 i1=1,101
      do 10 i2=1,101
      do 10 imode=1,imdmax
      apfwd2(i1,i2,imode)=-1.0d50
   10 continue
      prefmn=1.0d50
      prefmx=-1.0d50
ccc specify beam axial velocity spread tz (=delpz/pz0)
c   (always define itzmax, start itz from 1, and end itz at itzmax)
c   (to search for the saturation pref for the maximum saturated output 
c    power with a beam velocity spread, it is much faster to search for
c    the saturation pref for the cold beam case first.)
c   (this program will take the saturation pref obtained for case "itz"
c    as the starting pref for the search of the saturation pref for case
c    "itz+1".)
      itzmax=1
      do 200 itz=1,itzmax
      tz=0.05d0*(itz-1)
      atz(itz)=tz
c
ccc specify no of electrons
c
      itmax=9
      ipzmax=4*idint(tz*100.0d0)+1
      ircmax=1
      ipsimx=1
      iphimx=11
c note: if all the modes have im=0, set ipsimx to 1.
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
      write(6,17) tz,p0,pp0,pz0,delpz,ppmin,ppmax,pzmin,pzmax,
     &v0,vp0,vz0,vpmin,vpmax,vzmin,vzmax
   17 format(' momentum and velocity are in gaussian units:'/
     &' <tz=delpz/pz0=',1pd10.3,'>'/' p0=',1pd10.3,
     &', pp0=',1pd10.3,', pz0=',1pd10.3,', delpz=',1pd10.3/
     &' ppmin=',1pd10.3,', ppmax=',1pd10.3,/
     &' pzmin=',1pd10.3,', pzmax=',1pd10.3/
     &' v0=',1pd10.3,', vp0=',1pd10.3,', vz0=',1pd10.3/
     &' vpmin=',1pd10.3,', vpmax=',1pd10.3/
     &' vzmin=',1pd10.3,', vzmax=',1pd10.3/)
      write(6,253) imdmax,imdrc,imdbz,imdg,imdp,imdbcl,
     &imdbcr,itera,icut
  253 format(' imdmax=',i2,', imdrc=',i2,', imdbz=',i2,
     &', imdg=',i2,', imdp=',i2,', imdbcl=',i2,', imdbcr=',i2/
     &' itera=',i2,', icut=',i2)
      write(6,19) itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
   19 format(' itmax=',i3,', ipzmax=',i3,', ircmax=',i3,
     &', ipsimx=',i3,', iphimx=',i3,', itotal=',i7/)
      write(6,37)
   37 format('**********************************************************
     &**********************')

      do 30 i1=1,101
      do 30 i2=1,101
      do 30 imode=1,imdmax
      apfwd1(i1,i2,imode)=0.0d0
      agfwd1(i1,i2,imode)=0.0d0
      aefwd1(i1,i2,imode)=0.0d0
   30 continue

ccc specify incident wave frequency.  
c   (always define iwmax, start iw from 1, and end iw at iwmax)
      iwmax=1
      do 100 iw=1,iwmax
      fhzmin=16.8d9
      fhzmax=17.2d9
      if(iwmax.gt.1) fstep=(fhzmax-fhzmin)/dfloat(iwmax-1)
      if(iwmax.eq.1) fstep=0.0d0
      fhz=fhzmin+fstep*dfloat(iw-1)
      if(iwmax.eq.1) fhz=16.925d9
      w=2.0d0*pi*fhz
c 
      aw(1)=w
      aw(2)=aw(1)*2.0d0
      afreq(iw,1)=aw(1)/(2.0*pi)
      afreq(iw,2)=aw(2)/(2.0*pi)
c
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
      phsin(1)=0.0d0+(iphase-1)*2.0d0*pi/dfloat(iphsmx)
      phsin(2)=0.0d0
c
      ipinmx=31
      ipinmx=1
      do 100 ipin=1,ipinmx
      pin=1.0d0*1.20d0**(ipin-1)
      if(itz.gt.1) pin=apin2(iw,itz-1,imdg)*1.2**(ipin-3)
      if(ipin.eq.1) pin=1.0d-3
      if(ipinmx.eq.1.and.itz.eq.1) pin=7.0d0
      if(ipinmx.eq.1.and.itz.eq.2) pin=14.0d0
c 
      apin(1)=pin
      apin(2)=0.0d0
c
      pref=apin(imdg)
      apin1(ipin,iw)=pref

c   the following statement will terminate do loop ipin when the calculated
c   output power is lower than the maximum power by a certain amount (specified
c   near the end of do loop 100). cancellation of the following statement
c   will allow completion of do loop ipin.
      if(iskip.eq.1) go to 100
c
      if(pref.gt.prefmx) prefmx=pref
      if(pref.lt.prefmn) prefmn=pref
      if(ipinmx.eq.1) prefmx=2.0d0*prefmn

c   calculate the no of cyclotron orbits and maximum guide wavelengths 
c   over the total length of the rf structure.
      ncyc=fcyc*(zmark(izmark)-zmark(1))/vz0  
      nwvmax=0.0d0
      do 165 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      is=iis(imode)
      xmn=axmn(iabs(im)+1,in)
      rwref=arwref(imode)
      wcref=xmn*realc/rwref
      if(w.lt.wcref) nwave(imode)=0.0d0
      if(w.lt.wcref) lamdag(imode)=1.0d50
      if(w.lt.wcref) go to 167
      kz2=(w/realc)**2-xmn**2/rwref**2
      kz=dsqrt(kz2)
      lamdag(imode)=2.0d0*pi/kz
      nwave(imode)=(zmark(izmark)-zmark(1))/lamdag(imode)
  167 continue
      if(nwave(imode).gt.nwvmax) nwvmax=nwave(imode)
  165 continue
      nlarge=dmax1(ncyc,nwvmax)

ccc specify no of steps for z-integration. 
      do 40 idiv=15,15
      izstep=idiv*idint(nlarge+1.0d0)
c
c   for a new problem, always check convergence with respect to izstep. 
c   do loop 40 is for such a test. "idiv" is the no of steps per cyclotron 
c   orbit or per wavelength (2*pi/kz), whichever makes izstep greater.

ccc guess complex reflection coefficient  
      cguess=dcmplx(1.0d-3,1.0d-3)
      if(ipin.eq.1) cg(1)=cguess
c   cg(1) is the complex reflection coefficient (when pin>0) or amplitude
c   of the backward wave at the left end (when pin=0) of mode imdbcl, which
c   has been formulated to be the root of the equation cbc(cg)=0 of mode 
c   imdbcr. cg(1) defined above is a guessed value to be supplied to
c   subroutine muller as a starting value for the search of the correct 
c   reflection coefficient or backward wave amplitude (referred to the z=0 
c   plane) of mode imdbcl which satisfies the boundary condition cbc(cg)=0
c   for mode imdbcr at the right end
      ep1=1.0d-4
      ep2=ep1
      iroot=1
      imaxit=50
      icont=0
      idiag=0
c ********
c     call coord(arc,apsi,app,aphi,apz,atime,agamma,weight)
c     call bdiag
c ********
      if(itera.eq.1) call muller(0,iroot,cg,imaxit,ep1,ep2,cbc,.false.)
      if(itera.ne.1) cg(1)=dcmplx(0.0d0,0.0d0)
      idiag=1
      value=cdabs(cbc(cg(1)))
c*************************
c     call bdiag
c*************************
c  see comments in subroutine pdiag for the the diagnostic functions of
c  this subroutine
c     call pdiag(21)
c*************************
c   muller is called to solve equation cbc(cg)=0 for the complex root cg.
c   'idiag' is an instruction for diagnostics. if (and only if) idiag=1,
c   function cbc will diagnose and store the profiles of wall radius, resis-
c   tivity, rf field, wave powers, ohmic losses, gain, final beam coordinates,
c   and conservation/consistency factors, and pass them through the common 
c   block to the main program.
c   'icont' monitors the no of times that function cbc is called for root 
c   searching ('icont' goes up by one upon each call to cbc).
c   'value' monitors the accuracy of the root returned by muller.
      gr=cg(1)
      gi=-ci*cg(1)
      gamp=cdabs(cg(1))
      if(gr.ne.0.0d0.or.gi.ne.0.0d0) gphase=datan2(gi,gr)
      if(gr.eq.0.0d0.and.gi.eq.0.0d0) gphase=0.0d0
c   diagnose gain, efficiencies, and power conservation factor, etc.
      etab=0.0d0
      do 90 it=1,itmax
      do 90 ipz=1,ipzmax
      do 90 irc=1,ircmax
      do 90 ipsi=1,ipsimx
      do 90 iphi=1,iphimx
      gamma=agamma(it,ipz,irc,ipsi,iphi)
      etab=etab+weight(it,ipz,irc,ipsi,iphi)*(gamma0-gamma)
     &/(gamma0-1.0d0)
   90 continue
c     
      izmax=izstep+1
      do 154 imode=1,imdmax
      pohmdb(imode)=-1.0d50 
      gfwd(imode)=-1.0d50
      ratio1(imode)=1.0d50
      ratio2(imode)=1.0d50
      do 154 iz=1,izmax
      agfwd(iz,imode)=-1.0d50 
      agbwd(iz,imode)=-1.0d50 
      agohm(iz,imode)=-1.0d50 
  154 continue
      delz=(zmark(izmark)-zmark(1))/dfloat(izstep)
      sump=0.0d0
      do 151 imode=1,imdmax
      etafwd(imode)=(apfwd(izstep+1,imode)-apin(imode))/(pbkw*1.0d3)
      etabwd(imode)=(apbwd(1,imode)-apbwd(izstep+1,imode))/(pbkw*1.0d3)
      etaohm(imode)=pohmsm(imode)/(pbkw*1.0d3)
      if(apfwd(izstep+1,imode).gt.0.0d0)
     &ratio1(imode)=pohmsm(imode)/apfwd(izstep+1,imode)
      if(apfwd(izstep+1,imode).gt.0.0d0)
     &ratio2(imode)=apbwd(izstep+1,imode)/apfwd(izstep+1,imode)
      sump=sump+apnet(izstep+1,imode)-apnet(1,imode)+pohmsm(imode)
      if(apfwd(izstep+1,imode).gt.0.0d0)
     &gfwd(imode)=10.0d0*dlog10(apfwd(izstep+1,imode)/pref)
      if(pohmsm(imode).gt.0.0d0)
     &pohmdb(imode)=10.0*dlog10(pohmsm(imode)/pref)
c
      if(apfwd(1,imode).gt.0.0d0)
     &agfwd(1,imode)=10.0d0*dlog10(apfwd(1,imode)/pref)
      if(apfwd(1,imode).eq.0.0d0)
     &agfwd(1,imode)=10.0d0*dlog10(0.5d0*apfwd(2,imode)/pref)
      if(apbwd(1,imode).gt.0.0d0)
     &agbwd(1,imode)=10.0d0*dlog10(apbwd(1,imode)/pref)
      if(apbwd(1,imode).eq.0.0d0)
     &agbwd(1,imode)=10.0d0*dlog10(0.5d0*apbwd(2,imode)/pref)
      if(apohmz(1,imode).gt.0.0d0)
     &agohm(1,imode)=10.0d0*dlog10(apohmz(1,imode)/pref)
      do 80 iz=1,izstep
      if(apfwd(iz+1,imode).gt.0.0d0)
     &agfwd(iz+1,imode)=10.0d0*dlog10(apfwd(iz+1,imode)/pref)
      agfwdz(iz+1,imode)=(agfwd(iz+1,imode)-agfwd(iz,imode))/delz
      if(apbwd(iz+1,imode).gt.0.0d0)
     &agbwd(iz+1,imode)=10.0d0*dlog10(apbwd(iz+1,imode)/pref)
      if(apohmz(iz+1,imode).gt.0.0d0)
     &agohm(iz+1,imode)=10.0d0*dlog10(apohmz(iz+1,imode)/pref)
   80 continue
  151 continue
      pconsv=(sump-etab*pbkw*1.0d3)/(etab*pbkw*1.0d3)
      do 152 imode=1,imdmax
      if(iplot.ne.1) go to 105
      if(imode.ne.1) write(6,37)
      write(6,155) imode,iim(imode),iin(imode),iis(imode),
     &imdmax,imdg,imdp,imdbcl,imdbcr
  155 format(1x,' imode=',i2,' (im=',i2,',in=',i2,',is=',i2,
     &'),imdmax=',i2,', imdg=',i2,', imdp=',i2,', imdbcl=',i2,
     &', imdbcr=',i2)
c   plot rf structure shape and wall resistivity
c   find rwmax and nrhomx
      rwmax=arw(1)
      nrhomx=anrho(1)
      do 69 iz=1,izstep
      if(arw(iz+1).gt.rwmax) rwmax=arw(iz+1)
      if(anrho(iz+1).gt.nrhomx) nrhomx=anrho(iz+1)
   69 continue
      call sscale(zmark(1),zmark(izmark),0.0d0,rwmax)
      call splot(az,arw,izstep+1,'rr  rw(cm) vs z(cm)$')
      if(nrhomx.eq.0.0d0) go to 106
      call sscale(zmark(1),zmark(izmark),0.0d0,nrhomx)
      call splot(az,anrho,izstep+1,'**  rho(ohm-m)/1.72e-8 vs z(cm)$')
  106 continue
c   plot rf field and power profiles.
c   find fampmx, pnetmx, pnetmn, pfwdmx, pfwdmn, gbwdmx, gbwdmn, 
c   gohmmx, gohmmn, gfwdmx, gfwdmn, gzmax, and gzmin
      fampmx=afamp(1,imode)
      gbwdmx=agbwd(1,imode)
      gbwdmn=agbwd(1,imode)
      pnetmx=apnet(1,imode)
      pnetmn=apnet(1,imode)
      gohmmx=agohm(1,imode)
      gohmmn=agohm(1,imode)
      gfwdmx=agfwd(1,imode)
      gfwdmn=agfwd(1,imode)
      gzmax=agfwdz(1,imode)
      gzmin=agfwdz(1,imode)
      do 70 iz=1,izstep
      if(afamp(iz+1,imode).gt.fampmx) fampmx=afamp(iz+1,imode)
      if(agbwd(iz+1,imode).gt.gbwdmx) gbwdmx=agbwd(iz+1,imode)
      if(agbwd(iz+1,imode).le.gbwdmn) gbwdmn=agbwd(iz+1,imode)
      if(apnet(iz+1,imode).gt.pnetmx) pnetmx=apnet(iz+1,imode)
      if(apnet(iz+1,imode).le.pnetmn) pnetmn=apnet(iz+1,imode)
      if(agohm(iz+1,imode).gt.gohmmx) gohmmx=agohm(iz+1,imode)
      if(agohm(iz+1,imode).le.gohmmn) gohmmn=agohm(iz+1,imode)
      if(agfwd(iz+1,imode).gt.gfwdmx) gfwdmx=agfwd(iz+1,imode)
      if(agfwd(iz+1,imode).le.gfwdmn) gfwdmn=agfwd(iz+1,imode)
      if(dabs(agfwdz(iz+1,imode)).gt.10.0d0) go to 70
      if(agfwdz(iz+1,imode).gt.gzmax) gzmax=agfwdz(iz+1,imode)
      if(agfwdz(iz+1,imode).le.gzmin) gzmin=agfwdz(iz+1,imode)
   70 continue
      if(pnetmn.ge.0.0d0) pnetmn=0.0d0
c     if(gfwdmn.ge.0.0d0) gfwdmn=0.0d0
c     if(gzmin.ge.0.0d0) gzmin=0.0d0
c     if(gbwdmn.ge.0.0d0) gbwdmn=0.0d0
c     if(gohmmn.ge.0.0d0) gohmmn=0.0d0
      if(gfwdmn.lt.(gfwdmx-50.0d0)) gfwdmn=gfwdmx-50.0d0
      if(gbwdmn.lt.(gbwdmx-50.0d0)) gbwdmn=gbwdmx-50.0d0
      if(gohmmn.lt.(gohmmx-50.0d0)) gohmmn=gohmmx-50.0d0

      call sscale(zmark(1),zmark(izmark),0.0d0,fampmx)
      call splot(az,afamp(1,imode),izstep+1,
     &'ff  amplitude of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),-pi,pi)
      call splot(az,afphse(1,imode),izstep+1,
     &'pp  phase(radian) of f vs z(cm)$')
      call sscale(zmark(1),zmark(izmark),pnetmn,pnetmx)
      call splot(az,apnet(1,imode),izstep+1, 
     &'nn  net wave power(w) vs z(cm)$') 
      write(6,21)
   21 format(/ ' note: forward wave power (hence the plot below) may be
     &inaccurately diagnosed'/'       in regions where w is slightly abo
     &ve cutoff or where the wall is tapered.'/
     &'       it was artificially set to 1.0d-10 watt where w < wcmn.')
      call sscale(zmark(1),zmark(izmark),gfwdmn,gfwdmx)
      call splot(az,agfwd(1,imode),izstep+1,
     &'gg  gain=10*log(forward wave power/pref) vs z(cm)$') 
      write(6,21)
      call sscale(zmark(1),zmark(izmark),gzmin,gzmax)
      call splot(az,agfwdz(1,imode),izstep+1,
     &'gg  forward wave growth rate(db/cm) vs z(cm)$')
      if(pohmsm(imode).eq.0.0d0) go to 107
      call sscale(zmark(1),zmark(izmark),gohmmn,gohmmx)
      call splot(az,agohm(1,imode),izstep+1,
     &'oo  10*log(ohmic loss per cm/pref) vs z(cm)$')
  107 continue
      write(6,22)
   22 format(/ ' note: backward wave power (hence the plot below) may be
     & inaccurately diagnosed'/'       in regions where w is slightly ab
     &ove cutoff or where the wall is tapered.'/
     &'       it was artificially set to 1.0d-10 watt where w < wcmn.')
      call sscale(zmark(1),zmark(izmark),gbwdmn,gbwdmx)
      call splot(az,agbwd(1,imode),izstep+1, 
     &'bb  10*log(backward wave power/pref) vs z(cm)$') 
      write(6,23)
   23 format(/)
  105 continue

      if(iwrite.ne.1) go to 108
c   write results
      write(6,155) imode,iim(imode),iin(imode),iis(imode),
     &imdmax,imdg,imdp,imdbcl,imdbcr
      write(6,26) afreq(iw,imode),apin(imode),imdg,pref
   26 format(' freq=',1pd12.5,' hz, pin(of this mode)=',1pd10.3,
     &' w, pref=pin(',i2,')=',1pd10.3,' w')
      write(6,25) ncyc,nwave(imode),lamdag(imode)
   25 format(' ncyc=',1pd9.2,', nwave=',1pd9.2,', lamdag=',1pd9.2,
     &' cm') 
      write(6,98) idiv,izstep,itotal,itera,ep1   
   98 format(' idiv=',i3,', izstep=',i7,', itotal=',i7,
     &', itera=',i2,', ep1=',1pd8.1)
      write(6,24) imdbcl,gamp,gphase,icont,value
   24 format(' gamma for mode',i3,': gamp=',1pd11.4,
     &', gphase(radian)=',1pd10.3,'(',i3,',',1pd8.1,')')
      write(6,97) phsin(imode),phsout(imode),delphs(imode)
   97 format(' phsin(radian)=',1pd10.3,', phsout(radian)=',1pd10.3,
     &', delphs(radian)=',1pd10.3)
      write(6,27) apnet(1,imode),apnet(izstep+1,imode),apfwd(1,imode),
     &apfwd(izstep+1,imode),apbwd(1,imode),apbwd(izstep+1,imode),
     &pohmsm(imode),pconsv
   27 format(' pnet (net wave power):     at left end=',1pd10.3,
     &' w, at right end=',1pd10.3,' w'/ 
     &' pfwd (forward wave power): at left end=',1pd10.3,
     &' w, at right end=',1pd10.3,' w'/
     &' pbwd (backward wave power):at left end=',1pd10.3,
     &' w, at right end=',1pd10.3,' w!'/
     &' pohmsm (sum of wall losses)=',1pd10.3,
     &' w,        (power conservation=',1pd9.2,')')
      write(6,28) gfwd(imode),gmcsv
   28 format(' gfwd=10*log(pfwd at right end/pref)= ',1pd10.3,' db',
     &'(gamma consistency=',1pd9.2,')')
      write(6,34) agbwd(1,imode),ptcsv
   34 format(' gbwd=10*log(pbwd at left end/pref)=',1pd10.3,
     &' db(ptheta conservation=',1pd9.2,')')
      write(6,35) ratio1(imode),ratio2(imode)
   35 format(' pohmsm/pfwd(right end)=',1pd10.3,
     &', pbwd(right end)/pfwd(right end)=',1pd10.3,'!')
      write(6,36) etafwd(imode),etab,etabwd(imode),etaohm(imode)
   36 format(' etafwd=',1pd11.4,', etab=',1pd11.4,
     &'  (etabwd=',1pd11.4,', etaohm=',1pd11.4,')')
      if(imode.ne.imdmax) write(6,252)
      if(imode.eq.imdmax) write(6,37)
  108 continue
  152 continue
   40 continue
c   for a fixed tz (beam velocity spread), store pfwd (forward wave power at
c   right end), gfwd (=10*log10(pfwd/pref)), etafwd (=(pfwd-pin)/beam power)
c   as functions of the specified pref (reference power) for different 
c   frequencies as specified. 
c   store gfwd at the minimum specified pref as a function of frequency for
c   different tz as specified.
c   find and store the maximized pfwd of mode imdp (maximized within the 
c   specified pref range) at the output end and the associated pref, gfwd,
c   etafwd, etabwd (=(pbwd at left end - pbwd at right end)/beam power),
c   and etaohm (total ohmic power/beam power), etc. of all modes as functions
c   of the specified frequency for different tz as specified.
      do 153 imode=1,imdmax
      apfwd1(ipin,iw,imode)=apfwd(izstep+1,imode)
      agfwd1(ipin,iw,imode)=agfwd(izstep+1,imode)
      aefwd1(ipin,iw,imode)=etafwd(imode)
      if(ipin.eq.1.and.apfwd(izstep+1,imode).gt.0.0d0)
     &agfwd0(iw,itz,imode)=10.0*dlog10(apfwd(izstep+1,imode)/pref)
      if(ipin.eq.1.and.apfwd(izstep+1,imode).le.0.0d0)
     &agfwd0(iw,itz,imode)=-1.0d50
      if(ipin.eq.1)
     &agbwd0(iw,itz,imode)=agbwd(1,imode)
  153 continue
      do 158 imode=1,imdmax
      if(apfwd1(ipin,iw,imdp).lt.(0.85d0*apfwd2(iw,itz,imdp))) iskip=1
      if(apfwd1(ipin,iw,imdp).lt.apfwd2(iw,itz,imdp)) go to 158
      apin2(iw,itz,imode)=apin(imode)
      agfwd2(iw,itz,imode)=agfwd1(ipin,iw,imode)
      aefwd2(iw,itz,imode)=aefwd1(ipin,iw,imode)
      aebwd2(iw,itz,imode)=etabwd(imode)
      aeohm2(iw,itz,imode)=etaohm(imode)
      agbwd2(iw,itz,imode)=agbwd(1,imode)
      icont2(iw,itz,imode)=icont
      avalu2(iw,itz,imode)=value
      apcnsv(iw,itz)=pconsv
      aetab(iw,itz)=etab
      apfwd2(iw,itz,imode)=apfwd1(ipin,iw,imode)
  158 continue
  100 continue
      if(ipinmx.eq.1) go to 110
c   plot pfwd, gfwd, and efwd at the output end as functions of pref for
c   fixed tz and different frequencies as specified.
c
      do 157 imode=1,imdmax
c   find pfwdmx, gfwdmx, gfwdmn, efwdmx, and efwdmn
      pfwdmx=-1.0d50
      pfwdmn=1.0d50
      gfwdmx=-1.0d50
      gfwdmn=1.0d50
      efwdmx=-1.0d50
      efwdmn=1.0d50
      do 57 iw=1,iwmax
      do 57 ipin=1,ipinmx
      if(apfwd1(ipin,iw,imode).gt.pfwdmx) pfwdmx=apfwd1(ipin,iw,imode)
      if(apfwd1(ipin,iw,imode).le.pfwdmn) pfwdmn=apfwd1(ipin,iw,imode)
      if(agfwd1(ipin,iw,imode).gt.gfwdmx) gfwdmx=agfwd1(ipin,iw,imode)
      if(agfwd1(ipin,iw,imode).le.gfwdmn) gfwdmn=agfwd1(ipin,iw,imode)
      if(aefwd1(ipin,iw,imode).gt.efwdmx) efwdmx=aefwd1(ipin,iw,imode)
      if(aefwd1(ipin,iw,imode).le.efwdmn) efwdmn=aefwd1(ipin,iw,imode)
   57 continue
      if(pfwdmn.eq.pfwdmx) pfwdmx=2.0d0*dabs(pfwdmn)
      if(gfwdmn.eq.gfwdmx) gfwdmx=2.0d0*dabs(gfwdmn)
      if(efwdmn.eq.efwdmx) efwdmx=2.0d0*dabs(efwdmn)
      write(6,155) imode,iim(imode),iin(imode),iis(imode),
     &imdmax,imdg,imdp,imdbcl,imdbcr
      call sscale(prefmn,prefmx,pfwdmn,pfwdmx)
      imax=iwmax-1
      if(iwmax.eq.1) go to 52
      do 51 i=1,imax
   51 call splot(apin1(1,i),apfwd1(1,i,imode),ipinmx,'*   $')
   52 call splot(apin1(1,iwmax),apfwd1(1,iwmax,imode),ipinmx,
     &'pp  pfwd(w) vs pref(w)$')
      call sscale(prefmn,prefmx,gfwdmn,gfwdmx)
      if(iwmax.eq.1) go to 54
      do 53 i=1,imax
   53 call splot(apin1(1,i),agfwd1(1,i,imode),ipinmx,'*   $')
   54 call splot(apin1(1,iwmax),agfwd1(1,iwmax,imode),ipinmx,
     &'gg  gfwd(db) vs pref(w)$')
      call sscale(prefmn,prefmx,efwdmn,efwdmx)
      if(iwmax.eq.1) go to 56
      do 55 i=1,imax
   55 call splot(apin1(1,i),aefwd1(1,i,imode),ipinmx,'*   $')
   56 call splot(apin1(1,iwmax),aefwd1(1,iwmax,imode),ipinmx,
     &'ee  etafwd vs pref(w)$')
      write(6,37)
  157 continue
      do 109 imode=1,imdmax
      if(itz.eq.itzmax) go to 109
      write(6,38)
   38 format(15x,'-- intermediate summary for do loop 200 --'/)
c   write pfwd at right end (maximized within the range of specified pref)
c   and the associated pref, gfwd, etafwd, etabwd, etaohm as functions of
c   the specified frequency for a fixed tz. numbers in parentheses are
c   gfwd for the minimum specified pref (which is the small signal gain if
c   prefmn is sufficiently small).
      write(6,155) imode,iim(imode),iin(imode),iis(imode),
     &imdmax,imdg,imdp,imdbcl,imdbcr
      write(6,39) atz(itz),itot(itz),itera,idiv,ep1
   39 format('tz=',1pd9.2,', itotal=',i7,', itera=',i2,
     &', idiv=',i4,', ep1=',1pd8.1/)
      write(6,41)
   41 format(' freq(hz):','  pfwd(w)','    etafwd ','  gfwd(db)',
     &' (at prefmn)  ',' gbwd(db)',' (at prefmn)')
      do 91 iw=1,iwmax
      write(6,42) afreq(iw,imode),apfwd2(iw,itz,imode),
     &aefwd2(iw,itz,imode),agfwd2(iw,itz,imode),agfwd0(iw,itz,imode),
     &agbwd2(iw,itz,imode),agbwd0(iw,itz,imode)
   42 format(1pd10.4,':',1pd9.2,' ',1pd9.2,' ',1pd9.2,' (',1pd9.2,
     &')  ',1pd9.2,' (',1pd9.2,')')
   91 continue
      write(6,43)
   43 format(/' freq(hz):','   pin(w) ','   etabwd ','   etaohm',
     &'     etab ','   icont','   value','     pconsv')
      do 92 iw=1,iwmax
      write(6,44) afreq(iw,imode),apin2(iw,itz,imode),
     &aebwd2(iw,itz,imode),aeohm2(iw,itz,imode),aetab(iw,itz),
     &icont2(iw,itz,imode),avalu2(iw,itz,imode),apcnsv(iw,itz)
   44 format(1pd10.4,':',1pd9.2,' ',1pd9.2,' ',1pd9.2,
     &' ',1pd10.3,' ',i4,'  ',1pd9.2,' ',1pd9.2)
   92 continue
      write(6,37)
  109 continue
  110 continue
  200 continue
      do 111 imode=1,imdmax
      if(iwmax.eq.1) go to 111
c   plot gfwd for the minimum specified pref (which is the small signal
c   gain if prefmn is sufficiently amall) as a function of frequency for 
c   different tz as specified.
c   plot pfwd at the output end (maximized w.r.t. pref in the specified pref 
c   range) and the associated pref, gfwd, etafwd, etabwd, and etaohm as
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
      if(apfwd2(iw,itz,imode).gt.pfwdmx) pfwdmx=apfwd2(iw,itz,imode)
      if(apfwd2(iw,itz,imode).le.pfwdmn) pfwdmn=apfwd2(iw,itz,imode)
      if(agfwd0(iw,itz,imode).gt.gfwdmx) gfwdmx=agfwd0(iw,itz,imode)
      if(agfwd2(iw,itz,imode).gt.gfwdmx) gfwdmx=agfwd2(iw,itz,imode)
      if(agfwd0(iw,itz,imode).le.gfwdmn) gfwdmn=agfwd0(iw,itz,imode)
      if(agfwd2(iw,itz,imode).le.gfwdmn) gfwdmn=agfwd2(iw,itz,imode)
      if(aefwd2(iw,itz,imode).gt.efwdmx) efwdmx=aefwd2(iw,itz,imode)
      if(aefwd2(iw,itz,imode).le.efwdmn) efwdmn=aefwd2(iw,itz,imode)
      if(aebwd2(iw,itz,imode).gt.ebwdmx) ebwdmx=aebwd2(iw,itz,imode)
      if(aebwd2(iw,itz,imode).le.ebwdmn) ebwdmn=aebwd2(iw,itz,imode)
      if(agbwd2(iw,itz,imode).gt.gbwdmx) gbwdmx=agbwd2(iw,itz,imode)
      if(agbwd0(iw,itz,imode).gt.gbwdmx) gbwdmx=agbwd0(iw,itz,imode)
      if(agbwd2(iw,itz,imode).le.gbwdmn) gbwdmn=agbwd2(iw,itz,imode)
      if(agbwd0(iw,itz,imode).le.gbwdmn) gbwdmn=agbwd0(iw,itz,imode)
      if(aeohm2(iw,itz,imode).gt.eohmmx) eohmmx=aeohm2(iw,itz,imode)
      if(aeohm2(iw,itz,imode).le.eohmmn) eohmmn=aeohm2(iw,itz,imode)
   59 continue
      if(pfwdmn.ge.0.0) pfwdmn=0.0
c     if(gfwdmn.ge.0.0d0) gfwdmn=0.0d0
c     if(gbwdmn.ge.0.0d0) gbwdmn=0.0d0
c     if(efwdmn.ge.0.0d0) efwdmn=0.0d0
c     if(ebwdmn.ge.0.0d0) ebwdmn=0.0d0
c     if(eohmmn.ge.0.0d0) eohmmn=0.0d0
      write(6,155) imode,iim(imode),iin(imode),iis(imode),
     &imdmax,imdg,imdp,imdbcl,imdbcr
      imax=itzmax-1
      call sscale(afreq(1,imode),afreq(iwmax,imode),gfwdmn,gfwdmx)
      if(itzmax.eq.1) go to 62
      do 61 i=1,imax
   61 call splot(afreq(1,imode),agfwd0(1,i,imode),iwmax,'*   $')
   62 call splot(afreq(1,imode),agfwd0(1,itzmax,imode),iwmax,
     &'ss  gfwd(in db, at pref=prefmn) vs freq(hz)$')
      call sscale(afreq(1,imode),afreq(iwmax,imode),pfwdmn,pfwdmx)
      if(itzmax.eq.1) go to 64
      do 63 i=1,imax
   63 call splot(afreq(1,imode),apfwd2(1,i,imode),iwmax,'*   $')
   64 call splot(afreq(1,imode),apfwd2(1,itzmax,imode),iwmax,
     &'pp  pfwd(in w, maximized in specified range of pref) vs freq(hz)
     &$')
      call sscale(afreq(1,imode),afreq(iwmax,imode),prefmn,prefmx)
      if(itzmax.eq.1) go to 66
      do 65 i=1,imax
   65 call splot(afreq(1,imode),apin2(1,i,imode),iwmax,'*   $')
   66 call splot(afreq(1,imode),apin2(1,itzmax,imode),iwmax,
     &'dd  pref(in w, for above pfwd) vs freq(hz)$')
      call sscale(afreq(1,imode),afreq(iwmax,imode),gfwdmn,gfwdmx)
      if(itzmax.eq.1) go to 68
      do 67 i=1,imax
   67 call splot(afreq(1,imode),agfwd2(1,i,imode),iwmax,'*   $')
   68 call splot(afreq(1,imode),agfwd2(1,itzmax,imode),iwmax,
     &'gg  gfwd(in db, for above pfwd) vs freq(hz)$')
      call sscale(afreq(1,imode),afreq(iwmax,imode),efwdmn,efwdmx)
      if(itzmax.eq.1) go to 72
      do 71 i=1,imax
   71 call splot(afreq(1,imode),aefwd2(1,i,imode),iwmax,'*   $')
   72 call splot(afreq(1,imode),aefwd2(1,itzmax,imode),iwmax,
     &'ee  etafwd(for above pfwd) vs freq(hz)$')
      if(eohmmx.le.0.0d0) go to 112
      call sscale(afreq(1,imode),afreq(iwmax,imode),eohmmn,eohmmx)
      if(itzmax.eq.1) go to 74
      do 73 i=1,imax
   73 call splot(afreq(1,imode),aeohm2(1,i,imode),iwmax,'*   $')
   74 call splot(afreq(1,imode),aeohm2(1,itzmax,imode),iwmax,
     &'oo  etaohm(for above pfwd) vs freq(hz)$')
  112 continue
      call sscale(afreq(1,imode),afreq(iwmax,imode),ebwdmn,ebwdmx)
      if(itzmax.eq.1) go to 76
      do 75 i=1,imax
   75 call splot(afreq(1,imode),aebwd2(1,i,imode),iwmax,'*   $')
   76 call splot(afreq(1,imode),aebwd2(1,itzmax,imode),iwmax,
     &'bb  etabwd(for above pfwd) vs freq(hz)$')
      call sscale(afreq(1,imode),afreq(iwmax,imode),gbwdmn,gbwdmx)
      if(itzmax.eq.1) go to 78
      do 77 i=1,imax
   77 call splot(afreq(1,imode),agbwd2(1,i,imode),iwmax,'*   $')
   78 call splot(afreq(1,imode),agbwd2(1,itzmax,imode),iwmax,
     &'bb  10*log(pbwd at left end(for above pfwd)/pref) vs freq(hz)$')
      call sscale(afreq(1,imode),afreq(iwmax,imode),gbwdmn,gbwdmx)
      if(itzmax.eq.1) go to 82
      do 81 i=1,imax
   81 call splot(afreq(1,imode),agbwd0(1,i,imode),iwmax,'*   $')
   82 call splot(afreq(1,imode),agbwd0(1,itzmax,imode),iwmax,
     &'ss  10*log(pbwd at left end(at pref=prefmn)/pref) vs freq(hz)$')
      write(6,37)
  111 continue
      do 159 imode=1,imdmax
c              -- final summary for do loop 200 --
c   write pfwd at right end (at a pref for maximum pfwd of the ref mode imdref)
c   and the associated pin, gfwd, etafwd, etabwd, etaohm as functions of
c   the specified frequency for different tz. numbers in parentheses are
c   gfwd for the minimum specified pref (which is the small signal gain if
c   prefmn is sufficiently small).
      do 93 itz=1,itzmax
      write(6,45)
   45 format(20x,'-- final summary for do loop 200 --'/)
      write(6,155) imode,iim(imode),iin(imode),iis(imode),
     &imdmax,imdg,imdp,imdbcl,imdbcr
      write(6,39) atz(itz),itot(itz),itera,idiv,ep1
      write(6,41)
      do 94 iw=1,iwmax
      write(6,42) afreq(iw,imode),apfwd2(iw,itz,imode),
     &aefwd2(iw,itz,imode),agfwd2(iw,itz,imode),agfwd0(iw,itz,imode),
     &agbwd2(iw,itz,imode),agbwd0(iw,itz,imode)
   94 continue
      write(6,43)
      do 95 iw=1,iwmax
      write(6,44) afreq(iw,imode),apin2(iw,itz,imode),
     &aebwd2(iw,itz,imode),aeohm2(iw,itz,imode),aetab(iw,itz),
     &icont2(iw,itz,imode),avalu2(iw,itz,imode),apcnsv(iw,itz)
   95 continue
      if(itz.eq.itzmax) write(6,37)
      if(itz.ne.itzmax) write(6,252)
   93 continue
  159 continue
 1000 continue
      stop
      end
c ******************************************************************
      function cbc(cg)
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension y(500000),dy(500000),q(500000)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension az(20001),arw(20001),anrho(20001)
      dimension afamp(20001,3),afphse(20001,3),apnet(20001,3),
     &apfwd(20001,3),apbwd(20001,3),apohmz(20001,3)
      dimension arc(21,41,1,1,23),apsi(21,41,1,1,23),
     &app(21,41,1,1,23),aphi(21,41,1,1,23),apz(21,41,1,1,23),
     &atime(21,41,1,1,23),agamma(21,41,1,1,23),
     &weight(21,41,1,1,23)
      dimension apti(21,41,1,1,23),aptf(21,41,1,1,23),
     &aptcsv(21,41,1,1,23),agmcsv(21,41,1,1,23)
      dimension axmn(9,8),aw(11),apin(11),iis(11),iim(11),iin(11)        
      dimension phsin(11),phsout(11),delphs(11)
      dimension pohmsm(11)
      dimension adphi(20001,301),apx(301,20001),
     &apy(301,20001),azz(20001)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode1/axmn,aw,apin,iis,iim,iin,imdmax,imdbcl,imdbcr
      common/mode2/w,pin,xmn,kcapmn,is,im
      common/phase/phsin,phsout,delphs
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      common/cutoff/icut
      common/diagf/az,arw,anrho,afamp,afphse,apnet,apfwd,apbwd,apohmz,
     &pohmsm,idiag,icont
      common/diagb/arc,apsi,app,aphi,apz,atime,agamma
      common/diagp/adphi,apx,apy,azz,izmax
      common/check/agmcsv,apti,aptf,aptcsv,gmcsv,ptcsv
      external difeq
      z1=zmark(1)
      rw=radius(z1)
c check whether modes are below cut-off of the first section.
      do 70 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      if((w/realc)**2.le.xmn**2/rw**2) fhz=w/(2.0d0*pi)
      if((w/realc)**2.le.xmn**2/rw**2) write(6,1) fhz,im,in
      if((w/realc)**2.le.xmn**2/rw**2) stop 100
    1 format(' **f='1pd10.3,' hz is below cut-off of te(',
     &i2,',',i2,') mode in the first section**')
   70 continue
c
      ie=7*itotal
      ieqlst=ie+5*imdmax+1
c initialize coordinate part of coordinate/field array y at left end 
c according to the initial electron beam distribution in real and momentum 
c spaces (electron injection times(atime) are evenly spaced over the wave
c period of the mode with the lowest frequency).
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
c note: when iterative procedure is selected (itera=1), the reflection
c	coefficient cg (if pin>0) or backward wave amplitude cg (if pin=0) 
c       of mode imdbcl at the left end is evaluated to satisfy the boundary
c	condition that cbc(cg) of mode imdbcr vanishes at the right end.
c       all modes except mode imdbcl are assumed to be pure forward waves
c       at the left end, while the outgoing-wave or evanescent-wave boundary
c	condition (cbc(cg)=0) is insured at the right end for mode imdbcr only
      do 80 imode=1,imdmax
      w=aw(imode)
      pin=apin(imode)
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      kcapmn=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
      ckz2=(w/realc)**2-xmn**2*closs(z1)/rw**2
      ckz=cdsqrt(ckz2)
      cdum1=cdexp(ci*ckz*z1)
      cdum2=cdexp(-ci*ckz*z1)
      cf0=dcmplx(1.0,0.0)*cdexp(ci*phsin(imode))
      do 20 i=1,2
      if(imode.eq.imdbcl.and.pin.ne.0.0d0) cdum3=cf0*(cdum1+cg*cdum2)
      if(imode.eq.imdbcl.and.pin.ne.0.0d0) cdum4=cf0*(cdum1-cg*cdum2)
      if(imode.eq.imdbcl.and.pin.eq.0.0d0) cdum3=cf0*cdum1+cg*cdum2
      if(imode.eq.imdbcl.and.pin.eq.0.0d0) cdum4=cf0*cdum1-cg*cdum2
      if(imode.ne.imdbcl) cdum3=cf0*cdum1
      if(imode.ne.imdbcl) cdum4=cf0*cdum1
      y(ie+5*(imode-1)+1)=cdum3
      y(ie+5*(imode-1)+2)=-ci*cdum3
      y(ie+5*(imode-1)+3)=ci*ckz*cdum4
      y(ie+5*(imode-1)+4)=ckz*cdum4
      y(ie+5*(imode-1)+5)=0.0d0
      cf=dcmplx(y(ie+5*(imode-1)+1),y(ie+5*(imode-1)+2))
      cfp=dcmplx(y(ie+5*(imode-1)+3),y(ie+5*(imode-1)+4))
      kzr=ckz
      cdum10=cdexp(ci*kzr*z1)
      cdum20=cdexp(-ci*kzr*z1)
      cdum5=ci*kzr*cf
      cdum6=2.0d0*ci*kzr
      dum17=w*xmn**2*kcapmn/8.0d7
      delta=dsqrt(realc**2*rho(z1)/(2.0d0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0d0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      cffwd=(cdum5+cfp)*cdum20/cdum6
      cfbwd=(cdum5-cfp)*cdum10/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z1)
c       the backward component of cf is cfbwd*cdexp(-ci*kzr*z1)
      cfstar=dcmplx(y(ie+5*(imode-1)+1),-y(ie+5*(imode-1)+2))
      cfpstr=dcmplx(y(ie+5*(imode-1)+3),-y(ie+5*(imode-1)+4))
      cdum30=cf*cfpstr-cfstar*cfp
      az(1)=z1
      arw(1)=radius(z1)
      anrho(1)=rho(z1)/1.72d-8
      afamp(1,imode)=
     &dsqrt(y(ie+5*(imode-1)+1)**2+y(ie+5*(imode-1)+2)**2)
      if(y(ie+5*(imode-1)+2).ne.0.0d0.or. 
     &y(ie+5*(imode-1)+1).ne.0.0d0) 
     &afphse(1,imode)=datan2(y(ie+5*(imode-1)+2),y(ie+5*(imode-1)+1))         
      if(y(ie+5*(imode-1)+2).eq.0.0d0.and. 
     &y(ie+5*(imode-1)+1).eq.0.0d0) afphse(1,imode)=0.0d0
      apfwd(1,imode)=kzr*dum17*cdabs(cffwd)**2
      apbwd(1,imode)=kzr*dum17*cdabs(cfbwd)**2
      apnet(1,imode)=0.5d0*ci*dum17*cdum30
      apohmz(1,imode)=dum17*dum18*dum19*
     &(y(ie+5*(imode-1)+1)**2+y(ie+5*(imode-1)+2)**2)
      if(i.eq.1) cf0=cf0*dsqrt(pin/apfwd(1,imode))
   20 continue
   80 continue
      y(ie+5*imdmax+1)=z1

c store ptheta of all electrons at left end (when all modes have im=0).
      imode=1
   52 if(iim(imode).ne.0) go to 51
      imode=imode+1
      if(imode.le.imdmax) go to 52
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
      call bfield(z1,bz,bzp)
      oe=e*bz/(me*realc)
      vp=pp/(gamma*me)
      oc=oe/gamma
      rl=vp/oc
      r=dsqrt(rc**2+rl**2+2.0d0*rc*rl*dsin(phi-psi))
      ptheta=0.5d0*me*oe*(rl**2-rc**2)
      do 53 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      kmn=xmn/radius(z1)
      cf=dcmplx(y(ie+5*(imode-1)+1),y(ie+5*(imode-1)+2))
      ptheta=ptheta
     &-(e/realc)*r*kmn*cf*bj(1,kmn*r)*cdexp(-ci*w*time)
   53 continue
      apti(it,ipz,irc,ipsi,iphi)=ptheta
      i=i+7
   50 continue
   51 continue
       
c integrate differential equations for coordinate/field array y from left 
c end to right end.
      do 10 i=1,ieqlst
   10 q(i)=0.0d0
      l=zmark(izmark)-z1
      delz=l/dfloat(izstep)
      do 100 iz=1,izstep
c advance coordinate/field array y by one step.
      call rkint(difeq,y,dy,q,1,ieqlst,delz)
      z=y(ie+5*imdmax+1)
      rw=radius(z)
c
      if(icut.ne.1) go to 124
      do 123 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      if((w/realc)**2.gt.(xmn/rw)**2) go to 123
      y(ie+5*(imode-1)+1)=0.0d0
      y(ie+5*(imode-1)+2)=0.0d0
      y(ie+5*(imode-1)+3)=0.0d0
      y(ie+5*(imode-1)+4)=0.0d0
  123 continue
  124 continue
c
      if(idiag.ne.1) go to 99
c store radius, resistivity, rf field, and powers as functions of z
c (wave power in w, wall loss per unit length in w/cm). 
c forward wave power (apfwd) and backward wave power (apbwd) are 
c artificially set to 1.0d-10 where w < wcmn.
c arrays az, arw, afamp, afphse, apfwd, apbwd, apnet, and apohmz are 
c passed to the main program through the common block.
      do 90 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      kcapmn=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
      cf=dcmplx(y(ie+5*(imode-1)+1),y(ie+5*(imode-1)+2))
      cfp=dcmplx(y(ie+5*(imode-1)+3),y(ie+5*(imode-1)+4))
      dum17=w*xmn**2*kcapmn/8.0d7
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0d0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      az(iz+1)=z
      arw(iz+1)=radius(z)
      anrho(iz+1)=rho(z)/1.72d-8
      afamp(iz+1,imode)=
     &dsqrt(y(ie+5*(imode-1)+1)**2+y(ie+5*(imode-1)+2)**2)
      if(y(ie+5*(imode-1)+2).ne.0.0d0.or. 
     &y(ie+5*(imode-1)+1).ne.0.0d0) afphse(iz+1,imode)
     &=datan2(y(ie+5*(imode-1)+2),y(ie+5*(imode-1)+1))         
      if(y(ie+5*(imode-1)+2).eq.0.0d0.and. 
     &y(ie+5*(imode-1)+1).eq.0.0d0) afphse(iz+1,imode)=0.0d0
      cfstar=dcmplx(y(ie+5*(imode-1)+1),-y(ie+5*(imode-1)+2))
      cfpstr=dcmplx(y(ie+5*(imode-1)+3),-y(ie+5*(imode-1)+4))
      cdum30=cf*cfpstr-cfstar*cfp
      apnet(iz+1,imode)=0.5d0*ci*dum17*cdum30
      apohmz(iz+1,imode)=dum17*dum18*dum19*
     &(y(ie+5*(imode-1)+1)**2+y(ie+5*(imode-1)+2)**2)
      if((w/realc)**2.lt.(xmn/rw)**2) go to 29
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      ckz=cdsqrt(ckz2)
      kzr=ckz
      cdum10=cdexp(ci*kzr*z)
      cdum20=cdexp(-ci*kzr*z)
      cdum5=ci*kzr*cf
      cdum6=2.0d0*ci*kzr
      cffwd=(cdum5+cfp)*cdum20/cdum6
      cfbwd=(cdum5-cfp)*cdum10/cdum6
c note: the forward component of cf is cffwd*cdexp(ci*kzr*z)
c       the backward component of cf is cfbwd*cdexp(-ci*kzr*z)
      apfwd(iz+1,imode)=kzr*dum17*cdabs(cffwd)**2
      apbwd(iz+1,imode)=kzr*dum17*cdabs(cfbwd)**2
      if(iz.ne.izstep) go to 39
      ffwdr=cffwd*cdum10
      ffwdi=-ci*cffwd*cdum10
      phsout(imode)=datan2(ffwdi,ffwdr)
c     phsout(imode)=afphse(iz+1)
      delphs(imode)=phsout(imode)-phsin(imode)
      go to 39
   29 continue
      apfwd(iz+1,imode)=1.0d-10
      apbwd(iz+1,imode)=1.0d-10
      if(iz.ne.izstep) go to 39
      phsout(imode)=1.0d10
      delphs(imode)=1.0d10
   39 continue
   90 continue
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
      do 54 imode=1,imdmax
      pohmsm(imode)=y(ie+5*(imode-1)+5)
   54 continue

      if(idiag.ne.1) go to 58
c store final coordinate arrays (arc,..agamma) at right end. these arrays
c are passed to the main program through the common block for diagnosis.
      gmcsv=0.0d0
      ptcsv=0.0d0
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
c check constancy of gamma 
      npp=pp/(me*realc)
      npz=pz/(me*realc)
      gamma2=dsqrt(1.0d0+npp**2+npz**2)
      agmcsv(it,ipz,irc,ipsi,iphi)=(gamma-gamma2)/gamma0
      gmcsv=gmcsv+dabs(agmcsv(it,ipz,irc,ipsi,iphi))/dfloat(itotal)

c check constancy of ptheta (when all modes have im=0, ptcsv is evaluated,
c otherwise, ptcsv is artificially set to 1.0d50)
      imode=1
   55 if(iim(imode).ne.0) ptcsv=1.0d50
      if(iim(imode).ne.0) go to 59
      imode=imode+1
      if(imode.le.imdmax) go to 55
      z=y(ie+5*imdmax+1)
      call bfield(z,bz,bzp)
      oe=e*bz/(me*realc)
      vp=pp/(gamma*me)
      oc=oe/gamma
      rl=vp/oc
      r=dsqrt(rc**2+rl**2+2.0d0*rc*rl*dsin(phi-psi))
      ptheta=0.5d0*me*oe*(rl**2-rc**2)
      do 57 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      kmn=xmn/radius(z)
      cf=dcmplx(y(ie+5*(imode-1)+1),y(ie+5*(imode-1)+2))
      ptheta=ptheta
     &-(e/realc)*r*kmn*cf*bj(1,kmn*r)*cdexp(-ci*w*time)
   57 continue
      aptf(it,ipz,irc,ipsi,iphi)=ptheta
      aptcsv(it,ipz,irc,ipsi,iphi)=(aptf(it,ipz,irc,ipsi,iphi)
     &-apti(it,ipz,irc,ipsi,iphi))/apti(it,ipz,irc,ipsi,iphi)
      ptcsv=ptcsv+dabs(aptcsv(it,ipz,irc,ipsi,iphi))/dfloat(itotal)
   59 continue
      i=i+7
   60 continue
   58 continue

c calculate cbc(cg) for mode imdbcr
      w=aw(imdbcr)
      im=iim(imdbcr)
      in=iin(imdbcr)
      xmn=axmn(iabs(im)+1,in)
      fr=y(ie+5*(imdbcr-1)+1)
      fi=y(ie+5*(imdbcr-1)+2)
      fpr=y(ie+5*(imdbcr-1)+3)
      fpi=y(ie+5*(imdbcr-1)+4)
      cf=dcmplx(fr,fi)
      cfp=dcmplx(fpr,fpi)
      zf=y(ie+5*imdmax+1)
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
      dimension weight(21,41,1,1,23)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100)
      dimension t(1000)
      dimension axmn(9,8),aw(11),apin(11),iis(11),iim(11),iin(11)        
      dimension cx(11),xr(11),xi(11)
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode1/axmn,aw,apin,iis,iim,iin,imdmax,imdbcl,imdbcr
      common/mode2/w,pin,xmn,kcapmn,is,im
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      common/cutoff/icut
      ie=7*itotal
      z=y(ie+5*imdmax+1)
      rw=radius(z)
      call bfield(z,bz,bzp)
      oe=e*bz/(me*realc)
      oep=e*bzp/(me*realc)
      oepoe=oep/oe
      do 30 imode=1,imdmax
      cx(imode)=dcmplx(0.0d0,0.0d0)
   30 continue
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
      dy(i+1)=-0.5d0*oepoe*(rc-rl*dsin(phi-psi))
      dy(i+2)=-0.5d0*oepoe*(rl/rc)*dcos(phi-psi)
      dy(i+3)=0.5d0*me*oep*(rl+rc*dsin(phi-psi))
      dy(i+4)=me*oe/pz+0.5*me*(oep/pp)*rc*dcos(phi-psi)
      dy(i+5)=-0.5d0*me*oep*(pp/pz)*(rl+rc*dsin(phi-psi))
      dy(i+6)=1.0d0/vz
      dy(i+7)=0.0d0
      do 10 imode=1,imdmax
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      is=iis(imode)
      xmn=axmn(iabs(im)+1,in)
      kcapmn=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
      kmn=xmn/rw
      wcmn=kmn*realc
      kmnrc=kmn*rc
      kmnrl=kmn*rl
      ism=is-im
      as=w*time-dfloat(is)*phi+dfloat(is-im)*psi
     &+(dfloat(im)-dfloat(is)/2.0d0)*pi
      cdum18=cdexp(-ci*as)
      cf=dcmplx(y(ie+5*(imode-1)+1),y(ie+5*(imode-1)+2))
      cfp=dcmplx(y(ie+5*(imode-1)+3),y(ie+5*(imode-1)+4))
      cdum1=w*cf+ci*vz*cfp
      duma=kmn*e/(pz*realc)
ccc the following two statements model the converter section by tapering 
c   the electron charge in the first and last section
c     if(z.gt.zmark(izmark-1)) 
c    &duma=duma*(zmark(izmark)-z)/(zmark(izmark)-zmark(izmark-1))
c     if(z.lt.zmark(2)) duma=duma*(z-zmark(1))/(zmark(2)-zmark(1))
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
      dy(i+1)=dy(i+1)+duma*ci*(-cdum1*dum23/oc
     &+0.5d0*kmnrl*cf*(dum45-dum67))*cdum18
      dy(i+2)=dy(i+2)+dumb*(-cdum1*dfloat(is-im)*dum113/(oc*kmnrc)
     &+0.5d0*kmnrl*cf*(dum45+dum67))*cdum18
      dy(i+3)=dy(i+3)+dumc*ci*cdum1*dum112*cdum18
      dy(i+4)=dy(i+4)+
     &dumd*(-dfloat(is)*cdum1/kmnrl+kmn*vp*cf)*dum113*cdum18
      dy(i+5)=dy(i+5)+dume*cfp*dum112*cdum18
      dy(i+7)=dy(i+7)+dumf*ci*cf*dum112*cdum18
      dumg=8.0d0*kmn*rib/(xmn**2*kcapmn*realc)
      cx(imode)=cx(imode)+dumg*weight(it,ipz,irc,ipsi,iphi)*(pp/pz)
     &*dum112*cdexp(ci*as)
   10 continue
      i=i+7
   60 continue
c
      do 20 imode=1,imdmax
ccc the following two statements model the converter section by tapering 
c   the electron charge in the first and last section
c     if(z.gt.zmark(izmark-1)) cx(imode)=
c    &cx(imode)*(zmark(izmark)-z)/(zmark(izmark)-zmark(izmark-1))
c     if(z.lt.zmark(2)) 
c    &cx(imode)=cx(imode)*(z-zmark(1))/(zmark(2)-zmark(1))
      w=aw(imode)
      im=iim(imode)
      in=iin(imode)
      is=iis(imode)
      xmn=axmn(iabs(im)+1,in)
c
      if(icut.ne.1) go to 123
      if((w/realc)**2.gt.(xmn/rw)**2) go to 123
      dy(ie+5*(imode-1)+1)=0.0d0
      dy(ie+5*(imode-1)+2)=0.0d0
      dy(ie+5*(imode-1)+3)=0.0d0
      dy(ie+5*(imode-1)+4)=0.0d0
      dy(ie+5*(imode-1)+5)=0.0d0
      go to 20
  123 continue
c
      kcapmn=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
      kmn=xmn/rw
      wcmn=kmn*realc
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      kz2r=ckz2
      kz2i=-ci*ckz2
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      dum17=w*xmn**2*kcapmn/8.0d7
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0d0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      xr(imode)=cx(imode)
      xi(imode)=-ci*cx(imode)
      dy(ie+5*(imode-1)+1)=y(ie+5*(imode-1)+3)
      dy(ie+5*(imode-1)+2)=y(ie+5*(imode-1)+4)
      dy(ie+5*(imode-1)+3)=
     &-kz2r*y(ie+5*(imode-1)+1)+kz2i*y(ie+5*(imode-1)+2)-xr(imode)
      dy(ie+5*(imode-1)+4)=
     &-kz2i*y(ie+5*(imode-1)+1)-kz2r*y(ie+5*(imode-1)+2)-xi(imode)
      dy(ie+5*(imode-1)+5)=dum17*dum18*dum19*
     &(y(ie+5*(imode-1)+1)**2+y(ie+5*(imode-1)+2)**2)
   20 continue
      dy(ie+5*imdmax+1)=1.0d0
      return
      end
c    ***********************************************************************
      subroutine cldtst(imode,fhzmin,fhzmax,iwmax)
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
      dimension axmn(9,8),aw(11),apin(11),iis(11),iim(11),iin(11)        
      common/const/ci,pi,realc,me,e
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/mode1/axmn,aw,apin,iis,iim,iin,imdmax,imdbcl,imdbcr
      common/mode2/w,pin,xmn,kcapmn,is,im
      common/local1/pbwdi,pfwdf
      external cbcckt
c 
      im=iim(imode)
      in=iin(imode)
      xmn=axmn(iabs(im)+1,in)
      kcapmn=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
c
      if(iwmax.eq.1) return
      z1=zmark(1)
      rw=radius(z1)
      wmin=fhzmin*2.0d0*pi
      if((wmin/realc)**2.le.xmn**2/rw**2) write(6,1) im,in
      if((wmin/realc)**2.le.xmn**2/rw**2) return
    1 format(' *fhzmin below cut-off of te(',i2,',',i2,
     &') mode in the first section, no cold test*')
c
c    specify no of steps for z-integration.
      izstep=2001
c    specify incident (forward) wave power in w at left end
      pin=1.0d0
      tlmin=0.0d0
      rlmin=0.0d0
      do 100 iw=1,iwmax
      fstep=(fhzmax-fhzmin)/dfloat(iwmax-1)
      fhz=fhzmin+fstep*(iw-1)
      afreq(iw)=fhz
      w=2.0d0*pi*fhz
c    angular frequency of the incident wave (in radian/sec, denoted by w)
c    is passed to subprograms through the common block.
c
c    guess complex reflection coefficient
      cguess=dcmplx(0.5d0,0.5d0)
      if(iw.eq.1) cg(1)=cguess

c    calculate reflection coefficient (referred to the z=0 plane)
      ep1=1.0d-5
      ep2=ep1
      iroot=1
      imaxit=50
      call muller(0,iroot,cg,imaxit,ep1,ep2,cbcckt,.false.)
      value=cdabs(cbcckt(cg(1)))
      tldb=10.0d0*dlog10(pfwdf/pin)
      rldb=10.0d0*dlog10(pbwdi/pin)
      atldb(iw)=tldb
      arldb(iw)=rldb
      if(atldb(iw).lt.tlmin) tlmin=atldb(iw)
      if(arldb(iw).lt.rlmin) rlmin=arldb(iw)
  100 continue
      write(6,2) im,in
    2 format(/30x,' *** te(',i2,i2,') mode ***') 
      call sscale(afreq(1),afreq(iwmax),tlmin,0.0d0)
      call splot(afreq,atldb,iwmax,
     &'tt  transmission loss(db) vs freq.(hz)$')
      write(6,2) im,in
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
      common/mode2/w,pin,xmn,kcapmn,is,im
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
      cf0=dcmplx(1.0d0,0.0d0)
      do 20 i=1,2
      cdum3=cf0*(cdum1+cg*cdum2)
      cdum4=cf0*(cdum1-cg*cdum2)
      y(1)=cdum3
      y(2)=-ci*cdum3
      y(3)=ci*ckz*cdum4
      y(4)=ckz*cdum4
      y(5)=0.0d0
      y(6)=z1
c
      cf=dcmplx(y(1),y(2))
      cfp=dcmplx(y(3),y(4))
      kzr=ckz
      cdum10=cdexp(ci*kzr*z1)
      cdum20=cdexp(-ci*kzr*z1)
      cdum5=ci*kzr*cf
      cdum6=2.0d0*ci*kzr
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
      cdum6=2.0d0*ci*kzr
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
      common/local1/pbwdi,pfwdf
      common/mode2/w,pin,xmn,kcapmn,is,im
      z=y(6)
      rw=radius(z)
      ckz2=(w/realc)**2-xmn**2*closs(z)/rw**2
      kz2r=ckz2
      kz2i=-ci*ckz2
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      wcmn=xmn*realc/rw
      dum17=w*xmn**2*kcapmn/8.0d7
      dum18=(delta/rw)*(wcmn/realc)**2
      dum19=1.0d0+(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w/wcmn)**2
      dy(1)=y(3)
      dy(2)=y(4)
      dy(3)=-kz2r*y(1)+kz2i*y(2)
      dy(4)=-kz2i*y(1)-kz2r*y(2)
      dy(5)=dum17*dum18*dum19*(y(1)**2+y(2)**2)
      dy(6)=1.0d0
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
      common/mode2/w,pin,xmn,kcapmn,is,im
      delta=dsqrt(realc**2*rho(z)/(2.0d0*pi*w*9.0d9))
      rw=radius(z)
      wcmn=xmn*realc/rw
      cdum1=(1.0d0+ci)*delta/rw
      cdum2=(dfloat(im**2)/(xmn**2-dfloat(im**2)))*(w**2/wcmn**2)
      closs=1.0d0-cdum1*(1.0d0+cdum2)
      return
      end
c **********************************************************************
      subroutine bfield(z,bz,bzp)
c this subroutine calculates the external b-field (on the axis) and its
c derivative at position z.
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      common/bext/bz0
      bz=bz0
      bzp=0.0d0
      return
      end
c ***********************************************************************
      subroutine coord(arc,apsi,app,aphi,apz,atime,agamma,weight)
      implicit  real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension arc(21,41,1,1,23),apsi(21,41,1,1,23),
     &app(21,41,1,1,23),aphi(21,41,1,1,23),apz(21,41,1,1,23),
     &atime(21,41,1,1,23),agamma(21,41,1,1,23),
     &weight(21,41,1,1,23)
      dimension axmn(9,8),aw(11),apin(11),iis(11),iim(11),iin(11)        
      common/const/ci,pi,realc,me,e
      common/mode1/axmn,aw,apin,iis,iim,iin,imdmax,imdbcl,imdbcr
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
c
c this subroutine sets initial electron coordinates (in real and momentum  
c spaces) and weighing factors
c electron injection times (atime) are evenly spaced over the wave period
c of the mode with the lowest frequency.
      w=aw(1)
      do 10 imode=1,imdmax
      if(w.lt.aw(imode)) w=aw(imode)
   10 continue
      suma=0.0d0
      tstep=2.0d0*pi/(w*dfloat(itmax))
      phistp=2.0d0*pi/dfloat(iphimx)
      psistp=2.0d0*pi/dfloat(ipsimx)
      if(iim(1).ne.0.and.imdmax.eq.1)
     &psistp=2.0d0*pi/dfloat(ipsimx*iim(1))
      pzstep=(pzmax-pzmin)/dfloat(ipzmax)
      rcstep=(rc2-rc1)/dfloat(ircmax)
      do 120 it=1,itmax
      do 120 ipz=1,ipzmax
      pz=pzmin+pzstep*(dfloat(ipz)-0.5d0)
      pp=dsqrt(p0-pz)*dsqrt(p0+pz)
      do 120 irc=1,ircmax
      rc=rc1+rcstep*(dfloat(irc)-0.5d0)
      if(rc.eq.0.0d0) rc=1.0d-7
      if(delpz.eq.0.0d0) fac2=rc
      if(delpz.ne.0.0d0)
     &fac2=dexp(-0.5d0*((pz-pz0)/delpz)**2)*rc
      do 120 ipsi=1,ipsimx
      do 120 iphi=1,iphimx
      arc(it,ipz,irc,ipsi,iphi)=rc
      if(ipsimx.eq.1) apsi(it,ipz,irc,ipsi,iphi)=0.0d0
      if(ipsimx.ne.1)
     &apsi(it,ipz,irc,ipsi,iphi)=psistp*(dfloat(ipsi)-0.5d0)
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
           real*8  aimag
      external fn
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
      if ( kount .gt. maxit ) write(6,1492)
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
      write(6,1) title
    1 format(1h1,///,15x,101a1)
      do 17 ny=1,51
      if(mod(ny-1 , 5) .ne. 0)go to 65
      ay = ymax-(ny-1)*dy
      if(dabs(ay) .le. .001*dabs((ymax-ymin)))ay = 0.0
      write(6, 9) ay, (arr(ny,n),n = 1,101)
      go to 17
   65 write(6,10)     (arr(ny,n),n = 1,101)
   17 continue
    9 format(1pd14.4,1x,101a1)
   10 format(15x,101a1)
      write(6,12)  (xx(n),n=1,101,20)
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
      write(6,1) (title(it),it=1,nxap)
    1 format(1h0,/,12x,101a1)
      do 17 ny=1,nyap
      if(mod(ny-1 , nyp) .ne. 0)go to 65
      ay = ymax-(ny-1)*dy
      if(dabs(ay) .le. .001*dabs((ymax-ymin)))ay = 0.0
      write(6, 9) ay, (arr(ny,n),n = 1,nxap)
      go to 17
   65 write(6,10)     (arr(ny,n),n = 1,nxap)
   17 continue
    9 format(1x,1pd11.4,101a1)
   10 format(12x,101a1)
      write(6,12)  (xx(n),n=1,nxap,nxp2)
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
c *******************************************************************
      subroutine bdiag
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension arc(21,41,1,1,23),apsi(21,41,1,1,23),
     &app(21,41,1,1,23),aphi(21,41,1,1,23),apz(21,41,1,1,23),
     &atime(21,41,1,1,23),agamma(21,41,1,1,23),
     &weight(21,41,1,1,23)
      dimension zmark(100),rwl(100),rwr(100),rhol(100),rhor(100) 
      dimension appy(1001),apzy(1001),akey(1001),
     &appx(1001),apzx(1001),akex(1001)
      dimension atmy(1001),atmx(1001)
      dimension axmn(9,8),aw(11),apin(11),iis(11),iim(11),iin(11)        
      common/ckt/zmark,rwl,rwr,rhol,rhor,izmark,izstep
      common/const/ci,pi,realc,me,e
      common/ebeam1/gamma0,p0,pz0,pzmin,pzmax,rc1,rc2,rib,
     &delpz,itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
      common/ebeam2/weight
      common/diagb/arc,apsi,app,aphi,apz,atime,agamma
      common/mode1/axmn,aw,apin,iis,iim,iin,imdmax,imdbcl,imdbcr
c
      pp0=dsqrt(p0**2-pz0**2)
      ppmin=1.0d10
      ppmax=-1.0d10
      pzmin=1.0d10
      pzmax=-1.0d10
      gmmin=1.0d10
      gmmax=-1.0d10
      sumpp=0.0d0
      sumpz=0.0d0
      sumke=0.0d0
      sumpp2=0.0d0
      sumpz2=0.0d0
      sumke2=0.0d0
      sumwt=0.0d0
      do 100 it=1,itmax
      do 100 ipz=1,ipzmax
      do 100 irc=1,ircmax
      do 100 ipsi=1,ipsimx
      do 100 iphi=1,iphimx
      pp=app(it,ipz,irc,ipsi,iphi)
      pz=apz(it,ipz,irc,ipsi,iphi)
      gm=agamma(it,ipz,irc,ipsi,iphi)
      wt=weight(it,ipz,irc,ipsi,iphi)
      sumpp2=sumpp2+(pp-pp0)**2*wt
      sumpz2=sumpz2+(pz-pz0)**2*wt
      sumke2=sumke2+(gm-gamma0)**2*wt
      sumwt=sumwt+wt
      sumpp=sumpp+pp*wt
      sumpz=sumpz+pz*wt
      sumke=sumke+(gm-1.0d0)*wt
      if(pp.gt.ppmax) ppmax=pp
      if(pp.lt.ppmin) ppmin=pp
      if(pz.gt.pzmax) pzmax=pz
      if(pz.lt.pzmin) pzmin=pz
      if(gm.gt.gmmax) gmmax=gm
      if(gm.lt.gmmin) gmmin=gm
  100 continue
      if(ppmax.ne.ppmin) go to 101
      ppmin=(1.0d0-0.000001d0)*ppmin
      ppmax=(1.0d0+0.000001d0)*ppmax
  101 continue
      if(pzmax.ne.pzmin) go to 102
      pzmin=(1.0d0-0.000001d0)*pzmin
      pzmax=(1.0d0+0.000001d0)*pzmax
  102 continue
      if(gmmax.ne.gmmin) go to 103
      gmmin=(1.0d0-0.000001d0)*gmmin
      gmmax=(1.0d0+0.000001d0)*gmmax
  103 continue
      sdvpp0=dsqrt(sumpp2/sumwt)/pp0
      sdvpz0=dsqrt(sumpz2/sumwt)/pz0
      sdvke0=dsqrt(sumke2/sumwt)/(gamma0-1.0d0)
      avgpp=sumpp/sumwt
      avgpz=sumpz/sumwt
      avgke=sumke/sumwt
      rtopp=avgpp/pp0
      rtopz=avgpz/pz0
      rtoke=avgke/(gamma0-1.0d0)
      sumpp2=0.0d0
      sumpz2=0.0d0
      sumke2=0.0d0
      do 110 it=1,itmax
      do 110 ipz=1,ipzmax
      do 110 irc=1,ircmax
      do 110 ipsi=1,ipsimx
      do 110 iphi=1,iphimx
      pp=app(it,ipz,irc,ipsi,iphi)
      pz=apz(it,ipz,irc,ipsi,iphi)
      gm=agamma(it,ipz,irc,ipsi,iphi)
      wt=weight(it,ipz,irc,ipsi,iphi)
      sumpp2=sumpp2+(pp-avgpp)**2*wt
      sumpz2=sumpz2+(pz-avgpz)**2*wt
      sumke2=sumke2+(gm-(avgke+1.0d0))**2*wt
  110 continue
      sdvpp=dsqrt(sumpp2/sumwt)/avgpp
      sdvpz=dsqrt(sumpz2/sumwt)/avgpz
      sdvke=dsqrt(sumke2/sumwt)/avgke
      write(6,5) sumwt,sdvpp0,sdvpp,rtopp,sdvpz0,sdvpz,rtopz,
     &sdvke0,sdvke,rtoke
    5 format(/' sumwt=',1pd11.4/
     &' (sdv of pp about pp0)/pp0=',1pd11.4/
     &' (sdv of pp about avgpp)/avgpp=',1pd11.4/
     &' avgpp/pp0=',1pd11.4/
     &' (sdv of pz about pz0)/pz0=',1pd11.4/
     &' (sdv of pz about avgpz)/avgpz=',1pd11.4/
     &' avgpz/pz0=',1pd11.4/
     &' (sdv of ke about ke0)/ke0=',1pd11.4/
     &' (sdv of ke about avgke)/avgke=',1pd11.4/
     &' avgke/ke0=',1pd11.4/)
      do 1000 ino=1,10
      no=dfloat(itotal)/(3.0d0**(ino-1))
      ipoint=idint(no)
      if(ipoint.lt.10) go to 1000
      if(ipoint.gt.1000) go to 1000
      write(6,1) ipoint,ino
    1 format(/' *** no of sampling points=',i8,', ino=',i5,' ***')
      ppstep=(ppmax-ppmin)/dfloat(ipoint)
      pzstep=(pzmax-pzmin)/dfloat(ipoint)
      gmstep=(gmmax-gmmin)/dfloat(ipoint)
      do 200 i=1,ipoint
      appy(i)=0.0d0
      apzy(i)=0.0d0
      akey(i)=0.0d0
  200 continue
      ppymin=1.0d10
      ppymax=-1.0d10
      pzymin=1.0d10
      pzymax=-1.0d10
      keymin=1.0d10
      keymax=-1.0d10
      ppxmin=ppmin/pp0
      ppxmax=ppmax/pp0
      pzxmin=pzmin/pz0
      pzxmax=pzmax/pz0
      kexmin=(gmmin-1.0d0)/(gamma0-1.0d0)
      kexmax=(gmmax-1.0d0)/(gamma0-1.0d0)
      do 300 i=1,ipoint
      pp=ppmin+ppstep*dfloat(i-1)
      pz=pzmin+pzstep*dfloat(i-1)
      gm=gmmin+gmstep*dfloat(i-1)
      appx(i)=pp/pp0
      apzx(i)=pz/pz0
      akex(i)=(gm-1.0d0)/(gamma0-1.0d0)
      do 400 it=1,itmax
      do 400 ipz=1,ipzmax
      do 400 irc=1,ircmax
      do 400 ipsi=1,ipsimx
      do 400 iphi=1,iphimx
      if(app(it,ipz,irc,ipsi,iphi).ge.pp.and.
     &app(it,ipz,irc,ipsi,iphi).lt.(pp+ppstep))
     &appy(i)=appy(i)+weight(it,ipz,irc,ipsi,iphi)
      if(apz(it,ipz,irc,ipsi,iphi).ge.pz.and.
     &apz(it,ipz,irc,ipsi,iphi).lt.(pz+pzstep))
     &apzy(i)=apzy(i)+weight(it,ipz,irc,ipsi,iphi)
      if(agamma(it,ipz,irc,ipsi,iphi).ge.gm.and.
     &agamma(it,ipz,irc,ipsi,iphi).lt.(gm+gmstep))
     &akey(i)=akey(i)+weight(it,ipz,irc,ipsi,iphi)
  400 continue
      if(appy(i).lt.ppymin) ppymin=appy(i)
      if(appy(i).gt.ppymax) ppymax=appy(i)
      if(apzy(i).lt.pzymin) pzymin=apzy(i)
      if(apzy(i).gt.pzymax) pzymax=apzy(i)
      if(akey(i).lt.keymin) keymin=akey(i)
      if(akey(i).gt.keymax) keymax=akey(i)
  300 continue
      call sscale(ppxmin,ppxmax,0.0d0,ppymax)
      call splot(appx,appy,ipoint,'pp  ppy vs pp$')
      call sscale(pzxmin,pzxmax,0.0d0,pzymax)
      call splot(apzx,apzy,ipoint,'zz  pzy vs pz$')
      call sscale(kexmin,kexmax,0.0d0,keymax)
      call splot(akex,akey,ipoint,'kk  key vs ke$')
 1000 continue
c
c******************************
c the following statement bypasses the (unsuccessful) current vs time
c diagnostics
      if(gamma0.gt.0.0d0) go to 2001
c******************************
      tmmin=1.0d10
      tmmax=-1.0d10
      do 600 it=1,itmax
      do 600 ipz=1,ipzmax
      do 600 irc=1,ircmax
      do 600 ipsi=1,ipsimx
      do 600 iphi=1,iphimx
      tm=atime(it,ipz,irc,ipsi,iphi)
      if(tm.gt.tmmax) tmmax=tm
      if(tm.lt.tmmin) tmmin=tm
      y=dsin(aphi(it,ipz,irc,ipsi,iphi))
      x=dcos(aphi(it,ipz,irc,ipsi,iphi))
      phi1=aphi(it,ipz,irc,ipsi,iphi)/(2.0d0*pi)
      aphi(it,ipz,irc,ipsi,iphi)=datan2(y,x)
      phi2=aphi(it,ipz,irc,ipsi,iphi)/(2.0d0*pi)
c     write(6,7) it,ipz,irc,ipsi,iphi,phi1,phi2
    7 format(5(i3,','),1pd11.4,',',1pd11.4)
  600 continue
      w=aw(1)
      do 700 imode=1,imdmax
      if(w.lt.aw(imode)) w=aw(imode)
  700 continue
      period=2.0d0*pi/w
      tratio=(tmmax-tmmin)/period
      tmmin=tmmin-1.0d-7*period
      do 2000 ino1=1,5
      iphino=2**(ino1+1)
      phistp=2.0d0*pi/dfloat(iphino)
      do 2000 ino2=1,7
      itno=2**(ino2+1)
      if(iphino.gt.iphimx) go to 2000
      if(itno.gt.itmax) go to 2000
      if(itno.gt.1000) go to 2000
      tmstep=period/dfloat(itno)
      write(6,2) tratio, iphino,itno
    2 format(/' *** (tmax-tmin)/period=',1pd11.4,
     &', iphino=',i4,', itno=',i4,' ***')
      do 800 it=1,itno
      atmy(it)=0.0d0
  800 continue
      tmymin=1.0d10
      tmymax=-1.0d10
      tmxmin=0.0d0
      tmxmax=1.0d0
      do 810 i=1,itno
      tm=tmmin+tmstep*dfloat(i-1)
      atmx(i)=(tm-tmmin)/period
      inmax=1+idint(tratio)
      do 820 it=1,itmax
      do 820 ipz=1,ipzmax
      do 820 irc=1,ircmax
      do 820 ipsi=1,1
      do 820 iphi=1,iphimx
      gm=agamma(it,ipz,irc,ipsi,iphi)
      z=zmark(izmark)
      call bfield(z,bz,bzp)
      oe=e*bz/(me*realc)
      vp=pp/(gm*me)
      oc=oe/gm
      rl=vp/oc
      do 820 in=1,inmax
      if(aphi(it,ipz,irc,ipsi,iphi).lt.0.0d0.or.
     &aphi(it,ipz,irc,ipsi,iphi).gt.phistp) go to 820
      if(atime(it,ipz,irc,ipsi,iphi).ge.
     &(tm+dfloat(in-1)*period).and.atime(it,ipz,irc,ipsi,iphi)
     &.lt.(tm+dfloat(in-1)*period+tmstep))
     &atmy(i)=atmy(i)+weight(it,ipz,irc,ipsi,iphi)*rl
  820 continue
      if(atmy(i).gt.tmymax) tmymax=atmy(i)
      if(atmy(i).lt.tmymin) tmymin=atmy(i)
  810 continue
      if(tmymin.ge.tmymax) tmymin=0.0d0
      call sscale(tmxmin,tmxmax,tmymin,tmymax)
      call splot(atmx,atmy,itno,'**  i vs t$')
 2000 continue
 2001 continue
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

