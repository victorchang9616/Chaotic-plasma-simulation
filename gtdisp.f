      program gtdisp
c
c   (all in double precion)
c this program calculates complex kz as a function of real w for the
c gyro-twt and carm amplifier by using the dispersion relation eq.(47)
c or eq(43) in 'gain and bandwidth of the gyro-twt and carm amplifiers'
c by k. r. chu and a. t. lin (ieee trans. on plasma science. vol. 16,
c pp.90-104 (1988)). 
c
c parameters in this program are normalized to the waveguide wall radius
c according to the scheme on p.2358 of 'k. r. chu, phys. fluids, vol. 12,
c pp.2354-2364 (1978)'. if the wall radius is given, some of the key
c physical parameters (e.g. frequency in hz, b-field in gauss, gain in
c db/cm) are also printed in parentheses.
c
c note: comment statements in the main program beginning with ccc 
c       call for action by the user.
c 
c date of first version: 1987
c date of this version: may 15, 1999
c
      implicit  real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension ckz(20),axmn(9,8)
      dimension  akzr(201,10),akzi(201,10),aw(201),akzr0(201,10), 
     &akzi0(201,10) 
      dimension  akzi1(201,10)
      dimension afhz(201), again1(201,10),again2(201,10)
      common ci,pi,e,e2,realc,realc2,me,
     &ppmax,ppmin,pzmax,pzmin,divpp,divpz,
     &nu,r1,r2,rc0,omegae, w,eps,a,pp0,pz0,p0,gamma0,delpp,delpz,
     &ppc,xmn,kmn,delta,hsm,tsm,usm,
     &is,im,in,ippmax,ipzmax,icont1
      external cmask0,cmask1
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
c
c universal constants
c
      ci=dcmplx(0.0d0,1.0d0)
      pi=3.1415926d0
      e=4.8032d-10
      e2=e*e
      realc=2.99792d10
      realc2=realc*realc
      me=9.1095d-28
c
c
ccc1  specify plot instruction (1:yes, 0:no)
c
      iplot=0
c
ccc2  specify print instruction
c
      iskip=1
c
c explanation: the program prints out kz vs w by skipping 'iskip' data
c              points for every data point printed. the skipped data
c              will be stored for the plots.
c
ccc3  specify operating mode
c
      ismain=1
      immain=1
      inmain=1
c
ccc4  specify cyclotron harmonic number 'is' of mode(s) to be calculated
c
      do 600 is=1,1
c
ccc5  specify azimuthal mode number 'im' of mode(s) to be calculated
c
      do 600 imm=2,2
      im=imm-1
c
ccc6  specify radial mode number 'in' of mode(s) to be calculated
c
      do 600 in=1,1
c
      xmain=axmn(iabs(immain)+1,inmain)
      xmn=axmn(iabs(im)+1,in)
c
ccc7  specify normalized electron guiding center radius 'rc0'
c     (normalized to wall radius)
c
      rc0=0.35
c     ipeak=1
c     rc0=axmn(iabs(ismain-immain)+1,ipeak)/xmain
c     if(ismain.eq.immain.and.ipeak.eq.1) rc0=0.0e0
c     if(ismain.eq.immain.and.ipeak.ne.1) rc0=axmn(1,ipeak-1)/xmain
c
ccc8  specify wall radius 'rwcm' in cm (set rwcm to any value if
c     only the normalized quantities are needed)
c
      rwcm=0.2654
c
      fchz=xmain*realc/(2.0*pi*rwcm)
      lamdac=realc/fchz
      print 10, rwcm,fchz,lamdac
   10 format(/' rw=',1pe12.5,' cm'/
     &' cutoff freq of operating mode=',1pe12.5,' hz'/
     &' cutoff wavelength of operating mode=',1pe12.5,' cm'/)
c
      print 15
   15 format(' --------------------'/)
      print 1, is,im,in,ismain,immain,inmain,xmain,xmn
    1 format('     is=',i3,',     im=',i3,',     in=',i3/
     &' ismain=',i3,', immain=',i3,', inmain=',i3,/
     &' xmain=',1pe11.4,', xmn=',1pe11.4/)
      print 15
c
ccc9  specify 'alpha' (=pp0/pz0)
c
      do 900 ialpha=4,4
      alpha=0.7+(ialpha-1)*0.1
c
      theta0=datan2(alpha,1.0d0)
        angle=theta0/pi*180.0e0
c
ccc10 specify beam voltage in kv (beam is assumed to be monoenergetic)
c
      do 900 ivb=1,1
      vbkv=102.96
c
      gamma0=1.0e0+vbkv/511.0e0
        p0=dsqrt(gamma0**2-1.0e0)
      pp0=p0*dsin(theta0)
        v0=p0/gamma0
        vp0=v0*dsin(theta0)
        vz0=v0*dcos(theta0)
        gamz0=1.0e0/dsqrt(1.0e0-vz0**2)
c
ccc11 specify beam current 'ribamp' in amp
c
      do 900 ib=1,1
      ribamp=1.00
c
ccc12 specify 'b0bg' (=b0/bg) for the operating mode
c
      do 900 ib0bg=1,1
      b0bg=1.00+0.02*(ib0bg-1)
c
      omegae=gamma0*xmain/(gamz0*ismain)
      omegae=b0bg*omegae
      bgauss=(me*realc2/e)*omegae/rwcm
      rimain=bgauss/401.936 
c 
c rimain is the main coil current of the nthu superconducting 
c magnet corresponding to a magnetic field of bgauss. 
c 
        omegac=omegae/gamma0
      rl0=pp0/omegae
      xmnrl0=xmn*rl0
      xmnrc0=xmn*rc0
      kmn=bj(im,xmn)**2*(1.0-dfloat(im)**2/xmn**2)
      call couple(is,im,xmnrc0,xmnrl0,hsm,tsm,usm)
c
c calculate the critical current by using the analytical formula of
c lau, chu, barnett, and granatstein in int. j. ir and mm waves,
c vol. 2, p.373 (1981).
c note: the analytical formula is approximately valid for delta=0
c       and b0/bg=1 only.
      vz02=vz0*vz0
      bb0=dsqrt(1.0-vz02)
      bb=b0bg*bb0
      dum1=1.0+8.0*vz02
      dum2=dsqrt(8.0*vz02*(1.0-bb**2)+64.0*vz02**2)
      ws=(bb+dum2)/dum1
      ks=(ws-bb)/(4.0*vz0)
      epsc=27.0*vz02*ks**4
      nuc=epsc*gamma0*kmn*xmn**2/(4.0*vp0**2*hsm)
      nc=me*realc2*nuc/e2
      ricamp=nc*e*vz0*realc/3.0e9
      print 13, b0bg,ricamp
   13 format(/'   ** b0/bg=',1pe11.4,' **'/
     &' ic (analytically calculated critical current)=',1pe10.3,
     &' amp'/4x,
     &'(note: ic is approx. valid for delta=0 and b0/bg=1 only)')
c
ccc13 specify normalized skin depth 'delta' (normalized to wall radius)
c
      do 900 idelta=1,1
      delta=0.0
c
      nu=3.0e9*ribamp*e/(vz0*me*realc**3)
c
c default values
      do 20 i1=1,201
      aw(i1)=0.0
      afhz(i1)=0.0
      do 20 i2=1,10 
      akzr(i1,i2)=0.0
      akzr0(i1,i2)=0.0
      akzi(i1,i2)=0.0
      akzi0(i1,i2)=0.0
      akzi1(i1,i2)=0.0
      again1(i1,i2)=0.0
      again2(i1,i2)=0.0
   20 continue
      g1max=-1.0e10
      g2max=-1.0e10
c
ccc14 specify beam axial velocity spread 'tz' (=(delta pz)/pz0)
c
      do 200 itz=1,2
      if(itz.eq.1) tz=0.0
      if(itz.eq.2) tz=0.05
      if(itz.eq.3) tz=0.1
      if(itz.eq.4) tz=0.15
      if(itz.eq.5) tz=0.2
c
c note: (1) always start from tz=0.0 and increase tz in small steps
c           (e.g. at least 3 steps from tz=0.0 to the maximum
c           meaningful value)
c       (2) the program searches all 4 roots if tz=0.0 , but searches
c           only one unstable root if tz.ne.0.0
c
      print 2,vbkv,gamma0,angle,alpha,tz,rc0,omegae,bgauss,rimain,
     &omegac,delta,b0bg 
    2 format(' vb=',1pe11.4,' kv, gamma0=',
     &1pe11.4,', angle=',1pe11.4,' degrees',/,' alpha=',1pe11.4,
     &' tz=delpz/pz0=',1pe10.3/
     &' rc0=',1pe11.4/
     &' omegae=',1pe11.4,'(',1pe11.4,
c    &' gauss, corresponding to main coil current of ',1pe11.4, 
c    &' amp in the nthu superconducting magnet)'/ 
     &' gauss, imain of nthuscm1 = ',1pe11.4,' amp)'/ 
     &' omegac=',1pe11.4/' delta=',1pe11.4/ 
     &' b0/bg=',1pe11.4/) 
c
c numerical parameters 'pzcut' and 'ipzmax' are used as follows:  
c pz-intergration is done between pz0-pzcut*(delta pz) and
c pz0+pzcut*(delta pz), subject to the condition 0 < pz < ptotal.
c the interval is divided into 'ipzmax' divisions for numerical
c intergration.
c
c ipzmax is now set at a reasonable value. it can be increased to 
c obtain more accuracy in case of small gain. it can be decreased 
c to save computing time
c 
c in the case of a cold beam, ipzmax is automatically set to 0.
c
      pzcut=2.5d0
      ipzmax=idint(1000.0*tz)+200 
c
      ippmax=1
      if(tz.eq.0.0) ipzmax=1
      if(tz.eq.0.0) iroot=4
      if(tz.ne.0.0) iroot=1
      print 3, pzcut, ipzmax
    3 format(' numerical parameters:  pzcut=',1pe9.2,
     &', ipzmax=',i4/)
c
c calculations
c
      p0=dsqrt(gamma0**2-1.0e0)
      pp0=p0*dsin(theta0)
      pz0=p0*dcos(theta0)
      delpz=tz*pz0
      v0=p0/gamma0
      vp0=v0*dsin(theta0)
      vz0=v0*dcos(theta0)
      n=me*realc2*nu/e2
      rl0=pp0/omegae
      xmnrl0=xmn*rl0
      xmnrc0=xmn*rc0
      pzmin=dmax1(0.0d0,pz0-pzcut*delpz)
      pzmax=dmin1(p0,pz0+pzcut*delpz)
      ppmax=dsqrt(p0**2-pzmin**2)
      ppmin=dsqrt(p0**2-pzmax**2)
      vpmin=ppmin/gamma0
      vpmax=ppmax/gamma0
      vzmin=pzmin/gamma0
      vzmax=pzmax/gamma0
      print 4, p0,pp0,pz0,delpz,
     &ppmin,ppmax,pzmin,pzmax,
     &v0,vp0,vz0,vpmin,vpmax,vzmin,vzmax,
     &n,nu,ribamp,rc0,xmnrc0,rl0,xmnrl0
    4 format(' p0=',1pe10.3,', pp0=',1pe10.3,', pz0=',1pe10.3/
     &' delpz=',1pe10.3/
     &' ppmin=',1pe9.2,', ppmax=',1pe9.2,/
     &' pzmin=',1pe10.3,', pzmax=',1pe10.3//
     &' v0=',1pe10.3,', vp0=',1pe10.3,', vz0=',1pe10.3/
     &' vpmin=',1pe10.3,', vpmax=',1pe10.3/
     &' vzmin=',1pe10.3,', vzmax=',1pe10.3//
     &' n=',1pe10.3,'/cm (nu=', 1pe10.3,'), ib=',1pe10.3,' amp'//
     &' rc0=',1pe11.4,', (xmn*rc0=',1pe10.3,')'/
     &' rl0=',1pe11.4,', (xmn*rl0=',1pe10.3,')'/)
c
c
      if(tz.ne.0.0) go to 411
      r1=dabs(rc0-rl0)
      r2=rc0+rl0
      if(r2.gt.1.0e0) print 21
      if(r2.gt.1.0e0) go to 600
   21 format(////' ***** beam hits wall ( rc0+rl>1 ) *****'////)
c
      kmn=bj(im,xmn)**2*(1.0e0-dfloat(im)**2/xmn**2)
      call couple(is,im,xmnrc0,xmnrl0,hsm,tsm,usm)
c
      print 19, vbkv,gamma0,omegae,omegac,b0bg,rc0,r1,r2,
     &xmnrc0,xmnrl0,hsm,tsm,usm,kmn
   19 format(' < cold beam >',' vb=',1pe11.4,'kv, gamma0=',1pe11.4,
     &/,15x,'omegae=',1pe10.3,', omegac=',1pe10.3,', b0bg=',1pe10.3/
     &15x,'rc0=',1pe10.3,', r1=',1pe10.3,', r2=',1pe10.3/
     &15x,'xmnrc0=',1pe10.3,', xmnrl0=',1pe10.3,
     &/,14x,' hsm=',1pe10.3,', tsm=',1pe10.3,
     &',usm=',1pe10.3,', kmn=',1pe10.3/)
  411 continue
c
c calculate normalization constant for the distribution function
c
      suma=0.0e0
      divpz=(pzmax-pzmin)/(2.0e0*ipzmax)
      do 120 iz=1,ipzmax
      pz=pzmin+divpz*(2*iz-1)
      fac2=dexp(-(pz-pz0)**2/(2.0*delpz**2+1.0e-20))
      do 130 ip=1,ippmax
      pp=dsqrt(p0**2-pz**2)
      suma=suma+   fac2
  130 continue
  120 continue
      a=1.0e0/(2.0e0*pi*suma)
      sumg=0.0e0
      do 140 iz=1,ipzmax
      pz=pzmin+divpz*(2*iz-1)
      fac2=dexp(-(pz-pz0)**2/(2.0*delpz**2+1.0e-20))
      do 150 ip=1,ippmax
      pp=dsqrt(p0**2-pz**2)
      sumg=sumg+2.0e0*pi *a*fac2
  150 continue
  140 continue
      print 5, sumg
    5 format(/' distribution function normalized to ',1pe11.4/)
c
      kzrmax=-1.0e10
      kzimax=-1.0e10
c
c calculate points of intersections between waveguide mode and
c beam-wave resonsnce line
c     wc = cut-off w (w denotes angular frequency)
c     wg = w of grazing intersection
c     kzg = kz of grazing intersection
c     w1 = w of lower intersecting point
c     w2 = w of upper intersecting point
c     kz1 = kz of lower intersecting point
c     kz2 = kz of upper intersecting point
c
      wc=xmn
      if(is.ne.ismain) go to 510
      if(im.ne.immain) go to 510
      if(in.ne.inmain) go to 510
      if(b0bg.gt.1.0e0) go to 510
      kzg=gamz0*vz0*xmn
      wg=gamz0*xmn
      w1=wg
      w2=wg
      kz1=kzg
      kz2=kzg
      go to 520
  510 f11=(xmn/(is*omegac*gamz0))**2
      if(f11.gt.1.0e0) go to 190
      kz1=is*omegac*gamz0**2*(vz0-dsqrt(1.0e0-f11))
      kz2=is*omegac*gamz0**2*(vz0+dsqrt(1.0e0-f11))
      w1 =is*omegac*gamz0**2*(1.0e0-vz0*dsqrt(1.0e0-f11))
      w2 =is*omegac*gamz0**2*(1.0e0+vz0*dsqrt(1.0e0-f11))
  520 continue
      print 8, w1*wc,w2*wc,kz1*wc,kz2*wc
    8 format(' points of intersection:   (w1,w2)= (', 
     &1pe11.4,',',1pe11.4,
     &')',/,27x,'(kz1, kz2)= (',1pe11.4,',',1pe11.4,')'/)
      print 9
c   9 format(/' root(s) of dispersion relation:'/)
    9 format(/' root(s) of dispersion relation:'/'    f(hz)',4x
     +,'    w/wc',4x,' db/lamdac',2x,'   db/cm',4x,'    kzi',4x,      
     +'    kzr',4x,/,52x,'delwr',3x,'   value',4x,' i')    
c    +'    kzr',4x,3x,'delwr',3x,'   value',4x,' i')    
      iwmin=1
c
ccc15 specify frequency range 'wmin' and 'wmax', and the number of data
c     points 'iwmax' in-between (if iwmax.eq.1, kz is evaluated at w=wmin)
c
      wmin=1.0*wc 
      wmax=1.07*wg
      iwmax=51
c
      if(iwmax.eq.1) iplot=0
c
      do 100 iw=iwmin,iwmax
      if(iwmin.ne.iwmax) wstep=(wmax-wmin)/dfloat(iwmax-1)
      if(iwmin.eq.iwmax) wstep=0.0
      w=wmin+wstep*(iw-1)
      wwc=w/wc
      aw(iw)=w
      fhz=w*realc/(2.0*pi*rwcm)
      afhz(iw)=fhz
c
c guess root(s) to feed in subroutine 'muller'
c
      cdum1=dcmplx(w**2-xmn**2,0.0d0)
      cdum2=cdsqrt(cdum1)
      if(iw.eq.iwmin) ckz(1)=cdum2-ci*0.10e0/dfloat(itz)
      if(iw.eq.iwmin) ckz(2)=cdum2+ci*0.10e0/dfloat(itz)
      if(iw.eq.iwmin) ckz(3)=cdum2
      if(iw.eq.iwmin) ckz(4)=-cdum2
      if(tz.ne.0.0) ckz(1)=dcmplx(akzr0(iw,1),akzi0(iw,1))
      if(tz.ne.0.0) ckz(2)=dcmplx(akzr0(iw,2),akzi0(iw,2))
c
c the following statement gives instruction to skip further calculation
c of this data point under the specified conditions
c
      if(itz.ne.1.and.dabs(akzi0(iw,1)).lt.1.0d-3) go to 100
      icont1=0
      ep1=1.0e-4
      ep1=1.0e-6
      ep2=ep1
      imaxit=50
      if(tz.eq.0.0)
     &call muller(0, iroot,ckz,imaxit,ep1,ep2,cmask0,.false.)
      if(tz.ne.0.0)
     &call muller(0,iroot,ckz,imaxit,ep1,ep2,cmask1,.false.)
c
      iarg1=1
      iarg2=iroot
      do 101 ipnt=1,iroot
      kzr=ckz(ipnt)
      kzi=-ci*ckz(ipnt)
      if(kzi+1.0e-4) 42,41,41
   41 akzr0(iw,iarg2)=kzr
      akzi0(iw,iarg2)=kzi
      iarg2=iarg2-1
      go to 43
   42 akzr0(iw,iarg1)=kzr
      akzi0(iw,iarg1)=kzi
      akzi1(iw,iarg1)=-kzi
      if(iarg1.eq.1) akzr(iw,itz)=kzr
      if(iarg1.eq.1) akzi(iw,itz)=kzi
      iarg1=iarg1+1
   43 continue
      if(kzr.gt.kzrmax) kzrmax=kzr
      if(-1.0e0*kzi.gt.kzimax) kzimax=-1.0e0*kzi
      if(kzi.ge.(-1.0e-5)) go to 101
      if(kzr.lt.0.0d0) go to 101
      gain1=-54.66e0*kzi/xmn
      gain2=gain1/lamdac
      again1(iw,itz)=gain1
      again2(iw,itz)=gain2
      if(gain1.gt.g1max) g1max=gain1
      if(gain2.gt.g2max) g2max=gain2
      if(mod(iw-1, iskip+1).ne.0) go to 101
      vres=(w -is*omegac)/(kzr+1.0e-30)
      vphase=w/kzr
      delwr=w-kzr*vz0-is*omegac
      eta=1.25e0*gamma0*w *delwr/((gamma0-1.0e0)*(w **2-kzr**2))
      kzohm=xmn**2*delta/(2.0e0*dabs(kzr)+1.0e-40)
      etaw=eta*dabs(kzi)/(dabs(kzi)+kzohm+1.0e-40)
      etaohm=eta*kzohm/(dabs(kzi)+kzohm+1.0e-40)
      pwkw=etaw*ribamp*vbkv
      pohmkw=etaohm*ribamp*vbkv
      pohmpw=pohmkw/pwkw
      if(tz.eq.0.0) pole=1.0e-50
      cdum=ckz(ipnt)
      if(tz.eq.0.0) cvalue=cmask0(cdum)
      if(tz.ne.0.0) cvalue=cmask1(cdum)
      valuer=cvalue
      valuei=-ci*cvalue
      if(dabs(valuer).lt.1.0e-30) valuer=0.0e0
      if(dabs(valuei).lt.1.0e-30) valuei=0.0e0
      value=dsqrt(valuer**2+valuei**2)
c     if(ipnt.eq.1) print 12
c  12 format('  ')
c     print 11, w,fhz,wwc,kzr,delwr,kzi,value,icont1,eps,
c    &vphase
c  11 format(' w=',1pe11.4,'(',1pe11.4,' hz), w/wc=',1pe11.4,
c    &', kz=',1pe10.3,'(',1pe10.3,'), ',1pe10.3,
c    &',(',1pe9.2,',',i3,',',1pe10.3,',',1pe10.3,')')
c
c the following statements give instructions to skip printing data
c of secondary interest (as specified in the statements)
c
      if(kzi.gt.0.0) go to 101
      if(eta.lt.0.0) go to 101
      if(eta.gt.1.0) go to 101
      if(kzr.le.0.0) go to 101
c     print 6,  gain1,gain2,eta,etaw,pwkw,pohmkw
c   6 format(9x,
c    &' gain=',1pe10.3,' db per lamdac (',
c    &1pe10.3,' db/cm), eta=',1pe9.2,', etaw=',1pe9.2,
c    &', pwkw=',1pe9.2,', pohmkw=',1pe9.2)
      print 6,  fhz,wwc,gain1,gain2,kzi,kzr,delwr,value,
     +icont1           
    6 format(1x,
     &4(1pe11.4,1x),2(1pe10.3,1x),/,49x,2(1pe10.3,1x),i3)
c    &4(1pe11.4,1x),4(1pe10.3,1x),i3/)
  101 continue
  100 continue
      if(iplot.ne.1) go to 199
c 
c graphics
c  note:  for large plots, use subroutines 'bscale' and 'bplot'.
c         for small plots, use subroutines 'sscale' and 'splot'.
c 
      kzrmin=-kzrmax
      print 27
 27   format ('1',' ')
      call sscale(wmin,wmax, kzrmin,kzrmax)
      call splot(aw,akzr0(1,1),iwmax,'a   $')
      call splot(aw,akzr0(1,2),iwmax,'b   $')
      call splot(aw,akzr0(1,3),iwmax,'c   $')
      call splot(aw,akzr0(1,4),iwmax,'dd  kzr vs w$')
c
      call sscale(wmin,wmax, kzrmin,kzrmax)
      call splot(aw,akzr0(1,1),iwmax,'aa  unstable kzr vs w$')
c
      print 28
 28   format ('1',' ')
      call sscale(wmin,wmax,-kzimax,0.0d0)
      call splot(aw,akzi0(1,1),iwmax,'1   $')
      call splot(aw,akzi0(1,2),iwmax,'2   $')
      call splot(aw,akzi0(1,3),iwmax,'3   $')
      call splot(aw,akzi0(1,4),iwmax,'44  kzi vs w$')
c
      print 16
      go to 199
  190 print 18
   18 format(' **** no intersections ****'/////)
  199 continue
  800 continue
  801 continue
  200 continue
      if(iplot.ne.1) go to 599
      if(tz.eq.0.0) go to 599
c
      print 29
 29   format ('1',' ')
      call sscale(wmin,wmax,kzrmin,kzrmax)
      call splot(aw,akzr(1,1),iwmax,'1   $')
      call splot(aw,akzr(1,2),iwmax,'2   $')
      call splot(aw,akzr(1,3),iwmax,'3   $')
      call splot(aw,akzr(1,4),iwmax,'4   $')
      call splot(aw,akzr(1,5),iwmax,'5   $')
      call splot(aw,akzr(1,6),iwmax,'66  kzr vs w$')
c
      call sscale(wmin,wmax,-kzimax,0.0d0)
      call splot(aw,akzi(1,1),iwmax,'1   $')
      call splot(aw,akzi(1,2),iwmax,'2   $')
      call splot(aw,akzi(1,3),iwmax,'3   $')
      call splot(aw,akzi(1,4),iwmax,'4   $')
      call splot(aw,akzi(1,5),iwmax,'5   $')
      call splot(aw,akzi(1,6),iwmax,'66  kzi vs w$')
c
  599 continue
c
      print 30
 30   format ('1',' ')
      call sscale(afhz(iwmin),afhz(iwmax),0.0d0,g1max)
      call splot(afhz,again1(1,1),iwmax,'1   $')
      call splot(afhz,again1(1,2),iwmax,'2   $')
      call splot(afhz,again1(1,3),iwmax,'3   $')
      call splot(afhz,again1(1,4),iwmax,'4   $')
      call splot(afhz,again1(1,5),iwmax,'5   $')
      call splot(afhz,again1(1,6),iwmax,'66  gain (in db per cutoff wave
     &length) vs frequency (in hz)$')
c
      call sscale(afhz(iwmin),afhz(iwmax),0.0d0,g2max)
      call splot(afhz,again2(1,1),iwmax,'1   $')
      call splot(afhz,again2(1,2),iwmax,'2   $')
      call splot(afhz,again2(1,3),iwmax,'3   $')
      call splot(afhz,again2(1,4),iwmax,'4   $')
      call splot(afhz,again2(1,5),iwmax,'5   $')
      call splot(afhz,again2(1,6),iwmax,'66  gain (in db/cm) vs frequenc
     &y (in hz)$')
c
      print 16
   16 format(/////////)
  900 continue
  600 continue
      end
      function cmask0(ckz)
      implicit  real*8 (a,b,d-h,j-z), complex*16 (c)
      common ci,pi,e,e2,realc,realc2,me,
     &ppmax,ppmin,pzmax,pzmin,divpp,divpz,
     &nu,r1,r2,rc0,omegae, w,eps,a,pp0,pz0,p0,gamma0,delpp,delpz,
     &ppc,xmn,kmn,delta,hsm,tsm,usm,
     &is,im,in,ippmax,ipzmax,icont1
c
      fac1=1.0+dfloat(im*im)*w*w/((xmn**2-dfloat(im*im))*xmn**2)
c
      c1=gamma0*w-ckz*pz0-is*omegae
      c2=( w*w-ckz**2)*hsm*pp0**2/gamma0
c     c3=(w-ckz*pz0/gamma0)*qsm*c1
      c3=((w-ckz*pz0/gamma0)*tsm-xmn*pp0*usm/gamma0)*c1
      c4=( w*w-ckz**2-xmn**2*(1.0e0-delta*(1.0e0+ci)*fac1))*c1*c1
      c5=-4.0e0*nu*(c2-c3)/kmn
      eps=cdabs(c3/c2)
c
      cmask0=c4-c5
c
      icont1=icont1+1
      return
      end
      function cmask1(ckz)
      implicit  real*8 (a,b,d-h,j-z), complex*16 (c)
      common ci,pi,e,e2,realc,realc2,me,
     &ppmax,ppmin,pzmax,pzmin,divpp,divpz,
     &nu,r1,r2,rc0,omegae, w,eps,a,pp0,pz0,p0,gamma0,delpp,delpz,
     &ppc,xmn,kmn,delta,hsm,tsm,usm,
     &is,im,in,ippmax,ipzmax,icont1
      csum1=dcmplx(0.0e0,0.0e0)
      fac1=1.0+dfloat(im*im)*w*w/((xmn**2-dfloat(im*im))*xmn**2)
      do 100 iz=1,ipzmax
      pz=pzmin+divpz*(2*iz-1)
      fac2=dexp(-(pz-pz0)**2/(2.0*delpz**2+1.0e-20))
      do 100 ip=1,ippmax
      pp=dsqrt(p0**2-pz**2)
      gamma=dsqrt(1.0e0+pp**2+pz**2)
      rl=pp/omegae
      call couple(is,im,xmn*rc0,xmn*rl,hsm,tsm,usm)
      c1=gamma*w-ckz*pz-is*omegae
      c2=(w**2-ckz**2)*hsm*pp**2/(gamma*c1**2)
c     c3=(w-ckz*pz/gamma)*qsm/c1
      c3=((w-ckz*pz/gamma)*tsm-xmn*pp*usm/gamma)/c1
      csum1=csum1+a*fac2* (-c2+c3)
  100 continue
      c4=w**2-ckz**2-xmn**2*(1.0e0-delta*(1.0e0+ci)*fac1)
      c5=8.0e0*pi*nu*csum1/kmn
      eps=cdabs(c5/ w**2)
c
      cmask1=c4-c5
c
      icont1=icont1+1
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
c ref. ralston & wilf 'mathematical methods for digital computers',p.117
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

