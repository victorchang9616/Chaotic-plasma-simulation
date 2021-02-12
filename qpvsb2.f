c     program qpvsb2 (adapted for qpvsb calculation for a single cavity)
c                   note: compared with program qpvsb, this program can
c                         treat a beam with velocity spread. it can also
c                         calculate the hot cavity resonant frequency,
c                         hence the beam current and cavity q need be 
c                         specified.
c
c the hot cavity resonant frequency is calculated from the linear
c solution of eq. (18) of k. r. chu et al., ieee trans. plasma science
c vol.13, pp. 424-434 (1985)         
c
c date of present  version: March 1, 1999
c
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      dimension axmn(9,8),alcav(20),arwcav(20),aqcav(20),
     &iin(20),aw0(20)
      dimension llamcv(20)
      dimension adshot(201,9),ab0kg(201),aqpbkw(201,9),adwhw0(201,9)
      dimension y1(1,201,21,1,1),y2(1,201,21,1,1),
     &y3(1,201,21,1,1),y4(1,201,21,1,1),y5(1,201,21,1,1),
     &y6(1,201,21,1,1),weight(1,201,21,1,1)
c dimensions of y1, etc. should be > or =(itmax,ipzmax,ircmax,ipsimx,iphimx)
      dimension adww0f(201),avbz02(201),abz02(201)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      common/ebeam/weight,y1,y2,y3,y4,y5,y6
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
c constants
c
      ci=dcmplx(0.0d0,1.0d0)
      pi=3.1415927d0
      realc=2.99792d10
c
c give instructions on plot and write(6,(1:yes, 0:no)
c note: set ipmore to 1 if more detailed print-out is desired. in case of
c       large data run, set ipmore to 0 to avoid large print-out
c
      ipmore=0
      iplot=1
c
c specify no of cavities (for this adapted version, always set icavmx=1)
c
      icavmx=1
c
c specify the cold resonant frequency of the 1st cavity in ghz (if only the
c normalized results are desired, set fghz at any value, e.g. 10.0)
c note: quantities are normalized to the radius of the first cavity (hence
c       normalized radius of first cavity=1.0). see k. r. chu, phys. of
c       fluids, vol. 21(1978), p.2358 for the normalization scheme
c
      fghz=33.1307
c
c specify the cavity no at and beyond which numerical orbit integrations
c (valid in both small and large signal regimes) will be used. analytical
c small signal calculaions will be used for preceeding cavities.
c
      icavnm=1
c
c
c specify im, in, il, and is
c assume im il, and is are the same for all cavities.
c
      im=1
      il=1
      is=1
      iin(1)=1
c
c specify q of each cavity
c
      aqcav(1)=300.0
c
c specify lengths of cavities and drift regions in terms of
c ratios with respect to the cold wavelength of the first cavity
c
      llamcv(1)=11.145
c
c calculate normalized lengths of cavities and drift regions, and normalized
c cold resonant frequency (w0fst) of first cavity
c
      arwcav(1)=1.0
      infst=iin(1)
      xmnfst=axmn(im+1,infst)
      alcav(1)=(pi/xmnfst)*dsqrt(4.0*llamcv(1)**2-dfloat(il)**2)
      kzfst=dfloat(il)*pi/alcav(1)
      w0fst=dsqrt(xmnfst**2/arwcav(1)**2+kzfst**2)
      aw0(1)=w0fst
      rw1cm=w0fst*realc/(2.0e9*pi*fghz)
c
c calculate normalized radius of each cavity
c
      do 10 i=1,icavmx
      in=iin(i)
      l=alcav(i)
      xmn=axmn(iabs(im)+1,in)
      kz=dfloat(il)*pi/l
      arwcav(i)=xmn/dsqrt(aw0(i)**2-kz**2)
   10 continue
      write(6,33) icavmx,fghz,rw1cm
   33 format(/" total no of cavities=",i3,
     &", cold resonant freq. of 1st cavity=",1pd12.5,
     &" ghz"/" radius of 1st cavity=",1pd12.5," cm")
      write(6,15)is,im,il
   15 format(/" is=",i2,", im=",i2,", il=",i2)
      write(6,25)(iin(i),i=1,icavmx)
   25 format(/"    in=",6(5x,i2,",",5x))
      write(6,26)(alcav(i),i=1,icavmx)
   26 format("     l=",6(1pd12.5,","))
      write(6,28)(arwcav(i),i=1,icavmx)
   28 format("    rw=",6(1pd12.5,","))
      write(6,29)(aqcav(i),i=1,icavmx)
   29 format("     q=",6(1pd12.5,","))
      write(6,30)(aw0(i),i=1,icavmx)
   30 format("    w0=",6(1pd12.5,","))
      write(6,17)
   17 format("-------------------------------------------------------")
c
c specify beam parameters
c assume all electrons have the same gamma.
c
      do 2000 ircavg=1,1
      drc=0.00
      rcavg=0.35
      rc1=rcavg-0.5*drc
      rc2=rcavg+0.5*drc
c
      do 2000 ivbkv=1,1
      vbkv=100.0
c
      do 2000 ialpha=1,1
      alpha0=0.85
      angle=datan2(alpha0,1.0d0)*180.0/pi
c
c the following statements, when activated, will override the previously
c specified alpha0 and angle
c
c     angle=63.435
c     alpha0=dsin(angle*pi/180.0)/dcos(angle*pi/180.0)
c
      do 2000 irib=1,1
      ribamp=1.0d-7
c
c specify normalized rf field amplitude in the first cavity at w=w0hot
c (qp, etc. will be independent of cbz00)
c
      cbz00=dcmplx(1.0e-3,0.0)
c
      do 60 i1=1,201
      ab0kg(i1)=0.0d0
      do 60 i2=1,9
      adshot(i1,i2)=0.0d0
      aqpbkw(i1,i2)=-1.0d10
      adwhw0(i1,i2)=-1.0d10
   60 continue
      dsmax=-1.0d10
      dsmin=1.0d10
      b0kgmx=-1.0d10
      b0kgmn=1.0d10
      qpbmax=-1.0d10
      dwhmax=-1.0d10
      dwhmin=1.0d10
c  
      itzmax=3
      do 1500 itz=1,itzmax
      tz=0.10e0*(itz-1)
c
c calculate other beam parameters.
c
      dz=2.5
      gamma0=1.0e0+vbkv/511.0e0
      pbkw=vbkv*ribamp
      write(6,1)vbkv,ribamp,pbkw,gamma0,alpha0,angle,
     &rcavg,drc,rc1,rc2
    1 format(/" vb=",1pd12.5," kv, ib=",1pd10.3," amp, pb=",
     &1pd10.3," kw",/" gamma0=",1pd13.6,
     &", alpha=",1pd10.3,"(",1pd11.4,
     &" degree)"/" rcavg=",1pd10.3,", drc=",1pd10.3,", rc1=",1pd10.3,
     &", rc2=",1pd10.3/)
      theta0=datan2(alpha0,1.0d0)
      p0=dsqrt(gamma0**2-1.0e0)
      pp0=p0*dsin(theta0)
      pz0=p0*dcos(theta0)
      delpz=tz*pz0
      v0=p0/gamma0
      vp0=v0*dsin(theta0)
      vz0=v0*dcos(theta0)
      pzmin=dmax1(0.01d0*delpz,pz0-dz*delpz)
      pzmax=dmin1(0.99d0*p0,pz0+dz*delpz)
      ppmax=dsqrt(p0**2-pzmin**2)
      ppmin=dsqrt(p0**2-pzmax**2)
      vpmin=ppmin/gamma0
      vpmax=ppmax/gamma0
      vzmin=pzmin/gamma0
      vzmax=pzmax/gamma0
      write(6,4) p0,pp0,pz0,delpz,tz,
     &ppmin,ppmax,pzmin,pzmax,
     &v0,vp0,vz0,vpmin,vpmax,vzmin,vzmax
    4 format(" p0=",1pd10.3,", pp0=",1pd10.3,", pz0=",1pd10.3,
     &", delpz=",1pd10.3/" tz=delpz/pz0=",1pd10.3/
     &" ppmin=",1pd10.3,", ppmax=",1pd10.3/
     &" pzmin=",1pd10.3,", pzmax=",1pd10.3/
     &" v0=",1pd10.3,", vp0=",1pd10.3,", vz0=",1pd10.3/
     &" vpmin=",1pd10.3,", vpmax=",1pd10.3/
     &" vzmin=",1pd10.3,", vzmax=",1pd10.3/)
c
c specify no of electrons
c
      itmax=1
      ipzmax=4*int(tz*100.0)+1
      ircmax=6
      ipsimx=1
      iphimx=1
      if(im.eq.0) ipsimx=1
      if(tz.eq.0.0) ipzmax=1
      if(drc.eq.0.0) ircmax=1
      itotal=itmax*ipzmax*ircmax*ipsimx*iphimx
      write(6,2)  itmax,ipzmax,ircmax,ipsimx,iphimx,itotal
    2 format(" itmax=",i3,", ipzmax=",i3,", ircmax=",i3,
     &", ipsimx=",i3,", iphimx=",i3," itotal=",i8)
c
c specify normalized omegae by specifying ds(with respect to w0fst) of first
c cavity
c note: if idsmax>1, iwmax should be set to 1 and vice versa (i.e. b-field 
c       and frequency can not be scanned at the same time).
c
      idsmax=31
      do 1500 ids=1,idsmax
      dscold=-10.0d0+1.0*(ids-1)+1.0e-4
      in=iin(1)
      l=alcav(1)
      rw=arwcav(1)
      q=aqcav(1)
      kz=dfloat(il)*pi/l
      xmn=axmn(iabs(im)+1,in)
      w0fst=aw0(1)
      tau=l/vz0
      omegae=(w0fst-kz*vz0-dscold/tau)*gamma0/dfloat(is)
c
c the following statement, when activated, will override the previously
c specified dscold and omegae
c
c     omegae=8.6867*(1.0+0.003*(ids-4))
c
      omegac=omegae/gamma0
      dscold=(w0fst-kz*vz0-dfloat(is)*omegac)*tau
      rl0=pp0/omegae
      rtobw0=omegac/w0fst
      b0kg=1.7045*omegae/rw1cm
      ab0kg(ids)=b0kg
      if(b0kg.gt.b0kgmx) b0kgmx=b0kg
      if(b0kg.lt.b0kgmn) b0kgmn=b0kg
      write(6,8) w0fst,dscold,omegae,b0kg,omegac,rtobw0,rl0
    8 format(" w0(1)=",1pd12.5,
     &", ds(1) at (w=w0(1))=",1pd10.3/
     &" omegae=",1pd12.5,"(",1pd12.5," kg), omegac=",1pd12.5/
     &" omegac/w0fst=",1pd12.5,", rl0=",1pd10.3)
      if((rl0+rc2).ge.1.0e0) write(6,5)
      if((rl0+rc2).ge.1.0e0) stop
      if(rc1.le.0.0) write(6,32)
      if(rc1.le.0.0) stop
   32 format(//" ***** rc1.le.0.0 *****"//)
    5 format(//" ***** beam hits wall *****"//)
c
c calculate w0hot and qhot
c
      call sety(w0fst,weight,y1,y2,y3,y4,y5,y6)
      call hotw01(w0hot,value)
      dwh=w0hot-w0fst
      dwhw0=dwh/w0fst
      dum21=bj(im,xmn)**2*(1.0e0-dfloat(im*im)/xmn**2)
      call sety(w0hot,weight,y1,y2,y3,y4,y5,y6)
      call anavg(w0hot,cavg10,cavg20)
      dum22=cavg10
      poutkw=5.46376e5*l*rw**4*w0hot*w0fst**2*cdabs(cbz00)**2*dum21/
     &(q*xmn**2)
      pinkw=(1.0e0/8.0e0)*pbkw*l*rw**2*w0hot*cdabs(cbz00)**2*dum22
     &/((gamma0-1.0e0)*xmn**2)
      cdum23=1.0e0/q-2.2878e-7*pbkw*cavg10
     &/(w0fst**2*rw**2*(gamma0-1.0e0)*dum21)
      dum23=cdum23
      qhot=1.0/dum23
      cnum=(w0hot/w0fst)**2-1.0+ci*cdum23
      write(6,16) w0hot,value,icont1,dwhw0,cbz00,poutkw,pinkw,q,qhot
   16 format(" w0hot(1)=",1pd10.3,"(",1pd9.2,",",i4,
     &"), (w0hot(1)-w0(1))/w0(1)=",1pd10.3/
     &" cbz00 (=cbz0(1) at w=w0hot(1)) =",1pd10.3,",",1pd10.3/
     &" poutkw(1)=",1pd10.3,", pinkw(1)=",1pd10.3/
     &" qcold(1)=",1pd10.3,", qhot=",1pd10.3/)
c
c calculate start oscilation conditions of all cavities (assume w0hot of
c all cavities are the same as that of the first cavity).
c
      do 200 icav=1,icavmx
      in=iin(icav)
      l=alcav(icav)
      rw=arwcav(icav)
      q=aqcav(icav)
      kz=dfloat(il)*pi/l
      xmn=axmn(iabs(im)+1,in)
      w0=aw0(icav)
      call anavg0(w0hot,cavg1)
      dum21=bj(im,xmn)**2*(1.0e0-dfloat(im*im)/xmn**2)
      dum26=cavg1
      qpbkw=4.3710e6*w0**2*rw**2*(gamma0-1.0e0)*dum21/dum26
      dshot =(w0hot-kz*vz0-dfloat(is)*omegac)*l/vz0
      dsphot=(w0hot+kz*vz0-dfloat(is)*omegac)*l/vz0
      adshot(ids,itz)=dshot
      aqpbkw(ids,itz)=qpbkw
      adwhw0(ids,itz)=dwhw0
      if(dshot.gt.dsmax) dsmax=dshot
      if(dshot.lt.dsmin) dsmin=dshot
      if(qpbkw.gt.qpbmax) qpbmax=qpbkw
      if(dwhw0.gt.dwhmax) dwhmax=dwhw0
      if(dwhw0.lt.dwhmin) dwhmin=dwhw0
      write(6,9) dshot,dsphot
    9 format(" ds(1) (at w=w0hot)=",1pd10.3,", dsp(1) (at w=w0hot)=",
     &1pd10.3)
      write(6,18)icav,qpbkw
   18 format(" icav=",i3,", qpbkw=",1pd10.3)
  200 continue
c
c set default values of arrays
c
      do 40 i1=1,201
      adww0f(i1)=0.0
      abz02(i1)=0.0
      avbz02(i1)=0.0
   40 continue
c
c specify normalized driver frequency (w) and number of points desired (iwmax)
c
      dwmin=-2.0e-3
      dwmax=3.0e-3
      iwmax=1
      if(iwmax.eq.1) iplot=0
      wmin=w0fst*(1.0+dwmin)
      wmax=w0fst*(1.0+dwmax)
      if(iwmax.gt.1) wstep=(wmax-wmin)/dfloat(iwmax-1)
      gainmx=-1.0e10
      bz02mx=-1.0e10
c
c scan spetral responses of first cavity (if iwmax > 1)
c
      do 1000 iw=1,iwmax
      if(iwmax.gt.1) w=wmin+wstep*(iw-1)
      if(iwmax.eq.1) w=w0hot
      ww0fst=w/w0fst
      dww0f=ww0fst-1.0
      adww0f(iw)=dww0f
c
      if(ipmore.eq.1)
     &write(6,6) dww0f,ww0fst
    6 format(/" (w-w0(1))/w0(1)= ",1pd10.3," (w/w0(1)=",1pd13.6,")")
c
c first cavity calculations
c
      in=iin(1)
      l=alcav(1)
      rw=arwcav(1)
      q=aqcav(1)
      kz=dfloat(il)*pi/l
      xmn=axmn(iabs(im)+1,in)
      w0=aw0(1)
      ww0=w/w0
      ds =(w-kz*vz0-dfloat(is)*omegac)*l/vz0
      dsp=(w+kz*vz0-dfloat(is)*omegac)*l/vz0
      if(ipmore.eq.1)
     &write(6,11) in,l,rw,q,xmn,w0,ww0,kz,ds,dsp
   11 format(" in=",i3,", l=",1pd10.3,", rw=",1pd13.6,", q=",1pd10.3/
     &" xmn=",1pd12.5,", w0=",1pd13.6,", w/w0=",1pd13.6,", kz=",1pd10.3/
     &" ds=",1pd10.3,", dsp=",1pd10.3)
      dum25=ww0**2-1.0e0
      dum21=bj(im,xmn)**2*(1.0e0-dfloat(im*im)/xmn**2)
      call sety(w,weight,y1,y2,y3,y4,y5,y6)
      call anavg(w,cavg1,cavg2)
      dum26=cavg1
      cdum24=1.0e0/q-2.2878e-7*pbkw*cavg1
     &/(w0**2*rw**2*(gamma0-1.0e0)*dum21)
      cbz0=cnum*cbz00/(dum25+ci*cdum24)
      bz02=cdabs(cbz0/cbz00)**2
      cvbz0=ci*cbz00/(q*(dum25+ci/q))
      vbz02=cdabs(cvbz0/cbz00)**2
      avbz02(iw)=vbz02
      abz02(iw)=bz02
c "v" in cvbz0 and vbz02 implies "under vacuum condition"
      if(vbz02.gt.bz02mx) bz02mx=vbz02
      if(bz02.gt.bz02mx) bz02mx=bz02
      if(ipmore.eq.1)
     &write(6,12) cbz0,bz02
   12 format(" cbz0=",1pd10.3,",",1pd10.3,", !cbz0/cbz00!**2=",1pd12.5)
      if(ipmore.eq.1)
     &write(6,13) cavg20,cavg2
   13 format(" (check: the following numbers should be close to zero:"/
     &1x,4(1pd10.3,","),")")
      write(6,17)
 1000 continue
c
      if(iplot.ne.1) go to 86
      call sscale(adww0f(1),adww0f(iwmax),0.0,bz02mx)
      call splot(adww0f,avbz02,iwmax,'v   $')
      call splot(adww0f,abz02,iwmax,'bb  bz0**2 vs (w-w0)/w0 in the
     & presence and absence of beam$')
      write(6,88)
   88 format(//)
   86 continue
 1500 continue
      qpbmax=1.0d4
      if(idsmax.eq.1) go to 1501
      call sscale(dsmin,dsmax,0.0d0,qpbmax)
      call splot(adshot(1,1),aqpbkw(1,1),idsmax,'1   $')
      call splot(adshot(1,2),aqpbkw(1,2),idsmax,'2   $')
      call splot(adshot(1,3),aqpbkw(1,3),idsmax,'3   $')
      call splot(adshot(1,4),aqpbkw(1,4),idsmax,'4   $')
      call splot(adshot(1,5),aqpbkw(1,5),idsmax,'5   $')
      call splot(adshot(1,6),aqpbkw(1,6),idsmax,'6   $')
      call splot(adshot(1,7),aqpbkw(1,7),idsmax,'7   $')
      call splot(adshot(1,8),aqpbkw(1,8),idsmax,'8   $')
      call splot(adshot(1,9),aqpbkw(1,9),idsmax,'99  qpbkw vs dshot
     &$')
      call sscale(b0kgmn,b0kgmx,0.0d0,qpbmax)
      call splot(ab0kg,aqpbkw(1,1),idsmax,'1   $')
      call splot(ab0kg,aqpbkw(1,2),idsmax,'2   $')
      call splot(ab0kg,aqpbkw(1,3),idsmax,'3   $')
      call splot(ab0kg,aqpbkw(1,4),idsmax,'4   $')
      call splot(ab0kg,aqpbkw(1,5),idsmax,'5   $')
      call splot(ab0kg,aqpbkw(1,6),idsmax,'6   $')
      call splot(ab0kg,aqpbkw(1,7),idsmax,'7   $')
      call splot(ab0kg,aqpbkw(1,8),idsmax,'8   $')
      call splot(ab0kg,aqpbkw(1,9),idsmax,'99  qpbkw vs b0kg$')

      call sscale(dsmin,dsmax,dwhmin,dwhmax)
      call splot(adshot(1,1),adwhw0(1,1),idsmax,'1   $')
      call splot(adshot(1,2),adwhw0(1,2),idsmax,'2   $')
      call splot(adshot(1,3),adwhw0(1,3),idsmax,'3   $')
      call splot(adshot(1,4),adwhw0(1,4),idsmax,'4   $')
      call splot(adshot(1,5),adwhw0(1,5),idsmax,'5   $')
      call splot(adshot(1,6),adwhw0(1,6),idsmax,'6   $')
      call splot(adshot(1,7),adwhw0(1,7),idsmax,'7   $')
      call splot(adshot(1,8),adwhw0(1,8),idsmax,'8   $')
      call splot(adshot(1,9),adwhw0(1,9),idsmax,'99  (w0hot-w0)/w0 vs
     & dshot$')
      call sscale(b0kgmn,b0kgmx,dwhmin,dwhmax)
      call splot(ab0kg,adwhw0(1,1),idsmax,'1   $')
      call splot(ab0kg,adwhw0(1,2),idsmax,'2   $')
      call splot(ab0kg,adwhw0(1,3),idsmax,'3   $')
      call splot(ab0kg,adwhw0(1,4),idsmax,'4   $')
      call splot(ab0kg,adwhw0(1,5),idsmax,'5   $')
      call splot(ab0kg,adwhw0(1,6),idsmax,'6   $')
      call splot(ab0kg,adwhw0(1,7),idsmax,'7   $')
      call splot(ab0kg,adwhw0(1,8),idsmax,'8   $')
      call splot(ab0kg,adwhw0(1,9),idsmax,'99  (w0hot-w0)/w0 vs b0kg
     &$')
 1501 continue
 2000 continue
      stop
      end
c *********************************************************************
      subroutine anavg(w,cavg1,cavg2)
c input w, weight,  and y1...y6, return cavg1 and cavg2.
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      dimension y1(1,201,21,1,1),y2(1,201,21,1,1),
     &y3(1,201,21,1,1),y4(1,201,21,1,1),y5(1,201,21,1,1),
     &y6(1,201,21,1,1),weight(1,201,21,1,1)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      common/ebeam/weight,y1,y2,y3,y4,y5,y6
      dum1=w*w-kz*kz
      cavg1=dcmplx(0.0e0,0.0e0)
      cavg2=dcmplx(0.0e0,0.0e0)
      do 100 it=1,itmax
      do 100 ipz=1,ipzmax
      do 100 irc=1,ircmax
      do 100 ipsi=1,ipsimx
      do 100 iphi=1,iphimx
      rc  =y1(it,ipz,irc,ipsi,iphi)
      psi =y2(it,ipz,irc,ipsi,iphi)
      pp  =y3(it,ipz,irc,ipsi,iphi)
      phi =y4(it,ipz,irc,ipsi,iphi)
      pz  =y5(it,ipz,irc,ipsi,iphi)
      time=y6(it,ipz,irc,ipsi,iphi)
      gamma=dsqrt(1.0e0+pp*pp+pz*pz)
      betap=pp/gamma
      betaz=pz/gamma
      tau=l/betaz
      omegac=omegae/gamma
      ds =(w-kz*betaz-dfloat(is)*omegac)*tau
      dsp=(w+kz*betaz-dfloat(is)*omegac)*tau
      if(dabs(ds)-1.0e-3) 11,11,12
   11 ds=1.0e-3
      dsp=ds+2.0*kz*betaz*tau
      go to 14
   12 if(dabs(dsp)-1.0e-3) 13,13,14
   13 dsp=1.0e-3
      ds=dsp-2.0*kz*betaz*tau
   14 continue
      rl=pp/omegae
      xmnrc=xmn*rc/rw
      xmnrl=xmn*rl/rw
      call couple(is,im,xmnrc,xmnrl,hsm,tsm,usm)
      cdum2=cdexp(ci*ds)
      cdum3=cdexp(ci*dsp)
      dum4=(betap*tau)**2*hsm/gamma
      dum5=kz/l
      dum6=tau/gamma
      dum7=xmn*betap*usm/rw
      ca11=dum4*(dum1*(ci*ds*(cdum2+1.0e0)-2.0e0*(cdum2-1.0e0))
     &-dum5*ds*(ci*ds*cdum2-cdum2+1.0e0))/ds**3
      ca12=-dum6*(ci*ds-cdum2+1.0e0)*((w-kz*betaz)*tsm-dum7)/ds**2
      ca21=-dum4*(dum1*(ci*ds*cdum3-(1.0e0+ds/dsp)*(cdum3-1.0e0))
     &-dum5*ds*(ci*ds*cdum3-(ds/dsp)*(cdum3-1.0e0)))/(ds**2*dsp)
      ca22=-dum6*(cdum3-1.0e0)*((w-kz*betaz)*tsm-dum7)/(ds*dsp)
      ca31=dum4*(dum1*(ci*dsp*(cdum3+1.0e0)-2.0e0*(cdum3-1.0e0))
     &+dum5*dsp*(ci*dsp*cdum3-cdum3+1.0e0))/dsp**3
      ca32=-dum6*(ci*dsp-cdum3+1.0e0)*((w+kz*betaz)*tsm-dum7)/dsp**2
      ca41=-dum4*(dum1*(ci*dsp*cdum2-(1.0e0+dsp/ds)*(cdum2-1.0e0))
     &+dum5*dsp*(ci*dsp*cdum2-(dsp/ds)*(cdum2-1.0e0)))/(dsp**2*ds)
      ca42=-dum6*(cdum2-1.0e0)*((w+kz*betaz)*tsm-dum7)/(ds*dsp)
      cavg1=cavg1+(ca11+ca12+ca21+ca22+ca31+ca32+ca41+ca42)
     &*weight(it,ipz,irc,ipsi,iphi)/betaz
      bs0=w*time-dfloat(is)*phi+dfloat(is-im)*psi
     &+(dfloat(im)-dfloat(is)/2.0e0)*pi
      cdum11=(cdum2-1.0e0)/ds-(cdum3-1.0e0)/dsp
      dum8=bj(is-im,xmnrc)*bjp(is,xmnrl)
      cavg2=cavg2+betap* dum8*cdum11*cdexp(ci*bs0)
     &*weight(it,ipz,irc,ipsi,iphi)/betaz
  100 continue
      return
      end
c *********************************************************************
      subroutine anavg0(w,cavg1)
c this program is valid only for unbunched input beam parameters (weight,
c y1, ... y6,etc.)
c input w, weight,  and y1...y6, return cavg1.
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      dimension y1(1,201,21,1,1),y2(1,201,21,1,1),
     &y3(1,201,21,1,1),y4(1,201,21,1,1),y5(1,201,21,1,1),
     &y6(1,201,21,1,1),weight(1,201,21,1,1)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      common/ebeam/weight,y1,y2,y3,y4,y5,y6
      dum1=w*w-kz*kz
      cavg1=dcmplx(0.0e0,0.0e0)
      it=1
      ipsi=1
      iphi=1
      do 100 ipz=1,ipzmax
      do 100 irc=1,ircmax
      rc  =y1(it,ipz,irc,ipsi,iphi)
      pp  =y3(it,ipz,irc,ipsi,iphi)
      pz  =y5(it,ipz,irc,ipsi,iphi)
      gamma=dsqrt(1.0e0+pp*pp+pz*pz)
      betap=pp/gamma
      betaz=pz/gamma
      tau=l/betaz
      omegac=omegae/gamma
      ds =(w-kz*betaz-dfloat(is)*omegac)*tau
      dsp=(w+kz*betaz-dfloat(is)*omegac)*tau
      if(dabs(ds)-1.0e-3) 11,11,12
   11 ds=1.0e-3
      dsp=ds+2.0*kz*betaz*tau
      go to 14
   12 if(dabs(dsp)-1.0e-3) 13,13,14
   13 dsp=1.0e-3
      ds=dsp-2.0*kz*betaz*tau
   14 continue
      rl=pp/omegae
      xmnrc=xmn*rc/rw
      xmnrl=xmn*rl/rw
      call couple(is,im,xmnrc,xmnrl,hsm,tsm,usm)
      cdum2=cdexp(ci*ds)
      cdum3=cdexp(ci*dsp)
      dum4=(betap*tau)**2*hsm/gamma
      dum5=kz/l
      dum6=tau/gamma
      dum7=xmn*betap*usm/rw
      ca11=dum4*(dum1*(ci*ds*(cdum2+1.0e0)-2.0e0*(cdum2-1.0e0))
     &-dum5*ds*(ci*ds*cdum2-cdum2+1.0e0))/ds**3
      ca12=-dum6*(ci*ds-cdum2+1.0e0)*((w-kz*betaz)*tsm-dum7)/ds**2
      ca21=-dum4*(dum1*(ci*ds*cdum3-(1.0e0+ds/dsp)*(cdum3-1.0e0))
     &-dum5*ds*(ci*ds*cdum3-(ds/dsp)*(cdum3-1.0e0)))/(ds**2*dsp)
      ca22=-dum6*(cdum3-1.0e0)*((w-kz*betaz)*tsm-dum7)/(ds*dsp)
      ca31=dum4*(dum1*(ci*dsp*(cdum3+1.0e0)-2.0e0*(cdum3-1.0e0))
     &+dum5*dsp*(ci*dsp*cdum3-cdum3+1.0e0))/dsp**3
      ca32=-dum6*(ci*dsp-cdum3+1.0e0)*((w+kz*betaz)*tsm-dum7)/dsp**2
      ca41=-dum4*(dum1*(ci*dsp*cdum2-(1.0e0+dsp/ds)*(cdum2-1.0e0))
     &+dum5*dsp*(ci*dsp*cdum2-(dsp/ds)*(cdum2-1.0e0)))/(dsp**2*ds)
      ca42=-dum6*(cdum2-1.0e0)*((w+kz*betaz)*tsm-dum7)/(ds*dsp)
      cavg1=cavg1+(ca11+ca12+ca21+ca22+ca31+ca32+ca41+ca42)
     &*weight(it,ipz,irc,ipsi,iphi)/betaz
  100 continue
      cavg1=cavg1*dfloat(itmax*ipsimx*iphimx)
      return
      end
c *********************************************************************
      subroutine ancord(cbz0,w,weight,y1,y2,y3,y4,y5,y6)
c input cbz0, w, weight, and y1,...y6, return updated y1,... y6.
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      dimension y1(1,201,21,1,1),y2(1,201,21,1,1),
     &y3(1,201,21,1,1),y4(1,201,21,1,1),y5(1,201,21,1,1),
     &y6(1,201,21,1,1),weight(1,201,21,1,1)
      dimension t(1000)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      do 100 it=1,itmax
      do 100 ipz=1,ipzmax
      do 100 irc=1,ircmax
      do 100 ipsi=1,ipsimx
      do 100 iphi=1,iphimx
      rc  =y1(it,ipz,irc,ipsi,iphi)
      psi =y2(it,ipz,irc,ipsi,iphi)
      pp  =y3(it,ipz,irc,ipsi,iphi)
      phi =y4(it,ipz,irc,ipsi,iphi)
      pz  =y5(it,ipz,irc,ipsi,iphi)
      time=y6(it,ipz,irc,ipsi,iphi)
      gamma=dsqrt(1.0e0+pp*pp+pz*pz)
      y4(it,ipz,irc,ipsi,iphi)=y4(it,ipz,irc,ipsi,iphi)+omegae*l/pz
      y6(it,ipz,irc,ipsi,iphi)=y6(it,ipz,irc,ipsi,iphi)+gamma*l/pz
      if(cdabs(cbz0).eq.0.0) go to 100
      betap=pp/gamma
      betaz=pz/gamma
      tau=l/betaz
      omegac=omegae/gamma
      eps =w-kz*betaz-dfloat(is)*omegac
      epsp=w+kz*betaz-dfloat(is)*omegac
      ds=eps*tau
      dsp=epsp*tau
      rl=pp/omegae
      xmnrc=xmn*rc/rw
      xmnrl=xmn*rl/rw
      bs0=w*time-dfloat(is)*phi+dfloat(is-im)*psi
     &+(dfloat(im)-dfloat(is)/2.0e0)*pi
      cdum8=cdexp(-ci*ds)
      cdum9=cdexp(-ci*bs0)
      cdum10=cdexp(-ci*dsp)
      cdum1=cbz0*rw*cdum9/(2.0*xmn)
      ism=is-im
      call bes(iabs(ism)+2,xmnrc,0,result,t)
      if(ism.ge.0) dum11=t(ism+1)
      if(ism.lt.0) dum11=(-1)**ism*t(-ism+1)
      if((ism-1).ge.0) dum4=t(ism)
      if((ism-1).lt.0) dum4=(-1)**(ism-1)*t(2-ism)
      if((ism+1).ge.0) dum6=t(ism+2)
      if((ism+1).lt.0) dum6=(-1)**(ism+1)*t(-ism)
      call bes(iabs(is)+2,xmnrl,0,result,t)
      if(is.ge.0) dum3=t(is+1)
      if(is.lt.0) dum3=(-1)**is*t(-is+1)
      if((is-1).ge.0) dum5=t(is)
      if((is-1).lt.0) dum5=(-1)**(is-1)*t(2-is)
      if((is+1).ge.0) dum7=t(is+2)
      if((is+1).lt.0) dum7=(-1)**(is+1)*t(-is)
      dum2=0.5*(dum4-dum6)
      dum12=0.5*(dum5-dum7)
      dum23=dum2*dum3
      dum45=dum4*dum5
      dum67=dum6*dum7
      dum113=dum11*dum3
      dum112=dum11*dum12
      dum16=w-kz*betaz
      dum17=w+kz*betaz
      crcp=-ci*cdum1/(gamma*eps)*((gamma/omegae)*dum16*dum23
     &-0.5e0*xmnrl*(dum45-dum67))*(cdum8-1.0)
      crcm=-ci*cdum1/(gamma*epsp)*((gamma/omegae)*dum17*dum23
     &-0.5e0*xmnrl*(dum45-dum67))*(cdum10-1.0)
      cpsip=-cdum1/(gamma*rc*eps)*((gamma/omegae)*dum16*dfloat(ism)
     &*dum113/xmnrc-0.5e0*xmnrl*(dum45+dum67))*(cdum8-1.0)
      cpsim=-cdum1/(gamma*rc*epsp)*((gamma/omegae)*dum17*dfloat(ism)
     &*dum113/xmnrc-0.5e0*xmnrl*(dum45+dum67))*(cdum10-1.0)
      cppp=ci*cdum1/eps*dum16*dum112*(cdum8-1.0)
      cppm=ci*cdum1/epsp*dum17*dum112*(cdum10-1.0)
      cphip=cdum1*kz*omegae*betap/(gamma**2*betaz*eps**2)
     &*dum112*(cdum8-1.0+ci*tau*eps)
     &-cdum1/(gamma*betap*eps)*(dfloat(is)*dum16/xmnrl
     &-xmn*betap/rw)*dum113*(cdum8-1.0)
      cphim=-cdum1*kz*omegae*betap/(gamma**2*betaz*epsp**2)
     &*dum112*(cdum10-1.0+ci*tau*epsp)
     &-cdum1/(gamma*betap*epsp)*(dfloat(is)*dum17/xmnrl
     &-xmn*betap/rw)*dum113*(cdum10-1.0)
      cpzp=ci*cdum1*kz*betap/eps*dum112*(cdum8-1.0)
      cpzm=-ci*cdum1*kz*betap/epsp*dum112*(cdum10-1.0)
      ctimep=-cdum1*betap/(gamma*betaz*eps**2)*(w*betaz-kz)
     &*dum112*(cdum8-1.0+ci*tau*eps)
      ctimem=-cdum1*betap/(gamma*betaz*epsp**2)*(w*betaz+kz)
     &*dum112*(cdum10-1.0+ci*tau*epsp)
      y1(it,ipz,irc,ipsi,iphi)=y1(it,ipz,irc,ipsi,iphi)+(crcp-crcm)
      y2(it,ipz,irc,ipsi,iphi)=y2(it,ipz,irc,ipsi,iphi)+(cpsip-cpsim)
      y3(it,ipz,irc,ipsi,iphi)=y3(it,ipz,irc,ipsi,iphi)+(cppp-cppm)
      y4(it,ipz,irc,ipsi,iphi)=y4(it,ipz,irc,ipsi,iphi)+(cphip-cphim)
      y5(it,ipz,irc,ipsi,iphi)=y5(it,ipz,irc,ipsi,iphi)+(cpzp-cpzm)
      y6(it,ipz,irc,ipsi,iphi)=y6(it,ipz,irc,ipsi,iphi)+(ctimep-ctimem)
  100 continue
      return
      end
c *********************************************************************
      subroutine sety(w,weight,y1,y2,y3,y4,y5,y6)
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      dimension y1(1,201,21,1,1),y2(1,201,21,1,1),
     &y3(1,201,21,1,1),y4(1,201,21,1,1),y5(1,201,21,1,1),
     &y6(1,201,21,1,1),weight(1,201,21,1,1)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
c
c set initial electron coordinates and weighing factors
c
      suma=0.0
      tstep=2.0*pi/(w*dfloat(itmax))
      phistp=2.0e0*pi/dfloat(iphimx*is)
      if(im.ne.0) psistp=2.0e0*pi/dfloat(ipsimx*im)
      pzstep=(pzmax-pzmin)/dfloat(ipzmax)
      rcstep=(rc2-rc1)/dfloat(ircmax)
      do 120 it=1,itmax
      do 120 ipz=1,ipzmax
      pz=pzmin+pzstep*(dfloat(ipz)-0.5)
      pp=dsqrt(p0**2-pz**2)
      do 120 irc=1,ircmax
      rc=rc1+rcstep*(dfloat(irc)-0.5)
      if(rc.eq.0.0d0) rc=1.0d-7
      fac2=dexp(-(pz-pz0)**2/(2.0*delpz**2+1.0e-10))*rc
      do 120 ipsi=1,ipsimx
      do 120 iphi=1,iphimx
      y1(it,ipz,irc,ipsi,iphi)=rc
      if(ipsimx.eq.1) y2(it,ipz,irc,ipsi,iphi)=0.0e0
      if(ipsimx.ne.1) y2(it,ipz,irc,ipsi,iphi)=psistp*(dfloat(ipsi)-0.5)        
      y3(it,ipz,irc,ipsi,iphi)=pp
      y4(it,ipz,irc,ipsi,iphi)=phistp*(dfloat(iphi)-0.5e0)
     &+y2(it,ipz,irc,ipsi,iphi)
      y5(it,ipz,irc,ipsi,iphi)=pz
      y6(it,ipz,irc,ipsi,iphi)=tstep*dfloat(it-1)
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
      subroutine hotw01(w0hot,value)
c
c input parameters are fed through common, return w0hot and value.
c
      implicit real*8 (a,b,d-h,j-z), complex*16 (c)
      dimension cw(20)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      external cfunc1
      cw(1)=dcmplx(w0fst,0.0d0)
      ep1=1.0e-6
      ep2=ep1
      imaxit=50
      iroot=1
      icont1=0
      call muller(0,iroot,cw,imaxit,ep1,ep2,cfunc1,.true.)
      w0hot=cw(1)
      value=cdabs(cfunc1(cw(1)))
      return
      end
c *********************************************************************
      function cfunc1(cw)
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      w0hot=cw
      dum1=bj(im,xmn)**2*(1.0-dfloat(im*im)/xmn**2)
      call anavg0(w0hot,cavg1)
      dum2=-ci*cavg1
      dum3=2.2878e-7*pbkw*dum2/(w0fst**2*rw**2*(gamma0-1.0)*dum1)
      cfunc1=(w0hot/w0fst)**2-1.0+dum3
      icont1=icont1+1
      return
      end
c *********************************************************************
      subroutine derivy(yy,dyy,ieqfst,ieqlst)
c calculates derivatives of beam coordinates
      implicit  real*8  (a,b,d-h,j-z), complex*16 (c)
      dimension yy(20),dyy(20)
      dimension t(1000)
      common ci,pi,l,rw,q,kz,xmn,w0fst,w0,omegae,
     &gamma0,p0,pz0,pzmin,pzmax,pbkw,delpz,rcavg,rc1,rc2,is,im,
     &itmax,ipzmax,ircmax,ipsimx,iphimx,icont1,icont2,icont3
      common/freq/w
      common/field/cbz0
      icont3=icont3+1
      rc=yy(1)
      psi=yy(2)
      pp=yy(3)
      phi=yy(4)
      pz=yy(5)
      time=yy(6)
      z=yy(7)
      gamma=dsqrt(1.0e0+pp*pp+pz*pz)
      betap=pp/gamma
      betaz=pz/gamma
      omegac=omegae/gamma
      rl=pp/omegae
      xmnrc=xmn*rc/rw
      xmnrl=xmn*rl/rw
      cdum1=cbz0*rw/(2.0*xmn*pz)
      ism=is-im
      bs0=w*time-dfloat(is)*phi+dfloat(ism)*psi
     &+(dfloat(im)-dfloat(is)/2.0e0)*pi
      as=bs0-kz*z
      asp=bs0+kz*z
      dumc=dcos(as)
      dums=dsin(as)
      dumcp=dcos(asp)
      dumsp=dsin(asp)
      cdum18=dcmplx(dumc,-dums)*cdum1
      cdum20=dcmplx(dumcp,-dumsp)*cdum1
      cdum22=dcmplx(dumc-dumcp,dums-dumsp)
      call bes(iabs(ism)+2,xmnrc,0,result,t)
      if(ism.ge.0) dum11=t(ism+1)
      if(ism.lt.0) dum11=(-1)**ism*t(-ism+1)
      if((ism-1).ge.0) dum4=t(ism)
      if((ism-1).lt.0) dum4=(-1)**(ism-1)*t(2-ism)
      if((ism+1).ge.0) dum6=t(ism+2)
      if((ism+1).lt.0) dum6=(-1)**(ism+1)*t(-ism)
      call bes(iabs(is)+2,xmnrl,0,result,t)
      if(is.ge.0) dum3=t(is+1)
      if(is.lt.0) dum3=(-1)**is*t(-is+1)
      if((is-1).ge.0) dum5=t(is)
      if((is-1).lt.0) dum5=(-1)**(is-1)*t(2-is)
      if((is+1).ge.0) dum7=t(is+2)
      if((is+1).lt.0) dum7=(-1)**(is+1)*t(-is)
      dum2=0.5*(dum4-dum6)
      dum12=0.5*(dum5-dum7)
      dum23=dum2*dum3
      dum45=dum4*dum5
      dum67=dum6*dum7
      dum113=dum11*dum3
      dum112=dum11*dum12
      dum16=w-kz*betaz
      dum17=w+kz*betaz
      dcp=-(dum16*dum23/omegac-0.5*xmnrl*(dum45-dum67))
      dcm=-(dum17*dum23/omegac-0.5*xmnrl*(dum45-dum67))
      dpsip=(dum16*dfloat(ism)*dum113/(omegac*xmnrc)
     &      -0.5*xmnrl*(dum45+dum67))/rc
      dpsim=(dum17*dfloat(ism)*dum113/(omegac*xmnrc)
     &      -0.5*xmnrl*(dum45+dum67))/rc
      dpp=gamma*dum16*dum112
      dpm=gamma*dum17*dum112
      dphip=(dfloat(is)*dum16/xmnrl-xmn*betap/rw)*dum113/betap
      dphim=(dfloat(is)*dum17/xmnrl-xmn*betap/rw)*dum113/betap
      dzpm=kz*pp*dum112
      dyy(1)=dcp*cdum18-dcm*cdum20
      dyy(2)=ci*(dpsip*cdum18-dpsim*cdum20)
      dyy(3)=dpp*cdum18-dpm*cdum20
      dyy(4)=omegae/pz+ci*(dphip*cdum18-dphim*cdum20)
      dyy(5)=dzpm*(cdum18+cdum20)
      dyy(6)=gamma/pz
      dyy(7)=1.0
      if(ieqlst.le.7) return
      cdx=-0.5*w*rw*pp*dum112*cdum22/(xmn*pz)
      dyy(8)=cdx
      dyy(9)=-ci*cdx
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
c **********************************************************************
      function bi(in,x)
c
c bi: bessel function i(x) of integer order.
c in: positive or negative order of bi (-40< in <40 ).
c x : positive or negative argument of bi (x=0 or 1.0e-10<x <300 ).
c
      implicit  real*8 (a-h,o-z)
      dimension t(1000)
      isign=1
c     sx=sngl(x)
c     modified by c.s.hsue
      xabs=dabs(x)
      inabs=iabs(in)
      call bes(inabs,xabs,1,sbi,t)
      if(x.lt.0.0d0) isign=(-1)**inabs
      bi=sbi*dfloat(isign)
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
c **********************************************************************
      function step(x)
      implicit  real*8 (a-h,o-z)
      if(x.ge.0.0d0) go to 1
      step=0.0d0
      return
    1 step=1.0d0
      return
      end

