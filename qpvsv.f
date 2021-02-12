      program qpvsv
c this program calculates q times the start oscillation beam power (qp) vs
c the beam voltage for a gyromonotron by assuming a rf field given by that
c of an ideal cavity. 
c reference 1: k. r. chu, phys. fluids 21, 2354 (1978).
c reference 2: k. r. chu and d. dialetis, infrared and millimeter
c              wave book series (edited by k. j. button),
c              vol. 13, 45 (1985)
c note: unless noted otherwise, all quantities are normalized to the
c       radius of the cavity according to the procedure of ref.1 above.
c
c date of this version: March 1, 1999   
c
      implicit real*8 (a-h,j-z)
c
c     dimension axmn(51,40)
      dimension axmn(9,8)
      dimension aqp(501,51),avb(501),aqib(501),aalpha(501),aqibvb(501)
c
      data pi,e,e2,rc,rc2,me/3.1415926d0,4.8032d-10,2.30707d-19,
     &2.99792d10,8.98752d20,9.1095d-28/
c
c specify alpha as a function of vbkv
c
      falpha(vbkv)=alphaf*(vbkv/vbkvf)**nalpha
c
c specify q times ib as a function of vbkv
c
      fqib(vbkv)=qibf*(vbkv/vbkvf)**nqib
c
c all input and output quantities are dimensionless unless
c noted otherwise
c
c read in roots of bessel functions
c     read(45,95) axmn
   95 format(10f8.3)
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
c specify operating mode, harmonic, and frequency in hertz.
c (ismain and is > or = 1, inmain and in > or = 1,
c ilmain and il > or = 1, immain and im >, <, or = 0)
c
      ismain=1
      immain=-6
      inmain=2
      ilmain=1
c
      fhz=35.19e9
c
c specify final beam voltage (in kv) and final alpha
c
      vbkvf=600.0d0
      gammaf=1.0d0+vbkvf/511.0e0
      alphaf=1.0d0
c
c specify q times final beam current (in amp)
c
      qibf=250.*1200.
c
c specify beam guiding center position
c
      ipeak=1
      xmain=axmn(iabs(immain)+1,inmain)
      rho0=axmn(iabs(ismain-immain)+1,ipeak)/xmain
      if(ismain.eq.immain.and.ipeak.eq.1) rho0=0.0d0
      if(ismain.eq.immain.and.ipeak.ne.1) rho0=axmn(1,ipeak-1)/xmain
      rho0=0.724d0
c
c specify cavity length
c
      do 200 ill=1,1
c     llamda=4.7+0.12d0*(ill-1)
c     l=pi*dsqrt(4.0d0*llamda**2-dfloat(ilmain)**2)/xmain
      l=2.5d0
      llamda=0.5d0*dsqrt(xmain**2*l**2/(pi**2)+dfloat(ilmain)**2)
c
c specify alpha and q*ib profile parameters
c
      do 200 iprof=1,1
      nalpha=1.0d0+1.0d0*(iprof-1)
      nqib=1.0d0+1.0d0*(iprof-1)
c
c
c specify magnetic field
c
      do 200 idf=1,1
      df=0.5d0+(idf-1)*3.5d0+1.0d-4
c
c calculated parameters
c
      betapf=dsqrt((1.0d0-1.0d0/gammaf**2)/(1.0d0+1.0d0/alphaf**2))
      betazf=betapf/alphaf
      ppf=gammaf*betapf
      pzf=gammaf*betazf
      tauf=l/betazf
      kz=ilmain*pi/l
      w=dsqrt(kz**2+xmain**2)
      rwcm=w*rc/(2.0d0*pi*fhz)
      print 4, llamda,l,fhz,rwcm
    4 format("   l/lamda=",1pd11.4,", l=",1pd11.4,", f=",1pd11.4,
     &" hz, rw=",1pd11.4," cm"/)
      omegae=(w-kz*betazf-df/tauf)*gammaf/dfloat(ismain)
      omegac=omegae/gammaf
      df=(w-kz*betazf-dfloat(ismain)*omegac)*tauf
      ocwcf=omegac/xmain
      bgauss=omegae*1.707d3/rwcm
      ncycf=omegac*tauf/(2.0d0*pi)
      rholf=betapf/omegac
      dpf=df+2.0d0*kz*betazf*tauf
c
      print 22,ismain,immain,inmain,ilmain,xmain,vbkvf,
     &gammaf,alphaf,betapf,betazf,ppf,pzf,
     &rho0,tauf,ncycf,omegae,bgauss,df,dpf,ocwcf
   22 format(" ismain=",i3,", immain=",i3,", inmain=",i3,", ilmain=",i3,
     &", xmn=",1pd11.4/" vbf=",1pd10.3,
     &"kv, gammaf=",1pd10.3,", alphaf=",1pd10.3/" betapf=",1pd10.3,
     &", betazf=",1pd10.3,", ppf=",1pd10.3,", pzf=",1pd10.3/
     &" rho0=",1pd10.3,", tauf=",1pd9.2,", ncycf=",1pd9.2/
     &" omegae=",1pd10.3,"(",1pd10.3," gauss)"/" df=",1pd10.3,", dpf=",
     &1pd10.3,", omegac/xmain=",1pd10.3/)
      print 21,nalpha,nqib
   21 format(" alpha profile parameter:nalpha=",1pd10.3/
     &" q*ib profile parameter: nqib=",1pd10.3)
c
      do 30 ivb=1,501
      avb(ivb)=0.0d0
      aqib(ivb)=0.0d0
      aalpha(ivb)=0.0d0
      aqibvb(ivb)=0.0d0
      do 30 ip=1,51
   30 aqp(ivb,ip)=0.0d0
c
c specify range of beam voltage to be scanned
c
      vbkvi=1.0d-3
c
      ivbmin=1
      ivbmax=241
      iprint=20
c
      qpmin=1.0d30
      vbstep=(vbkvf-vbkvi)/dfloat(ivbmax-1)
      do 100 ivb=ivbmin,ivbmax
      vbkv=vbkvi+vbstep*(ivb-1)
      alpha=falpha(vbkv)
      qib=fqib(vbkv)
      aqib(ivb)=qib
      aalpha(ivb)=alpha
      aqibvb(ivb)=qib*vbkv
      gamma=1.0d0+vbkv/511.0d0
      avb(ivb)=vbkv
      betap=dsqrt((1.0d0-1.0d0/gamma**2)/(1.0d0+1.0d0/alpha**2))
      betaz=betap/alpha
      pp=gamma*betap
      pz=gamma*betaz
      tau=l/betaz
      kz=ilmain*pi/l
      w=dsqrt(kz**2+xmain**2)
      omegac=omegae/gamma
      d=(w-kz*betaz-dfloat(ismain)*omegac)*tau
      ocwc=omegac/xmain
      ncyc=omegac*tau/(2.0d0*pi)
      rhol=betap/omegac
      dp=d+2.0d0*kz*betaz*tau
c
      if(mod(ivb-1,iprint).eq.0) print 1,ismain,immain,inmain,ilmain,
     &xmain,vbkv,
     &gamma,alpha,betap,betaz,pp,pz,
     &rho0,l,tau,ncyc,omegae,bgauss,d,dp,ocwc
    1 format(/" ***************************************************"/
     &" ismain=",i3,", immain=",i3,", inmain=",i3,", ilmain=",i3,
     &", xmn=",1pd11.4/" vb=",1pd10.3,
     &"kv, gamma=",1pd10.3,", alpha=",1pd10.3/" betap=",1pd10.3,
     &", betaz=",1pd10.3,", pp=",1pd10.3,", pz=",1pd10.3/
     &" rho0=",1pd10.3,", l=",1pd9.2,", tau=",1pd9.2,", ncyc=",1pd9.2/
     &" omegae=",1pd10.3,"(",1pd10.3," gauss)"/" d=",1pd10.3,", dp=",
     &1pd10.3,", omegac/xmain=",1pd10.3/
     &" ***************************************************")
c
c specify modes and harmonics to be scanned
c
      ipar=1
      do 400 is=1,1
      do 400 imm=1,61
      im=imm-31
      do 400 in=1,15
      do 400 il=1,1
c
c calculated parameters
c
      kz=il*pi/l
      xmn=axmn(iabs(im)+1,in)
      w=dsqrt(kz**2+xmn**2)
      d=(w-kz*betaz-is*omegac)*tau
      dp=d+2.0d0*kz*betaz*tau
      if(d.gt.(-pi).and.d.lt.(2.0d0*pi)) go to 20
      if(dp.gt.(-pi).and.dp.lt.(2.0d0*pi)) go to 20
      go to 400
   20 continue
c
c calculating gain functions
c
      call couple(is,im,xmn*rho0,xmn*rhol,hsm,tsm,usm)
      fac2=hsm*betap**2*tau
      fac3=w**2-kz**2
      sind2=dsin(d/2.0d0)**2
      sind=dsin(d)
      sindp2=dsin(dp/2.0d0)**2
      sindp=dsin(dp)
      qqq=xmn*betap*usm
      f1a=fac2*(-fac3*(4.0d0*sind2-d*sind)
     &-kz*d*(d*sind-2.0d0*sind2)/l)/d**3
      f1b  =2.0d0*(tsm*(w-kz*betaz)-qqq)*sind2/d**2
      alpha1=f1a+f1b
      f2a=fac2*(fac3*(2.0d0*(1.0d0+d/dp)*sindp2-d*sindp)
     &+kz*d*(d*sindp-2.0d0*d*sindp2/dp)/l)/(dp*d**2)
      f2b  =-2.0d0*(tsm*(w-kz*betaz)-qqq)*sindp2/(d*dp)
      alpha2=f2a+f2b
      f3a=fac2*(-fac3*(4.0d0*sindp2-dp*sindp)
     &+kz*dp*(dp*sindp-2.0d0*sindp2)/l)/dp**3
      f3b  =2.0d0*(tsm*(w+kz*betaz)-qqq)*sindp2/dp**2
      alpha3=f3a+f3b
      f4a=fac2*(fac3*(2.0d0*(1.0d0+dp/d)*sind2-dp*sind)
     &-kz*dp*(dp*sind-2.0d0*dp*sind2/d)/l)/(d*dp**2)
      f4b  =-2.0d0*(tsm*(w+kz*betaz)-qqq)*sind2/(d*dp)
      alpha4=f4a+f4b
      fac4=bj(im,xmn)**2*(1.0d0-dfloat(im*im)/xmn**2)
      fac1= 2.0d0*tau**2/(gamma*fac4*w)
      f1=fac1*alpha1
      f2=fac1*alpha2
      f3=fac1*alpha3
      f4=fac1*alpha4
      f=f1+f2+f3+f4
c note: f here corresponds to f/nu in the manuscript
      fac5=me**2*rc**5/(e2*1.0d10)
      fac6a=f1a+f2a+f3a+f4a
      fac6b=f1b+f2b+f3b+f4b
      fac6=fac6a+fac6b
      if(fac6.eq.0.0d0) qp=1.0d50
      if(fac6.ne.0) qp=-w**2*betaz**2*gamma*(gamma-1.0d0)*fac4*fac5
     &/(2.0d0*l*fac6)
      if(fac6a.ne.0.0d0) rtoab=fac6b/fac6a
      if(fac6a.eq.0.0d0) rtoab=1.0d50
      if(qp.lt.0.0d0) go to 400
      if(is.ne.ismain.or.in.ne.inmain) go to 23
      if(im.ne.immain.or.il.ne.ilmain) go to 23
      aqp(ivb,51)=qp
      if(qp.lt.qpmin) qpmin=qp
      go to 24
   23 continue
      if(ipar.le.50) aqp(ivb,ipar)=qp
      ipar=ipar+1
   24 continue
c
c specify maximum qp (in kw) to be printed
c
      qpcut=3.0d8
c
      if(qp.gt.qpcut) go to 400
      if(mod(ivb-1,iprint).eq.0) print 3, is,im,in,il,xmn,
     &w,kz,rho0,rhol,hsm,tsm,usm
    3 format(/" is=",i3,", im=",i3,", in=",i3,", il=",i3,
     &", xmn=",1pd12.5,
     &", w=",1pd10.3,", kz=",1pd10.3/" rho0=",1pd10.3,
     &", rhol=",1pd10.3,", hsm=",1pd10.3,",tsm=",1pd10.3,
     &",usm=",1pd10.3)
      if(mod(ivb-1,iprint).eq.0) print 2, d,dp,f,f1,f2,f3,f4,rtoab,qp
    2 format( " d=",1pd10.3,", dp=",1pd9.2,
     &",f/nu=",1pd9.2,"(",4(1pd8.1,","),")"/" fac6b/fac6a=",
     &1pd8.1,", qp=",1pd10.3," kw")
  400 continue
  100 continue
      call sscale(avb(ivbmin),avb(ivbmax),0.0d0,alphaf)
      call splot(avb,aalpha,ivbmax,"aa  alpha vs vb(in kv)$")
      call sscale(avb(ivbmin),avb(ivbmax),0.0d0,qibf)
      call splot(avb,aqib,ivbmax,"ii  q*ib(in amp) vs vb(in kv)$")
      qpmax=5.0d0*qpmin
      qpmax=2.25d8
      call sscale(avb(ivbmin),avb(ivbmax),0.0d0,qpmax)
      call splot(avb,aqibvb,ivbmax,"b   $")
      do 450 ipar=1,50
  450 call splot(avb,aqp(ivbmin,ipar),ivbmax,"*   $")
      call splot(avb,aqp(ivbmin,51),ivbmax,"00  qp (kw) vs vb(in kv)$")
      print 7
    7 format(///)
  200 continue
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
c *********************************************************************
      function step(x)
      implicit  real*8 (a-h,o-z)
      if(x.ge.0.0d0) go to 1
      step=0.0d0
      return
    1 step=1.0d0
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

