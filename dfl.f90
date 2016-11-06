!Module version 1.0.0
!Initial release
!Copyleft by Dmitriy Samosvat

module dfl
  doubleprecision,save :: pi=3.1415926535897932384626433832795028841971693993751d0
  doubleprecision,save :: e=4.803242d-10
  doubleprecision,save :: m0=9.10938188d-28
  doubleprecision,save :: hbar=1.05457159642d-27
  doubleprecision,save :: kB=1.3806503d-16
  doubleprecision,save :: evolt=1.602176462d-12
  doubleprecision,save :: angstrem=1d-8
  doubleprecision,save :: c_light=2.99792458d10
  doubleprecision,save :: z_=0d0
  doublecomplex,save   :: im=(0d0,1d0)
  doubleprecision,allocatable :: ksi(:)
  doubleprecision,allocatable :: an(:)
  doubleprecision,allocatable :: factorial(:)
  doubleprecision,allocatable :: gamma2n(:)
  doubleprecision,save :: x_sav
  doubleprecision,save :: nu_sav
  doubleprecision,save :: a_sav,b_sav,c_sav
  doublecomplex,save :: z_s,z_sav
  integer,save :: n_sav,l_s,m_s
CONTAINS

doubleprecision function dint(f,a,b,m)
  doubleprecision dx,a,b,sum,x1,x2
  integer m,i,n
  doubleprecision f
  external f
  dx=(b-a)/m
  i=0
  sum=0d0
  do while (i<m)
    x1=a+i*dx
    x2=x1+dx
    sum=sum+gint(f,x1,x2)
    i=i+1
  end do
  dint=sum
end function dint

doublecomplex function cint(f,a,b,m)
  doubleprecision dx,a,b,x1,x2
  integer m,i,n
  doublecomplex sum
  doublecomplex f
  external f
  dx=(b-a)/m
  i=0
  sum=(0d0,0d0)
  do while (i<m)
    x1=a+i*dx
    x2=x1+dx
    sum=sum+gcint(f,x1,x2)
    i=i+1
  end do
!  print *,'ci',sum
  cint=sum
end function cint

doubleprecision function dint_double(f,a,b,c,d,m)
  doubleprecision dx,dy,a,b,c,d,sum,x1,x2,y1,y2
  integer m,i,j,n
  doubleprecision f
  external f
  dx=(b-a)/m
  dy=(d-c)/m
  i=0
  j=0
  sum=0d0
  do while (i<m)
    do while(j<m)
      x1=a+i*dx
      x2=x1+dx
      y1=c+j*dy
      y2=y1+dy
      sum=sum+gint_double(f,x1,x2,y1,y2)
      j=j+1
    end do
    j=0
    i=i+1
  end do
  dint_double=sum
end function dint_double

doubleprecision function dint_triple(f,a,b,c,d,p,q,m)
  doubleprecision dx,dy,dz,a,b,c,d,p,q,sum,x1,x2,y1,y2,z1,z2
  integer m,i,j,k,n
  doubleprecision f
  external f
  dx=(b-a)/m
  dy=(d-c)/m
  dz=(q-p)/m
  i=0
  j=0
  k=0
  sum=0d0
  do i=0,m-1
    do j=0,m-1
      do k=0,m-1
        x1=a+i*dx
        x2=x1+dx
        y1=c+j*dy
        y2=y1+dy
        z1=p+k*dz
        z2=z1+dz
        sum=sum+gint_triple(f,x1,x2,y1,y2,z1,z2)
      end do
    end do
  end do
  dint_triple=sum
end function dint_triple


doubleprecision function dint_six(f,a,b,c,d,p,q,a1,b1,c1,d1,p1,q1,m)
  doubleprecision dx,dy,dz,a,b,c,d,p,q,sum,x1,x2,y1,y2,z1,z2
  doubleprecision a1,b1,c1,d1,p1,q1
  doubleprecision u1,u2,v1,v2,w1,w2,du,dv,dw
  integer m,i,j,k,i1,j1,k1,n
  doubleprecision f
  external f
  dx=(b-a)/m
  dy=(d-c)/m
  dz=(q-p)/m
  du=(b1-a1)/m
  dv=(d1-c1)/m
  dw=(q1-p1)/m
  
  i=0
  j=0
  k=0
  i1=0
  j1=0
  k1=0
  sum=0d0
  do i=0,m-1
    do j=0,m-1
      do k=0,m-1
        do i1=0,m-1
          do j1=0,m-1
            do k1=0,m-1
              x1=a+i*dx
              x2=x1+dx
              y1=c+j*dy
              y2=y1+dy
              z1=p+k*dz
              z2=z1+dz
              u1=a1+i1*du
              u2=u1+du
              v1=c1+j1*dv
              v2=v1+dv
              w1=p1+k1*dw
              w2=w1+dw
              sum=sum+gint_six(f,x1,x2,y1,y2,z1,z2,u1,u2,v1,v2,w1,w2)
            end do
          end do
        end do
      end do
    end do
  end do
  dint_six=sum
end function dint_six

doubleprecision function cint_six(f,a,b,c,d,p,q,a1,b1,c1,d1,p1,q1,m)
  doubleprecision dx,dy,dz,a,b,c,d,p,q,x1,x2,y1,y2,z1,z2
  doubleprecision a1,b1,c1,d1,p1,q1
  doubleprecision u1,u2,v1,v2,w1,w2,du,dv,dw
  integer m,i,j,k,i1,j1,k1,n
  doublecomplex sum
  doublecomplex f
  external f
  dx=(b-a)/m
  dy=(d-c)/m
  dz=(q-p)/m
  du=(b1-a1)/m
  dv=(d1-c1)/m
  dw=(q1-p1)/m
  
  i=0
  j=0
  k=0
  i1=0
  j1=0
  k1=0
  sum=0d0
  do i=0,m-1
    do j=0,m-1
      do k=0,m-1
        do i1=0,m-1
          do j1=0,m-1
            do k1=0,m-1
              x1=a+i*dx
              x2=x1+dx
              y1=c+j*dy
              y2=y1+dy
              z1=p+k*dz
              z2=z1+dz
              u1=a1+i1*du
              u2=u1+du
              v1=c1+j1*dv
              v2=v1+dv
              w1=p1+k1*dw
              w2=w1+dw
              sum=sum+gcint_six(f,x1,x2,y1,y2,z1,z2,u1,u2,v1,v2,w1,w2)
            end do
          end do
        end do
      end do
    end do
  end do
  cint_six=sum
end function cint_six

doubleprecision function cint_triple(f,a,b,c,d,p,q,m)
  doubleprecision dx,dy,dz,a,b,c,d,p,q,x1,x2,y1,y2,z1,z2
  doublecomplex sum
  integer m,i,j,k,n
  doublecomplex f
  external f
  dx=(b-a)/m
  dy=(d-c)/m
  dz=(q-p)/m
  i=0
  j=0
  k=0
  sum=0d0
  do i=0,m-1
    do j=0,m-1
      do k=0,m-1
        x1=a+i*dx
        x2=x1+dx
        y1=c+j*dy
        y2=y1+dy
        z1=p+k*dz
        z2=z1+dz
        sum=sum+gcint_triple(f,x1,x2,y1,y2,z1,z2)
      end do
    end do
  end do
  cint_triple=sum
end function cint_triple

doublecomplex function cint_double(f,a,b,c,d,m)
  doubleprecision dx,dy,a,b,c,d,x1,x2,y1,y2
  integer m,i,j,n
  doublecomplex sum
  doublecomplex f
  external f
  dx=(b-a)/m
  dy=(d-c)/m
  i=0
  j=0
  sum=(0d0,0d0)
  do while (i<m)
    do while(j<m)
      x1=a+i*dx
      x2=x1+dx
      y1=c+j*dy
      y2=y1+dy
      sum=sum+gcint_double(f,x1,x2,y1,y2)
      j=j+1
    end do
    i=i+1
  end do
!  print *,'ci',sum
  cint_double=sum
end function cint_double

doubleprecision function mkint(f,a,b,m)
  doubleprecision a,b,omega,summa
  doubleprecision x1,ran
  integer i,m
  doubleprecision f
  external f
  omega=(b-a)
  summa=0d0
  do i=1,m
    call random_number(ran)
    x1=a+(b-a)*ran
    summa=summa+f(x1)
  end do
  mkint=summa*omega/m
end function mkint

doubleprecision function mkint_six(f,a,b,c,d,p,q,a1,b1,c1,d1,p1,q1,m)
  doubleprecision a,b,c,d,p,q,a1,b1,c1,d1,p1,q1,omega,summa
  doubleprecision x1,x2,x3,x4,x5,x6,ran
  integer i,m
  doubleprecision f
  external f
  omega=(b-a)*(d-c)*(q-p)*(b1-a1)*(d1-c1)*(q1-p1)
  summa=0d0
  do i=1,m
    call random_number(ran)
    x1=a+(b-a)*ran
    call random_number(ran)
    x2=c+(d-c)*ran
    call random_number(ran)
    x3=p+(q-p)*ran
    call random_number(ran)
    x4=a1+(b1-a1)*ran
    call random_number(ran)
    x5=c1+(d1-c1)*ran
    call random_number(ran)
    x6=p1+(q1-p1)*ran
    summa=summa+f(x1,x2,x3,x4,x5,x6)
  end do
  mkint_six=summa*omega/m  
end function mkint_six


doubleprecision recursive function Gamma(x) result(rv1)
  doubleprecision x,x1,xm,rv
  integer n
  x_sav=x
  rv=0d0
  n=200
  if(x<=0d0) then
    rv=Gamma(x+1d0)/x
  else if(x<=7d0) then
    rv=Gamma(x+1d0)/x
  else
    x1=max(1d0,abs(x*log(x)))
    x1=max(x1,abs(x))
    xm=500d0
    rv=rv+dint(toint1,0d0,x1,n)
    rv=rv+dint(toint1,x1,xm,n)
  end if
  rv1=rv
end function Gamma

doublecomplex function toint_l(t)
  doublecomplex ret
  doubleprecision t
  ret=((im)**(m_s))*(z_s+((z_s**2-1d0)**(0.5d0))*cos(t))**l_s*cos(m_s*t)
!  print *,'ret=',ret,'z_s=',z,'t=',t
  toint_l=ret
end function toint_l

doubleprecision function LegendrePlm(l,m,x)
  integer m,l
  doubleprecision x
  doublecomplex rv
  l_s=l
  m_s=m
  z_s=x
  rv=(Gamma(l+m+1d0)/(pi*Gamma(l+1d0)))*cint(toint_l,0d0,pi,100)
  LegendrePlm=dble(rv)
end function LegendrePlm

doublecomplex function SphericalHarmonicY(l,m,theta,phi)
  doublecomplex rv
  integer l,m
  doubleprecision theta,phi
  rv=exp(im*m*phi)*sqrt((2d0*l+1d0)*Fact(dble(l-m))/(4d0*pi*Fact(dble(l+m))))*LegendrePlm(l,m,cos(theta))
  SphericalHarmonicY=rv
end function SphericalHarmonicY

doubleprecision function tobes(t)
  doubleprecision t
  tobes=sin(t)**(2d0*nu_sav)*cos(x_sav*cos(t))
end function tobes

doublecomplex function ctobes(t)
  doubleprecision t
  ctobes=sin(t)**(2d0*nu_sav)*cos(z_sav*cos(t))
end function ctobes

doubleprecision function toy1(t)
  doubleprecision t
  toy1=sin(x_sav*sin(t))*(cos(t)**(2d0*nu_sav))
end function toy1

doubleprecision function toy2(t)
  doubleprecision t,rv
  rv=exp(-x_sav*sinh(t))*(cosh(t)**(2d0*nu_sav))
  if(isnan(rv))then
    rv=0d0
  end if
  toy2=rv
end function toy2


doublecomplex function ctoy1(t)
  doubleprecision t
  ctoy1=sin(z_sav*sin(t))*(cos(t)**(2d0*nu_sav))
end function ctoy1

doublecomplex function ctoy2(t)
  doubleprecision t
  doublecomplex rv
  rv=exp(-z_sav*sinh(t))*(cosh(t)**(2d0*nu_sav))
  if(isnan(Real(rv)).or. isnan(AImag(rv))) then
    rv=(0d0,0d0)
  end if
  ctoy2=rv
end function ctoy2

doubleprecision recursive function BesselY(nu,x) result(rv1)
!doubleprecision function BesselY(nu,x)
  doubleprecision nu,x,res
  integer n
  nu_sav=nu
  x_sav=x
  n=200
  if(nu<0d0) then
    rv1=2d0*(nu+1d0)*BesselY(nu+1d0,x)/x-BesselY(nu+2d0,x)
  else
    res=dint(toy1,0d0,pi/2,n)-(dint(toy2,0d0,1d0,n)+dint(toy2,1d0,10d0,n)+dint(toy2,10d0,100d0,n))
    !print *,'res=',dint(toy2,10d0,100d0,n)
    res=2d0*res*((x/2d0)**nu)/(Gamma(nu+1d0/2d0)*Gamma(1d0/2d0))
    rv1=res
  end if
end function BesselY


doublecomplex recursive function CBesselY(nu,z) result(rv1)
!doubleprecision function BesselY(nu,x)
  doubleprecision nu
  doublecomplex z,res
  integer n
  nu_sav=nu
  z_sav=z
  n=200
  if(nu<0d0) then
    rv1=2d0*(nu+1d0)*CBesselY(nu+1d0,z)/z-CBesselY(nu+2d0,z)
  else
    res=cint(ctoy1,0d0,pi/2,n)-(cint(ctoy2,0d0,1d0,n)+cint(ctoy2,1d0,10d0,n)+cint(ctoy2,10d0,100d0,n))
    !print *,'res=',dint(toy2,10d0,100d0,n)
    res=2d0*res*((z/2d0)**nu)/(Gamma(nu+1d0/2d0)*Gamma(1d0/2d0))
    rv1=res
  end if
end function CBesselY

doublecomplex function HankelH1(nu,x)
  doubleprecision nu,x
  HankelH1=BesselJ(nu,x)+im*BesselY(nu,x)
end function HankelH1

doublecomplex function HankelH2(nu,x)
  doubleprecision nu,x
  HankelH2=BesselJ(nu,x)-im*BesselY(nu,x)
end function HankelH2

doublecomplex function CHankelH1(nu,x)
  doubleprecision nu
  doublecomplex x
  CHankelH1=CBesselJ(nu,x)+im*CBesselY(nu,x)
end function CHankelH1

doublecomplex function CHankelH2(nu,x)
  doubleprecision nu
  doublecomplex x
  CHankelH2=CBesselJ(nu,x)-im*CBesselY(nu,x)
end function CHankelH2

doubleprecision recursive function BesselJ(nu,x) result(rv1)
  doubleprecision x,nu,res
  integer n,err
  nu_sav=nu
  x_sav=x
  n=200
  if(nu<0d0) then
    rv1=2d0*(nu+1d0)*BesselJ(nu+1d0,x)/x-BesselJ(nu+2d0,x)
  else
!    res=dint(toint2,-1d0,1d0,n)*((x/2d0)**nu)
    res=dint(tobes,0d0,pi,n)*((x/2d0)**nu)
    res=res/(Gamma(0.5d0)*Gamma(nu+0.5d0))
    !Additional check
!    if(x==0d0) then
!      if(nu<0d0) then
!        res=1d0/z_
!      else
!        if(nu==0d0) then
!          res=1d0
!        else
!          res=0d0
!        end if
!      end if
!    end if
  rv1=res
  end if
end function BesselJ

doublecomplex recursive function CBesselJ(nu,z) result(rv1)
  doubleprecision nu
  doublecomplex res,z
  integer n,err
  nu_sav=nu
  z_sav=z
  n=200
  if(nu<0d0) then
    rv1=2d0*(nu+1d0)*CBesselJ(nu+1d0,z)/z-CBesselJ(nu+2d0,z)
  else
!    res=dint(toint2,-1d0,1d0,n)*((x/2d0)**nu)
    res=cint(ctobes,0d0,pi,n)*((z/2d0)**nu)
    res=res/(Gamma(0.5d0)*Gamma(nu+0.5d0))
    !Additional check
!    if(x==0d0) then
!      if(nu<0d0) then
!        res=1d0/z_
!      else
!        if(nu==0d0) then
!          res=1d0
!        else
!          res=0d0
!        end if
!      end if
!    end if
  rv1=res
  end if
end function CBesselJ

doubleprecision function SphericalBesselJ(nu,x)
  doubleprecision nu,x
  SphericalBesselJ=Sqrt(pi/(2d0*x))*BesselJ(nu+1d0/2d0,x)
end function SphericalBesselJ

doubleprecision function SphericalBesselY(nu,x)
  doubleprecision nu,x
  SphericalBesselY=Sqrt(pi/(2d0*x))*BesselY(nu+1d0/2d0,x)
end function SphericalBesselY

doublecomplex function SphericalHankelH1(nu,x)
  doubleprecision nu,x
  SphericalHankelH1=Sqrt(pi/(2d0*x))*HankelH1(nu+1d0/2d0,x)
end function SphericalHankelH1

doublecomplex function SphericalHankelH2(nu,x)
  doubleprecision nu,x
  SphericalHankelH2=Sqrt(pi/(2d0*x))*HankelH2(nu+1d0/2d0,x)
end function SphericalHankelH2

!CSpherical functions
doublecomplex function CSphericalBesselJ(nu,x)
  doubleprecision nu
  doublecomplex x
  CSphericalBesselJ=Sqrt(pi/(2d0*x))*CBesselJ(nu+1d0/2d0,x)
end function CSphericalBesselJ

doublecomplex function CSphericalBesselY(nu,x)
  doubleprecision nu
  doublecomplex x
  CSphericalBesselY=Sqrt(pi/(2d0*x))*CBesselY(nu+1d0/2d0,x)
end function CSphericalBesselY

doublecomplex function CSphericalHankelH1(nu,x)
  doubleprecision nu
  doublecomplex x
  CSphericalHankelH1=Sqrt(pi/(2d0*x))*CHankelH1(nu+1d0/2d0,x)
end function CSphericalHankelH1

doublecomplex function CSphericalHankelH2(nu,x)
  doubleprecision nu
  doublecomplex x
  CSphericalHankelH2=Sqrt(pi/(2d0*x))*CHankelH2(nu+1d0/2d0,x)
end function CSphericalHankelH2


doubleprecision function BesselK(nu,x)
  doubleprecision x,nu,res
  integer n,err
  nu_sav=nu
  x_sav=x
  n=200
!  res=dint(toint4,0d0,10d0,n)+dint(toint4,10d0,100d0,n)+dint(toint4,100d0,1000d0,n)
!  res=(1/2d0)*res*(x/2d0)**nu
  res=dint(tointk,0d0,10d0,n)+dint(tointk,10d0,100d0,n)+dint(tointk,100d0,1000d0,n)
!  res=(1/2d0)*res*(x/2d0)**nu
  BesselK=res
end function BesselK

doubleprecision function tointk(t)
  doubleprecision rv,t
  tointk=exp(-x_sav*cosh(t))*cosh(nu_sav*t)
!  print *,'DZEBUG=',rv,'t=',t,'nu_sav=',nu_sav,'x_sav=',x_sav
end function tointk


doubleprecision function toint1(t)
  doubleprecision t
  toint1=exp(-t)*(t**(x_sav-1d0))
end function toint1

doubleprecision function toint2(t)
  doubleprecision t
  toint2=((1d0-t**2)**(nu_sav-0.5d0))*cos(x_sav*t)
end function toint2

doublecomplex function toint3(t)
  doubleprecision t
  doublecomplex res
  res=t**(a_sav-1d0)*(1-t)**(c_sav-a_sav-1d0)*(1-x_sav*t)**(-b_sav)
!  print *,'r=',res
  toint3=res
end function toint3

doubleprecision function toint4(t)
  doubleprecision t
  toint4=exp(-t-x_sav**2/(4d0*t))*t**(-nu_sav-1d0)
end function toint4
!Produce gauss integrate
!INPUT n - order
!doubleprecision function gint(f,a,b,n)

doubleprecision function gint(f,a,b)
  integer i,n
  doubleprecision a,b,x,eta,rv
  doubleprecision f
  external f
  rv=0d0
  n=n_sav
  do i=1,n
    x=ksi(i)*(b-a)/2d0+(a+b)/2d0
    eta=f(x)*(b-a)/2d0
    rv=rv+an(i)*eta
  end do
  gint=rv
end function gint

doublecomplex function gcint(f,a,b)
  integer i,n
  doubleprecision a,b,x
  doublecomplex eta,rv
  doublecomplex f
  external f
  rv=(0d0,0d0)
  n=n_sav
  do i=1,n
    x=ksi(i)*(b-a)/2d0+(a+b)/2d0
    eta=f(x)*(b-a)/2d0
    rv=rv+an(i)*eta
  end do
!  print *,'gc=',rv
  gcint=rv
end function gcint

doubleprecision function gint_double(f,a,b,c,d)
  integer i,j,n
  doubleprecision a,b,c,d,x,y,eta,rv
  doubleprecision f
  external f
  rv=0d0
  n=n_sav
  do i=1,n
    do j=1,n
      x=ksi(i)*(b-a)/2d0+(a+b)/2d0
      y=ksi(j)*(d-c)/2d0+(c+d)/2d0
      eta=f(x,y)*(b-a)*(d-c)/4d0
!      zeta=f(x,y)*(d-c)/2d0
      rv=rv+an(i)*an(j)*eta
    end do
  end do
  gint_double=rv
end function gint_double

doubleprecision function gint_triple(f,a,b,c,d,p,q)
  integer i,j,k,n
  doubleprecision a,b,c,d,p,q,x,y,z,eta,rv
  doubleprecision f
  external f
  rv=0d0
  n=n_sav
  do i=1,n
    do j=1,n
      do k=1,n
        x=ksi(i)*(b-a)/2d0+(a+b)/2d0
        y=ksi(j)*(d-c)/2d0+(c+d)/2d0
        z=ksi(k)*(q-p)/2d0+(q+p)/2d0
        eta=f(x,y,z)*(b-a)*(d-c)*(q-p)/8d0
!      zeta=f(x,y)*(d-c)/2d0
        rv=rv+an(i)*an(j)*an(k)*eta
      end do
    end do
  end do
  gint_triple=rv
end function gint_triple

doubleprecision function gint_six(f,a,b,c,d,p,q,a1,b1,c1,d1,p1,q1)
  integer i,j,k,i1,j1,k1,n
  doubleprecision a,b,c,d,p,q,a1,b1,c1,d1,p1,q1,x,y,z,u,v,w,eta,rv
  doubleprecision f
  external f
  rv=0d0
  n=n_sav
  do i=1,n
    do j=1,n
      do k=1,n
        do i1=1,n
          do j1=1,n
            do k1=1,n
              x=ksi(i)*(b-a)/2d0+(a+b)/2d0
              y=ksi(j)*(d-c)/2d0+(c+d)/2d0
              z=ksi(k)*(q-p)/2d0+(q+p)/2d0
              u=ksi(i1)*(b1-a1)/2d0+(a1+b1)/2d0
              v=ksi(j1)*(d1-c1)/2d0+(c1+d1)/2d0
              w=ksi(k1)*(q1-p1)/2d0+(q1+p1)/2d0
              eta=f(x,y,z,u,v,w)*(b-a)*(d-c)*(q-p)*(b1-a1)*(d1-c1)*(q1-p1)/64d0
              rv=rv+an(i)*an(j)*an(k)*an(i1)*an(j1)*an(k1)*eta
            end do
          end do
        end do
      end do
    end do
  end do
  gint_six=rv
end function gint_six

doublecomplex function gcint_six(f,a,b,c,d,p,q,a1,b1,c1,d1,p1,q1)
  integer i,j,k,i1,j1,k1,n
  doubleprecision a,b,c,d,p,q,a1,b1,c1,d1,p1,q1,x,y,z,u,v,w,eta
  doublecomplex rv
  doublecomplex f
  external f
  rv=0d0
  n=n_sav
  do i=1,n
    do j=1,n
      do k=1,n
        do i1=1,n
          do j1=1,n
            do k1=1,n
              x=ksi(i)*(b-a)/2d0+(a+b)/2d0
              y=ksi(j)*(d-c)/2d0+(c+d)/2d0
              z=ksi(k)*(q-p)/2d0+(q+p)/2d0
              u=ksi(i1)*(b1-a1)/2d0+(a1+b1)/2d0
              v=ksi(j1)*(d1-c1)/2d0+(c1+d1)/2d0
              w=ksi(k1)*(q1-p1)/2d0+(q1+p1)/2d0
              eta=f(x,y,z,u,v,w)*(b-a)*(d-c)*(q-p)*(b1-a1)*(d1-c1)*(q1-p1)/64d0
              rv=rv+an(i)*an(j)*an(k)*an(i1)*an(j1)*an(k1)*eta
            end do
          end do
        end do
      end do
    end do
  end do
  gcint_six=rv
end function gcint_six



doublecomplex function gcint_double(f,a,b,c,d)
  integer i,j,n
  doubleprecision a,b,c,d,x,y
  doublecomplex eta,zeta,rv
  doublecomplex f
  external f
  rv=(0d0,0d0)
  n=n_sav
  do i=1,n
    do j=1,n
      x=ksi(i)*(b-a)/2d0+(a+b)/2d0
      y=ksi(j)*(d-c)/2d0+(c+d)/2d0
      eta=f(x,y)*(b-a)*(d-c)/4d0
!      zeta=f(x,y)*(d-c)/2d0
      rv=rv+an(i)*an(j)*eta
    end do
  end do
!  print *,'gc=',rv
  gcint_double=rv
end function gcint_double

doublecomplex function gcint_triple(f,a,b,c,d,p,q)
  integer i,j,k,n
  doubleprecision a,b,c,d,p,q,x,y,z
  doublecomplex eta,zeta,rv
  doublecomplex f
  external f
  rv=(0d0,0d0)
  n=n_sav
  do i=1,n
    do j=1,n
      do k=1,n
        x=ksi(i)*(b-a)/2d0+(a+b)/2d0
        y=ksi(j)*(d-c)/2d0+(c+d)/2d0
        z=ksi(k)*(q-p)/2d0+(q+p)/2d0
        eta=f(x,y)*(b-a)*(d-c)*(q-p)/8d0
!       zeta=f(x,y)*(d-c)/2d0
        rv=rv+an(i)*an(j)*an(k)*eta
      end do
    end do
  end do
!  print *,'gc=',rv
  gcint_triple=rv
end function gcint_triple


subroutine istart(n)
  doubleprecision sep,x0,x1
  integer i,n,nf
  allocate(ksi(n))
  allocate(an(n))
  i=1
  n_sav=n
  x0=-1d0
  x1=1d0
  sep=1d0/(n**2)
  do while (i<=n)
    ksi(i)=izerolim(tst,x0,x1,sep,i,nf)
    an(i)=2d0/((1d0-ksi(i)**2)*(DLegendreP(n,ksi(i))**2))
    i=i+1
  end do
end subroutine istart

subroutine istop()
  deallocate(ksi)
  deallocate(an)
end subroutine istop

doubleprecision function tst(x)
  doubleprecision x
  tst=LegendreP(n_sav,x)
end function tst

doubleprecision function LegendreP(n,x)
  integer n,i
  doubleprecision x,r,t
  doubleprecision arr(3)
  arr(1)=1d0
  arr(2)=x
  i=1
  do while(i<n)
    arr(3)=((2d0*i+1d0)*x*arr(2)-i*arr(1))/(i+1d0)
    arr(1)=arr(2)
    arr(2)=arr(3)
    i=i+1
  end do
  r=arr(2)
  if(n==0) then
    r=arr(1)
  end if
  LegendreP=r
end function LegendreP

doubleprecision function DLegendreP(n,x)
  doubleprecision x
  integer n
  DLegendreP=(LegendreP(n+1,x)-x*LegendreP(n,x))*(n+1d0)/(x**2-1d0)
end function DLegendreP

doubleprecision function Fact(x)
  doubleprecision x
  Fact=Gamma(x+1d0)
end function Fact

doubleprecision function ClebschGordan(j1,m1,j2,m2,j3,m3)
  doubleprecision j1,m1,j2,m2,j3,m3,rval,a,b,n,xmin,xmax
  integer nmin,nmax,i,v
  a=0d0
  b=0d0
  if((m1+m2/=m3).or.(abs(m1)>j1).or.(abs(m2)>j2).or.(abs(m3)>j3)) then
    rval=0d0
  else
    if((abs(j1-j2)<=j3).and.(abs(j1+j2)>=j3)) then
      xmin=0d0
      xmin=max(xmin,j2-j3-m1)
      xmin=max(xmin,j1-j3+m2)
      xmax=min(j2+m2,j1+j2-j3)
      xmax=min(xmax,j1-m1)
      nmin=idnint(xmin)
      nmax=idnint(xmax)
      a=sqrt((2d0*j3+1d0)*Fact(j1+j2-j3)*Fact(j3+j1-j2)*Fact(j3+j2-j1)/Fact(j1+j2+j3+1))
      b=0d0
      do i=nmin,nmax
        n=i
        b=b+((-1)**i)*sqrt(Fact(j1+m1)*Fact(j1-m1)*Fact(j2+m2)*Fact(j2-m2)*Fact(j3+m3)*Fact(j3-m3))/                     (Fact(n)*Fact(j1+j2-j3-i)*Fact(j1-m1-i)*Fact(j2+m2-i)*Fact(j3-j2+m1+i)*Fact(j3-j1-m2+i))
      end do
      rval=a*b
    else
      rval=0
    end if
  end if
  ClebschGordan=rval
end function ClebschGordan

doubleprecision function to_ai(t)
  doubleprecision t
  to_ai=(1d0/pi)*cos(t**3/3+x_sav*t)
end function to_ai

doubleprecision function AiryAi(x)
  doubleprecision rv,x
  doubleprecision a,b,mx,my
  integer n,i
  x_sav=x
  n=100000
  a=0d0
  b=1d0
  rv=0d0
  if(x<0d0) then
    mx=-x
    my=(2d0/3d0)*(mx**(3d0/2d0))
    rv=(sqrt(mx)/3d0)*(BesselJ(-1d0/3d0,my)+BesselJ(1d0/3d0,my))
  else
    if(x==0d0) then
      rv=1d0/(3d0**(2d0/3d0)*Gamma(2d0/3d0))
    else
      rv=(1d0/pi)*sqrt(x/3d0)*BesselK(1d0/3d0,(2d0/3d0)*(x**(3d0/2d0)))
    end if

  end if
!  rv=dint(to_ai,0d0,1d0,n)+dint(to_ai,1d0,10d0,n)!+dint(to_ai,10d0,100d0,n)+dint(to_ai,100d0,300d0,n)
  AiryAi=rv
end function AiryAi

subroutine gen_fact(n)
  integer n,i
  if(.not. allocated(factorial)) then
    allocate(factorial(n+1))
    do i=0,n
      factorial(i+1)=Fact(dble(i))
    end do
  end if
  if(.not. allocated(gamma2n)) then
    allocate(gamma2n(2*n+1))
    do i=1,2*n+1
      gamma2n(i)=Gamma(i/2d0)
    end do
  end if
end subroutine gen_fact

doubleprecision function Fact1(n)
  integer n
  Fact1=factorial(n+1)
end function Fact1

doubleprecision function ClebschGordan1(j1,m1,j2,m2,j3,m3)
integer j1,m1,j2,m2,j3,m3,nmin,nmax,i,v
doubleprecision rval,a,b,xmin,xmax
  a=0d0
  b=0d0
  if((m1+m2/=m3).or.(abs(m1)>j1).or.(abs(m2)>j2).or.(abs(m3)>j3)) then
    rval=0d0
  else
    if((abs(j1-j2)<=j3).and.(abs(j1+j2)>=j3)) then
      nmin=0
      nmin=max(nmin,j2-j3-m1)
      nmin=max(nmin,j1-j3+m2)
      nmax=min(j2+m2,j1+j2-j3)
      nmax=min(nmax,j1-m1)
 !     nmin=idnint(xmin)
 !     nmax=idnint(xmax)
      a=sqrt((2d0*j3+1d0)*Fact1(j1+j2-j3)*Fact1(j3+j1-j2)*Fact1(j3+j2-j1)/Fact1(j1+j2+j3+1))
      b=0d0
      do i=nmin,nmax
        b=b+((-1)**i)*sqrt(Fact1(j1+m1)*Fact1(j1-m1)*Fact1(j2+m2)*Fact1(j2-m2)*Fact1(j3+m3)*Fact1(j3-m3))/(Fact1(i)*Fact1(j1+j2-j3-i)*Fact1(j1-m1-i)*Fact1(j2+m2-i)*Fact1(j3-j2+m1+i)*Fact1(j3-j1-m2+i))
      end do
      rval=a*b
    else
      rval=0d0
    end if
  end if
  ClebschGordan1=rval
end function ClebschGordan1

doubleprecision function SixJSymbol(j1,j2,j12,j3,j,j23)
integer j1,j2,j12,j3,j,j23,m1,m2,m3,m12,m23,m,m0
doubleprecision result,s0,s1
result=0d0
if (abs(j1-j2)>j12 .or. ((j1+j2)<j12) .or. (abs(j3-j)>j12) .or. ((j3+j)<j12) .or. (abs(j3-j2)>j23) .or. ((j3+j2)<j23) .or. (abs(j1-j)>j23) .or. ((j1+j)<j23))then
  result=0d0
  SixJSymbol=result
else
s0=sqrt(2d0*j+1)*sqrt(2d0*j+1)*sqrt(2d0*j12+1)*sqrt(2d0*j23+1)*(-1)**(j1+j2+j3+j)
s1=0
do m1=-j1,j1
  do m2=-j2,j2
    do m3=-j3,j3
      m12=m1+m2
      m=m12+m3
      m23=m2+m3
      s1=s1+ClebschGordan1(j1,m1,j2,m2,j12,m12)*ClebschGordan1(j12,m12,j3,m3,j,m)*ClebschGordan1(j1,m1,j23,m23,j,m)*ClebschGordan1(j2,m2,j3,m3,j23,m23)
    end do
  end do
end do
SixJSymbol=s1/s0
end if
end function SixJSymbol

              
doublecomplex function Hypergeometric2F1(a,b,c,x)
  doubleprecision a,b,c,x
  doublecomplex res
  a_sav=a
  b_sav=b
  c_sav=c
  x_sav=x
  if(abs(x)<1d0) then
    res=cint(toint3,0d0,1d0,100)
    res=res*Gamma(c)/(Gamma(a)*Gamma(c-a))
  else
    res=0d0
  end if
  Hypergeometric2F1=res
end function Hypergeometric2F1

doubleprecision function ak(a,k)    
  doubleprecision a,rv
  integer k,i
  rv=1d0
  do i=0,(k-1)
    rv=rv*(a+i)
  end do
  ak=rv
end function ak

doubleprecision function AppelF4(a,b,c,d,x,y)
  doubleprecision a,b,c,d,x,y,rv,s,eps
  integer m,n
  logical cont
  cont=.true.
  m=0
  n=0
  rv=0d0
  s=0d0
  eps=1d-14
!  print *,'a=',a,'b=',b,'c=',c,'d=',d
  do m=0,20
    do n=0,20
      s=(ak(a,m+n)*ak(b,m+n)/(ak(c,m)*ak(d,n)*Fact1(m)*Fact1(n)))*(x**m)*(y**n)
!      s=(ak(a,m+n)*ak(b,m+n)/(ak(c,m)*ak(d,n)*Fact1(m)*Fact1(n)))*(xpow(m+1))*(ypow(n+1))
!      print *,'m=',m,'n=',n,'dzebug=',s
      rv=rv+s
    end do
  end do
!  print *,'Appel=',rv
  AppelF4=rv
end function AppelF4

integer function KroneckerDelta(n1,n2)
  integer n1,n2,rv
  if(n1==n2) then
    rv=1
  else
    rv=0
  end if
  KroneckerDelta=rv
end function KroneckerDelta

doubleprecision function to_ei(x)
  doubleprecision x
  to_ei=exp(x)/x
end function to_ei

doubleprecision function to_k(x)
  doubleprecision x
  to_k=1d0/(sqrt(1-x_sav*sin(x)**2))
end function to_k

doubleprecision function EllipticK(x)
  doubleprecision x,res
  x_sav=x
  res=dint(to_k,0d0,pi/2d0,100)
  EllipticK=res
end function EllipticK

doubleprecision function to_e(x)
  doubleprecision x
  to_e=sqrt(1-x_sav*sin(x)**2)
end function to_e

doubleprecision function EllipticE(x)
  doubleprecision x,res
  x_sav=x
  res=dint(to_e,0d0,pi/2d0,100)
  EllipticE=res
end function EllipticE

doubleprecision function Ei(x)
  doubleprecision rval,epsilon,x
  epsilon=1d-10
  if(x<0d0) then
    rval=dint(to_ei,-1000d0,x,10000)
  else
    if(x>0d0) then
      rval=dint(to_ei,-1000d0,-1d0,10000)+dint(to_ei,-1d0,-epsilon,10000)+dint(to_ei,epsilon,1d0,10000)+dint(to_ei,1d0,x,10000)
    else
      rval=z_/z_
    end if
  end if
  Ei=rval
end function Ei                      

integer function sign_old(x)
  doubleprecision x
  integer rv
  if(x>0d0) then
    rv=1
  end if
  if(x<0d0) then
    rv=-1
  end if
  if(x==0d0) then
    rv=0
  end if
  sign_old=rv
end function sign_old

doublecomplex function csjn(n,x)
  integer n
  doubleprecision nu
  doublecomplex x,result
  nu=dble(n)+0.5d0
  result=sqrt(pi/(2d0*x))*CBesselJ(nu,x)
  csjn=result
end function csjn

doublecomplex function csyn(n,x)
  integer n
  doubleprecision nu
  doublecomplex x,result
  nu=dble(n)+0.5d0
  result=sqrt(pi/(2d0*x))*(1d0/sin(nu*pi))*(-CBesselJ(-nu,x))
  csyn=result
end function csyn

doublecomplex function csh1(n,x)
  integer n
  doublecomplex x
  csh1=csjn(n,x)+im*csyn(n,x)
end function csh1

doublecomplex function csh2(n,x)
  integer n
  doublecomplex x
  csh2=csjn(n,x)-im*csyn(n,x)
end function csh2

doubleprecision function skm(x)
  doubleprecision x
  skm=(pi/2d0)*(exp(-x))/x
end function skm

doubleprecision function sk0(x)
  doubleprecision x
  sk0=(pi/2d0)*(exp(-x))/x
end function sk0

doubleprecision function sk1(x)
  doubleprecision x
  sk1=(pi/2d0)*(exp(-x)*(x+1d0))/(x**2)
end function sk1

doubleprecision function skn(j,x)
  integer j,i
  doubleprecision x,result
  doubleprecision,allocatable :: prev(:)
  allocate(prev(j+1))
  if(j==-1) then
    skn=skm(x)
  else
    do i=0,j
      if(i==0) then
        prev(i+1)=sk0(x)
      else if(i==1) then
        prev(i+1)=sk1(x)
      else
        prev(i+1)=((2d0*i-1d0)/x)*prev(i)+prev(i-1)
      end if
    end do
    result=prev(j+1)
    deallocate(prev)
    skn=result
  end if
end function skn


doubleprecision function sjm(x)
  doubleprecision x
  sjm=cos(x)/x
end function sjm

doubleprecision function sj0(x)
  doubleprecision x
  if(x==0d0) then
    sj0=1d0
  else
    sj0=sin(x)/x
  end if
end function sj0

doubleprecision function sj1(x)
  doubleprecision x
  if(x==0d0) then
    sj1=0d0
  else
    sj1=(sin(x)-x*cos(x))/(x**2)
  end if
end function sj1

doubleprecision function sjn(j,x)
  integer j,s,k,i
  doubleprecision x,epsilon,result,add
  doubleprecision,allocatable :: prev(:)
  allocate(prev(j+1))
  if(j==-1) then
    sjn=sjm(x)
  else if(x==0d0) then
    sjn=0d0
  else if(x*1.3d0>j) then
    do i=0,j
      if(i==0) then
        prev(i+1)=sj0(x)
      else if(i==1) then
        prev(i+1)=sj1(x)
      else
        prev(i+1)=((2d0*i-1d0)/x)*prev(i)-prev(i-1)
      end if
    end do
    result=prev(j+1)
    deallocate(prev)
    sjn=result
  else
    epsilon=1d-30
    result=0d0
    add=1
    k=0
    s=1!Sign
    do while ((abs(add)>epsilon).or.(k<50))
      add=s/x
      do i=0,j
        add=add*x/(2d0*i+1d0)
      end do
      do i=1,k
        add=add*(x**2)/(2d0*i*(2d0*(i+j)+1d0))
        end do
      s=-s
      k=k+1
      result=result+add
    end do
    sjn=result
  end if
end function sjn

doubleprecision function zero(f,a,b,errcode)
  doubleprecision z1,z2,z3,znew,zold
  doubleprecision f1,f2,f3
  doubleprecision a,b,epsilon,er_rel
  doubleprecision qj,Aj,Bj,Cj
  integer errcode,maxiter,cnt,cont_iter
  doubleprecision f
  external f
  maxiter=100
  cnt=0
  errcode=0
  epsilon=1d-18
  cont_iter=1
!Addition check for NaN
  if(f(a)*f(b)>=0d0 .or. f(a)*f(b)<0d0) then
    if(f(a)*f(b)>=0) then
      errcode=1
      if(f(a)==0d0) then
        errcode=0
        zero=a
      else if(f(b)==0d0) then
        errcode=0
        zero=b
      else
        errcode=1
        zero=(z_/z_) !Nothing found
      end if
    else
      z1=a
      z2=b
      z3=(a+b)/2d0
      zold=z1
      do while (cnt<maxiter .and. cont_iter==1)
        f1=f(z1)
        f2=f(z2)
        f3=f(z3)
        qj=(z3-z2)/(z2-z1)
        Aj=qj*f3-qj*(1d0+qj)*f2+f1*qj**2
        Bj=(2d0*qj+1d0)*f3-((1d0+qj)**2)*f2+f1*qj**2
        Cj=f3*(1d0+qj)
        znew=z3-sign_old(Bj)*(z3-z2)*(2d0*Cj)/(abs(Bj)+sqrt(Bj**2-4d0*Aj*Cj))
        if(.not. (znew>=0d0 .or. znew<0d0)) then
          cont_iter=0
        else
          if(f1*f3>0d0) then
            z1=z2
            z2=z3
            z3=znew
          else
            z2=z3
            z3=znew
          end if
        end if
        er_rel=abs((znew-zold)/znew)
        if(er_rel<epsilon) then
!          print *,'root found. cnt=',cnt,'eps=',er_rel
          cont_iter=0
        end if
        cnt=cnt+1
        zold=znew
      end do
!      print *,'Max iter counted. cnt=',cnt,'eps=',er_rel
      zero=z2
    end if
  else
    errcode=1
    zero=(z_/z_)
  end if
end function zero

logical function isnan(x)
  logical rv
  doubleprecision x
  rv=.false.
  if(x>0d0) then
    rv=.false.
  else
    if(x<=0d0) then
      rv=.false.
    else
      rv=.true.
    end if
  end if
  isnan=rv
end function isnan

doubleprecision function izerolim(f,x0,x1,sep,i,nf)
  doubleprecision x0,x1,sep,x,rv
  doubleprecision a,b,epsilon
  integer i,cnt,nf,err
  doubleprecision f
  external f
  !First check x0 point
!  print *,'i=',i
  cnt=0
  nf=0
  epsilon=1d-10
  a=x0+sep*epsilon
  b=x0+sep+sep*epsilon
  if(i==0) then
    print *,'WARNING!! YOU WANT TO CHECK x0 point. IT MAY CONTAIN SINGULARITY'
    nf=1
    if(f(x0)==0d0) then
      nf=0
    end if
    izerolim=x0
  else
    do while ((cnt<i) .and. (a<x1))
!Note that i cannot call function f in points greater than x1
      if(b>=x1) then
        b=x1-sep*epsilon!Note that point x1 may contain singularity
        if(f(a)*f(b)<0) then
          cnt=cnt+1
        end if
      else
        if(f(a)*f(b)<0) then
          cnt=cnt+1
        end if
        b=b+sep
      end if
      a=a+sep!Always increment a point
    end do
    if(cnt<i) then
      nf=1!Nothing found
      izerolim=x1
    else
      if(a<x1) then
        b=b-sep!Usually. Else b=x1-sep*epsilon
      end if
      a=a-sep
      izerolim=zero(f,a,b,err)
    end if
  end if
end function izerolim

subroutine mmul(N,H,A)
  integer N,i,j,k
  double precision H(N,N),A(N,N),TMP(N,N)
  do i=1,N
    do j=1,N
      TMP(i,j)=0
      do k=1,N
        TMP(i,j)=TMP(i,j)+H(i,k)*A(k,j)
      end do
    end do
  end do
  do i=1,N
    do j=1,N
      A(i,j)=TMP(i,j)
    end do
  end do
end subroutine mmul

subroutine mmul_left(N,H,A)
  integer N,i,j,k
  double precision H(N,N),A(N,N),TMP(N,N)
  do i=1,N
    do j=1,N
      TMP(i,j)=0
      do k=1,N
        TMP(i,j)=TMP(i,j)+H(i,k)*A(k,j)
      end do
    end do
  end do
  do i=1,N
    do j=1,N
      H(i,j)=TMP(i,j)
    end do
  end do
end subroutine mmul_left

subroutine mvmul(N,A,V)
  integer N,i,j,k
  double precision A(N,N),V(N),TMP(N)
  do i=1,N
    TMP(i)=0
    do j=1,N
      TMP(i)=TMP(i)+A(i,j)*V(j)
    end do
  end do
  do i=1,N
    V(i)=TMP(i)
  end do
end subroutine mvmul

subroutine trans(N,A)
  integer N,i,j
  double precision A(N,N),tmp
  do i=1,N
    do j=1,i
      tmp=A(i,j)
      A(i,j)=A(j,i)
      A(j,i)=tmp
    end do
  end do
end subroutine trans

subroutine QR(N,Asav,E,A)!E=R,A=Q
  integer N,i,j,k
  doubleprecision A(N,N),Asav(N,N),E(N,N),H(N,N),V(N)
  doubleprecision norm
  do i=1,N
    do j=1,N
      A(i,j)=Asav(i,j)
    end do
  end do
  do i=1,N
    do j=1,i
      E(i,j)=0
      E(j,i)=0
      E(i,i)=1
    end do
  end do
  do k=1,N-1
    norm=0
    do i=1,N
      V(i)=0
      if(i>=k)then
        V(i)=A(i,k)
      end if
      norm=norm+V(i)*V(i)
    end do
    norm=sqrt(norm)
    if(V(k)>0) then
      V(k)=V(k)+norm
    else
      V(k)=V(k)-norm
    end if
!He-he we need to calculate norm two times!!
    norm=0
    do i=1,N
      norm=norm+V(i)*V(i)
    end do
    do i=1,N
      do j=1,i
        H(i,j)=-2*V(i)*V(j)/norm
        H(j,i)=-2*V(i)*V(j)/norm
        if(i==j)then
          H(i,j)=H(i,j)+1;!Здесь добавляется матрица E
        end if
      end do
    end do
    call mmul(N,H,A)!A:=H*A то есть в конце разложения здесь получится матрица R
    call mmul(N,H,E)!E:=H*E то есть в конце разложения здесь получится матрица Q
  end do
  call trans(N,E)
!Конец разложения
end subroutine QR

subroutine eigen(N,Asav,E,L)
  integer N,i,j,k
  doubleprecision Asav(N,N),E(N,N),L(N,N),Q(N,N),R(N,N)
  do i=1,N
    do j=1,N
      L(i,j)=Asav(i,j)
      if(i==j) then
        E(i,j)=1
      else
        E(i,j)=0
      end if
    end do
  end do
!  print *,'--------L----------'
!  call outer(N,L)
  do i=1,100*N
    call QR(N,L,Q,R)
!    print *,'=====Q====='
!    call outer(N,Q)
!    print *,'=====R====='
!    call outer(N,R)
    call mmul_left(N,E,Q)
    call mmul(N,R,Q)
    do j=1,N
      do k=1,N
        L(j,k)=Q(j,k)
      end do
    end do
  end do
end subroutine eigen

subroutine outer(N,A)
  integer i,j,N
  double precision A(N,N)
  do 70, i=1,N                         
    print 80,(A(i,j), j=1,N)
70  continue
80  format(8F8.4)
end subroutine outer

subroutine outvec(N,V)
  integer i,N
  double precision V(N)
    print 80,(V(i), i=1,N)
70  continue
80  format(8F8.4)
end subroutine outvec

!Differential solver
!Subroutine call convention:
!call F(NEQ,t,y,ret_val)
subroutine rk4(F,y,NEQ,tstart,dt,nsteps)
  external F
  integer NEQ,i,nsteps
  doubleprecision k1(NEQ)
  doubleprecision k2(NEQ)
  doubleprecision k3(NEQ)
  doubleprecision k4(NEQ)
  doubleprecision y(NEQ)
  doubleprecision y1(NEQ)
  doubleprecision rv(NEQ)
  doubleprecision tstart,dt,t
  t=tstart
  do i=1,nsteps
!    t=tstart+(i-1)*dt
    call F(NEQ,t,y,rv)
    k1=dt*rv
    y1=y+0.5d0*k1
    call F(NEQ,t+0.5d0*dt,y1,rv)
    k2=dt*rv
    y1=y+0.5d0*k2
    call F(NEQ,t+0.5d0*dt,y1,rv)
    k3=dt*rv
    y1=y+k3
    call F(NEQ,t+dt,y1,rv)
    k4=dt*rv
    y=y+(1d0/6d0)*(k1+2d0*k2+2d0*k3+k4)
    t=t+dt
  end do
end subroutine rk4


!Subroutine call convention:
!call F(NEQ,t,y,ret_val)
subroutine rk4_step(F,y,NEQ,t,dt)
  external F
  integer NEQ,i,nsteps
  doubleprecision k1(NEQ)
  doubleprecision k2(NEQ)
  doubleprecision k3(NEQ)
  doubleprecision k4(NEQ)
  doubleprecision y(NEQ)
  doubleprecision y1(NEQ)
  doubleprecision rv(NEQ)
  doubleprecision dt,t
  call F(NEQ,t,y,rv)
  k1=dt*rv
  y1=y+0.5d0*k1
  call F(NEQ,t+0.5d0*dt,y1,rv)
  k2=dt*rv
  y1=y+0.5d0*k2
  call F(NEQ,t+0.5d0*dt,y1,rv)
  k3=dt*rv
  y1=y+k3
  call F(NEQ,t+dt,y1,rv)
  k4=dt*rv
  y=y+(1d0/6d0)*(k1+2d0*k2+2d0*k3+k4)
  t=t+dt
end subroutine rk4_step


subroutine rkf45_step(F,y,NEQ,t,dt)
  external F
  integer NEQ,i,nsteps
  doubleprecision k1(NEQ)
  doubleprecision k2(NEQ)
  doubleprecision k3(NEQ)
  doubleprecision k4(NEQ)
  doubleprecision k5(NEQ)
  doubleprecision k6(NEQ)
  doubleprecision y(NEQ)
  doubleprecision y1(NEQ)
  doubleprecision rv(NEQ)
  doubleprecision dt,t
  call F(NEQ,t,y,rv)
  k1=dt*rv
  
  y1=y+(1d0/4d0)*k1
  call F(NEQ,t+(1d0/4d0)*dt,y1,rv)
  k2=dt*rv
  
  y1=y+(3d0/32d0)*k1+(9d0/32d0)*k2
  call F(NEQ,t+(3d0/8d0)*dt,y1,rv)
  k3=dt*rv

  y1=y+(1932d0/2197d0)*k1-(7200d0/2197d0)*k2+(7296d0/2197d0)*k3
  call F(NEQ,t+(12d0/13d0)*dt,y1,rv)
  k4=dt*rv
  
  y1=y+(439d0/216d0)*k1-8d0*k2+(3680d0/513d0)*k3-(845d0/4104d0)*k4
  call F(NEQ,t+dt,y1,rv)
  k5=dt*rv
  
  y1=y-(8d0/27d0)*k1+2d0*k2-(3544d0/2565d0)*k3+(1859d0/4104d0)*k4-(11d0/40d0)*k5
  call F(NEQ,t+(1d0/2d0)*dt,y1,rv)
  k6=dt*rv
  
  y=y+(16d0/135d0)*k1+(6656d0/12825d0)*k3+(28561d0/56430d0)*k4-(9d0/50d0)*k5+(2d0/55d0)*k6
  t=t+dt
end subroutine rkf45_step


subroutine rkf45(F,y,NEQ,tstart,dt,nsteps)
  external F
  integer NEQ,i,nsteps
  doubleprecision k1(NEQ)
  doubleprecision k2(NEQ)
  doubleprecision k3(NEQ)
  doubleprecision k4(NEQ)
  doubleprecision k5(NEQ)
  doubleprecision k6(NEQ)
  doubleprecision y(NEQ)
  doubleprecision y1(NEQ)
  doubleprecision rv(NEQ)
  doubleprecision dt,t,tstart
  t=tstart
  do i=1,nsteps
    call F(NEQ,t,y,rv)
    k1=dt*rv
  
    y1=y+(1d0/4d0)*k1
    call F(NEQ,t+(1d0/4d0)*dt,y1,rv)
    k2=dt*rv
  
    y1=y+(3d0/32d0)*k1+(9d0/32d0)*k2
    call F(NEQ,t+(3d0/8d0)*dt,y1,rv)
    k3=dt*rv

    y1=y+(1932d0/2197d0)*k1-(7200d0/2197d0)*k2+(7296d0/2197d0)*k3
    call F(NEQ,t+(12d0/13d0)*dt,y1,rv)
    k4=dt*rv
  
    y1=y+(439d0/216d0)*k1-8d0*k2+(3680d0/513d0)*k3-(845d0/4104d0)*k4
    call F(NEQ,t+dt,y1,rv)
    k5=dt*rv
  
    y1=y-(8d0/27d0)*k1+2d0*k2-(3544d0/2565d0)*k3+(1859d0/4104d0)*k4-(11d0/40d0)*k5
    call F(NEQ,t+(1d0/2d0)*dt,y1,rv)
    k6=dt*rv
  
    y=y+(16d0/135d0)*k1+(6656d0/12825d0)*k3+(28561d0/56430d0)*k4-(9d0/50d0)*k5+(2d0/55d0)*k6
    t=t+dt
  end do
end subroutine rkf45


end module dfl