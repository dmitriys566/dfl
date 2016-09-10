program example
use dfl
implicit none
doublecomplex z
doubleprecision nu
call istart(11)!This is need to gauss integrate. You must call this in program every time you want to integrate.
call gen_fact(100)
z=(1d0,1d0)
nu=5d0
print *,'BesselJ=',CBesselJ(nu,z)
print *,'Integral=',dint(f,0d0,1d0,10)
CONTAINS
doubleprecision function f(x)
doubleprecision x
f=x**3
end function f
end program example