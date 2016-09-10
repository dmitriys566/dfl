program ode
use dfl
implicit none
integer NEQ,i
doubleprecision,allocatable :: y(:)
doubleprecision t
NEQ=1
allocate(y(NEQ))
y(1)=0d0
!y(2)=1d0
t=0d0
print *,t,y(1),abs(y(1)-tan(t))
do while(t<1.4d0)
  call rkf45_step(F,y,NEQ,t,0.0001d0)
  print *,t,y(1),abs(y(1)-tan(t))
!  print *,t,y(1),abs(y(1)-exp(t))
end do
CONTAINS
subroutine F(NEQ,t,y,rv)
  integer NEQ
  doubleprecision y(NEQ),rv(NEQ),t
  rv(1)=y(1)**2+1d0
!  rv(1)=y(1)
end subroutine F
end program ode