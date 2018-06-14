program test_compute_cl
   implicit none

   interface 
  subroutine compute_cl(X,ttmin,ttmax,cl_tt,cl_te,cl_ee,cl_bb)
  use CAMB
  implicit none
  double precision, intent(in) :: X(:)
  integer, intent(in) :: ttmin,ttmax
  real(8), dimension(:), intent(out)  :: cl_tt,cl_te,cl_ee,cl_bb
  end subroutine compute_cl
  end interface

  integer :: i
  double precision, dimension(6) :: X 
  integer, parameter :: ttmin=1,ttmax=2500 
  double precision, dimension(ttmin:ttmax) ::  cl_tt,cl_te,cl_ee,cl_bb
  double precision :: tmp1=0.0d0, tmp2=0.0d0

  X  = (/0.022640d0, 0.11380d0,70.0d0,2.41d0,0.9720d0,0.0890d0/)

 call compute_cl(X,ttmin,ttmax,cl_tt,cl_te,cl_ee,cl_bb)

 do i  = 2, ttmax 
  write(*,'(I6,5(E14.6))') i,cl_tt(i),cl_ee(i),cl_te(i),tmp1,tmp2
 end do


end program test_compute_cl
