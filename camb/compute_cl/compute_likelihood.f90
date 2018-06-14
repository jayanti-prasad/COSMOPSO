subroutine compute_likelihood_new(p,y)
  use defs 
  Use wmap_likelihood_9yr
  Use wmap_options
  Use wmap_util
  implicit none 
  double precision, intent(in) :: p(:)
  double precision, intent(out) :: y 
!  integer, parameter :: tmin=2,tmax=2500
  

  interface 
     function test_wmap9(cl_tt, cl_te, cl_ee, cl_bb)
       Use wmap_likelihood_9yr
       Use wmap_options
       Use wmap_util
       implicit none
       double precision, intent(in) :: cl_tt(:), cl_te(:), cl_ee(:), cl_bb(:)
       double precision :: test_wmap9
     end function test_wmap9
     
     function planck(Cl_TT,Cl_EE,Cl_BB,Cl_TE)
       use clik
       use defs
       implicit none
       double precision, intent(in) :: Cl_TT(:),Cl_EE(:),Cl_BB(:),Cl_TE(:)
       double precision :: planck
     end function planck
     
     
     subroutine compute_cl(p,tmin,tmax,tt,te,ee,bb)
       use CAMB
       implicit none
       integer, intent(in) :: tmin,tmax
       real(8), intent(in) :: p(:)
       real(8), dimension(:), intent(out)  :: tt,te,ee,bb
     end   subroutine compute_cl
  end interface
  
  double precision, Dimension(:), Allocatable :: cl_tt, cl_te, cl_ee, cl_bb, cl_pp
  double precision, dimension(1:tmax)   :: c_tt,c_te,c_ee,c_bb
  double precision :: y_wmap9,y_planck 
  integer :: l
!  double precision, dimension(6) :: x0

  y_wmap9 = 0.d0
  y_planck = 0.d0 

  
!   x0=(/0.02264d0,0.1138d0,70.0d0,2.41d0,0.9720d0,0.089d0/)
 
   do l = 1, tmax
         c_tt(l) = 0.d0
         c_ee(l) = 0.d0
         c_te(l) = 0.d0
         c_bb(l)  = 0.d0
   end do
 
   call compute_cl(p,tmin,tmax,c_tt,c_te,c_ee,c_bb)
   
   if (USE_WMAP) then 
      
      Allocate(cl_tt(ttmin:ttmax))
      Allocate(cl_te(ttmin:ttmax))
      Allocate(cl_ee(ttmin:ttmax))
      Allocate(cl_bb(ttmin:ttmax))
      Allocate(cl_pp(ttmin:ttmax))
      
      
      do l=ttmin,ttmax
         cl_tt(l)=c_tt(l)
         cl_ee(l)=c_ee(l)
         cl_te(l)=c_te(l)
         cl_bb(l)=c_bb(l)
      end do
      
      y_wmap9 = test_wmap9(cl_tt, cl_te, cl_ee, cl_bb) 
      
      deallocate(cl_tt)
      deallocate(cl_ee)
      deallocate(cl_bb)
      deallocate(cl_te)
      
   end if
   
   if (USE_PLANCK) then 
      y_planck=planck(C_TT,C_EE,C_BB,C_TE)
   end if
   
   
   y =   y_wmap9 - y_planck
   
   !write(*,'(9(F18.6))')y_wmap9,y_planck,y,(p(l),l=1,6)
   
 end subroutine compute_likelihood_new
 
 
 

