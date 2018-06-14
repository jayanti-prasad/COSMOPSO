function planck(Cl_TT,Cl_EE,Cl_BB,Cl_TE)
  use clik
  use defs  
  implicit none
  double precision, intent(in) :: Cl_TT(:),Cl_EE(:),Cl_BB(:),Cl_TE(:)  
  double precision :: planck
  character(len=max_len) ::  CAMspec,commander,lowLike
  character(len=12), dimension(:), pointer :: names
  integer :: numnames, nl, i, j, l,ncl,n1,total_numcls
  integer, dimension(6) :: has_cl,lmax  
  double precision, dimension(:), allocatable :: cl_and_pars
  double precision, dimension(14) :: p
  double precision, dimension(2)  :: lkl
  double precision ::  lkl_camspec,lkl_commander,lkl_Lowl 
  
  interface 
     subroutine CAMspecpara(p)
       implicit none
       double precision, intent(out)  :: p(:)
     end  subroutine CAMspecpara
  end interface

  
  CAMspec=trim(planck_data)//'CAMspec_v6.2TN_2013_02_26_dist.clik/'
  commander=trim(planck_data)//'commander_v4.1_lm49.clik/'
  lowLike=trim(planck_data)//'lowlike_v222.clik/'

  ! comppute for commander 
  
!  call clik_init(pself,trim(CAMspec))

  call clik_get_has_cl(pself_camspec,has_cl)
  call clik_get_lmax(pself_camspec,lmax)
 
 ! write(*,*)"camspec: has_cl", has_cl 
 ! write(*,*)"camspec: lmax", lmax  

 
  numnames=clik_get_extra_parameter_names(pself_camspec,names)

  !write(*,*)"camspec extra params",numnames 

  ncl = total_numcls(lmax)

  nl = numnames + ncl
  
  allocate(cl_and_pars(0:nl))
  
  cl_and_pars = 0.d0
  
  cl_and_pars(0) = -1.d0
  cl_and_pars(1) = -1.d0
  
  do l = 2, ncl
     cl_and_pars(l) = dble(Cl_TT(l) *(2.0*pi/(l*(l+1))))
  end do
  
  call CAMspecpara(p)
  
  do  i = 1, 14
     cl_and_pars(ncl+i)=dble(p(i))
  end do
  
  lkl_camspec  = clik_compute(pself_camspec,cl_and_pars)
  
  deallocate(cl_and_pars)
  
! compute commander likelihood

  call clik_get_has_cl(pself_commander,has_cl)
  call clik_get_lmax(pself_commander,lmax)

 ! write(*,*)"commander: has_cl", has_cl
 ! write(*,*)"commander: lmax", lmax 
  
  numnames=clik_get_extra_parameter_names(pself_commander,names)
  ncl = total_numcls(lmax)
  nl = numnames + ncl
  
 ! write(*,*)"commande extra params",numnames 
  allocate(cl_and_pars(0:nl))
  
  cl_and_pars = 0.d0
  
  cl_and_pars(0) = -1.d0
  cl_and_pars(1) = -1.d0
  
  do l = 2, ncl
     cl_and_pars(l) = dble(Cl_TT(l) *(2.0*pi/(l*(l+1))))
  end do
  
  lkl_commander  = clik_compute(pself_commander,cl_and_pars)
  
  deallocate(cl_and_pars)
  
  ! compute lowl Likelihood 
  
  call clik_get_has_cl(pself_lowl,has_cl)
  call clik_get_lmax(pself_lowl,lmax)
  numnames=clik_get_extra_parameter_names(pself_lowl,names)
  
  ncl = total_numcls(lmax)
  nl = numnames + ncl
  allocate(cl_and_pars(0:nl))
  
  cl_and_pars = 0.d0
  
  do l =  0, 3
     
     cl_and_pars(0+l*lmax(1)) = -1.d0
     cl_and_pars(1+l*lmax(1)) = -1.d0
     
  end do
  
  do l = 2, lmax(1)
     cl_and_pars(l) = dble(Cl_TT(l) *(2.0*pi/(l*(l+1))))
     cl_and_pars(l+lmax(1)) = dble(Cl_EE(l) *(2.0*pi/(l*(l+1))))
     cl_and_pars(l+lmax(1)+lmax(2)) = 0.d0 !dble(Cl_BB(l) *(2.0*pi/(l*(l+1))))
     cl_and_pars(l+lmax(1)+lmax(2)+lmax(3)) = dble(Cl_TE(l) *(2.0*pi/(l*(l+1))))
  end do
  
  lkl_lowl = 0.d0 !clik_compute(pself_lowl,cl_and_pars)
  
  
  deallocate(cl_and_pars)
  
  !  write(*,*)"low Like: has_cl", has_cl
  !  write(*,*)"low lmax:", lmax

  
  !  write(*,*)"lowlike extra params",numnames 
  
  planck =  lkl_camspec +  lkl_commander+lkl_lowl 
  
  write(*,'(3(A12,F14.6))')"commander=",lkl_commander,"CAMSPEC=",lkl_camspec,"Lowl",lkl_lowl 
  
end function planck

integer function total_numcls(lmax)
  integer, intent(in), dimension(6) :: lmax 
  integer  i, ncl; 
 
  ncl  = 0
  do i =  1, 6
     if ( lmax(i) .gt. 0) then 
        ncl = ncl + lmax(i)
     end if
  end do
  total_numcls = ncl 
  
end function total_numcls

subroutine CAMspecpara(x)
  implicit none 
  double precision, intent(out)  :: x(:) 
  
  x(1)= 152.d0
  x(2)= 63.3d0
  x(3)= 117.0d0
  x(4)= 8.0d0
  x(5)= 27.2d0
  x(6)= 6.80d0
  x(7)= 0.916d0
  x(8)= 0.406d0
  x(9)= 0.6d0
  x(10)= 0.1000711D+01
  x(11)= 0.9962765D+00
  x(12)= 0.8832844D+00
  x(13)= 0.9D0
  x(14)= 0.4262047D+00
  
  !01  * A_ps_100  : the point source contribution at 100Ghz, \muK^2 at l=3000
  !02  * A_ps_143  : the point source contribution at 143Ghz, \muK^2 at l=3000
  !03  * A_ps_217  : the point source contribution at 217Ghz, \muK^2 at l=3000
  !04  * A_cib_143 : the clustered cib contribution at 143Ghz, \muK^2 at l=3000
  !05  * A_cib_217 : the clustered cib contribution at 217Ghz, \muK^2 at l=3000
  !06  * A_sz      : the tSZ contriibution at 143Ghz
  !07  * r_ps      : the correlation between A_ps_143 and A_ps_217
  !08  * r_cib     : the correlation between A_cib_143 and A_cib_217
  !09  * n_Dl_cib  : the slope of the clustered cib spectrum (Cl^cib ~ l^(2+n_Dl_cib)
  !10  * cal_100   : the relative calibration of the 100Ghz channel to the 143Ghz
  !11  * cal_217   : the relative calibration of the 217Ghz channel to the 143Ghz
  !12  * xi_sz_cib : the correlation between SZ and cib
  !13  * A_ksz     : The amplitude of the kSZ
  !14  * Bm_1_1    : The first eigenmode of the beam error for the 100Ghz channel
  
end subroutine CAMspecpara
