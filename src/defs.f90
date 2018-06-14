module defs 
  use camb 
  use clik
  integer, parameter :: tmin=2,tmax=2500,max_len=512,ntimes_max=1000,ndim_max=50,npart_max=100 
  double precision, parameter :: minusone = -1.d0
  
  double precision, dimension(npart_max,ndim_max) :: x, v  
  double precision, dimension(npart_max) :: f
  double precision, dimension(npart_max,ndim_max) ::  xpbest 
  double precision, dimension(npart_max) :: pbest 
  double precision, dimension(ntimes_max) :: gbest 
  double precision, dimension(ntimes_max,ndim_max) :: xgbest 
  double precision, dimension(ntimes_max,ndim_max) :: xcbest   
  double precision, dimension(ndim_max) :: xmin,xmax,vmin,vmax,xdef 
 
  integer :: npart, ndim, nstops
  type(clik_object)  :: pself_camspec,pself_commander, pself_lowl
  character(len=max_len) :: planck_data
  integer :: seed_uran,ierror
  double precision :: eps,w,c_gbest,c_pbest,c_cbest,c_vmax 
  LOGICAL :: use_WMAP, use_PLANCK
    
end module defs
