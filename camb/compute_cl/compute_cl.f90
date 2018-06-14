subroutine compute_cl(X,ttmin,ttmax,cl_tt,cl_te,cl_ee,cl_bb)
  use CAMB
  implicit none
  double precision, intent(in) :: X(:)
  integer, intent(in) :: ttmin,ttmax
  real(8), dimension(:), intent(out)  :: cl_tt,cl_te,cl_ee,cl_bb
  real :: Cls(ttmin:ttmax,1:4)
  type(CAMBparams)  P 
  integer l,error,il
  logical :: GC_conventions
  real(dl) ::  output_factor    
  
  call CAMB_SetDefParams(P)
  
 ! P%omegab  = 0.02264 ! actually ombh2
 ! P%omegac  = 0.1138   !actually och2 
 ! P%H0      = 70.0 ! 67.04346
 ! P%InitPower%ScalarPowerAmp(1) = 2.41e-9
 ! P%InitPower%an(1)=0.9720
 ! P%Reion%optical_depth=0.089  


 P%omegab  = X(1)
 P%omegac  = X(2) 
 P%H0      = X(3) 
 P%InitPower%ScalarPowerAmp(1) = X(4)*1.0e-9

! convert log[10^10A_s] to 1.0e-9
!  P%InitPower%ScalarPowerAmp(1) = (exp(x(4))/10.0) * 1.0e-9 

 P%InitPower%an(1)=X(5)
 P%Reion%optical_depth=X(6)

! do not change anything below ! 

  output_factor = 7.42835025e12

  P%InitPower%k_0_scalar = 0.002
  P%InitPower%k_0_tensor = 0.002

  P%omegab  = P%omegab /(P%H0/100)**2
  P%omegac  = P%omegac /(P%H0/100)**2
  P%omegan  = p%omegan /(P%H0/100)**2

  P%omegav  = 1.d0 -   (P%omegab+  P%omegac+ p%omegan)

  P%AccuratePolarization = .true.

  P%Reion%redshift = 10
  P%Reion%delta_redshift = 0.5_dl
  P%TCMB = 2.72548
  P%yhe  =  0.24
  P%Num_Nu_massless=3.04
 
  P%Max_l=2600
  P%Max_eta_k=4000
  P%Max_l_tensor=1500
  P%Max_eta_k_tensor=3000
 
   P%Reion%Reionization = .true.
   P%Reion%use_optical_depth = .TRUE.
   P%AccurateReionization = .true.
   P%AccurateBB=.true.

   P%WantCls=.TRUE.
   P%WantTransfer=.FALSE.
   P%WantScalars=.TRUE.
   P%WantTensors=.FALSE.
   P%WantVectors=.FALSE.
   P%DoLensing=.TRUE.

  call CAMB_GetResults(P, error)

  call CAMB_GetCls(Cls,ttmax,1,GC_conventions)

  do il=1, ttmax
     cl_tt(il)=0.d0
     cl_ee(il)=0.d0
     cl_bb(il)=0.d0
     cl_te(il)=0.d0
  end do 
 
  do il=ttmin,ttmax
     cl_tt(il)=output_factor * Cls(il,1)
     cl_ee(il)=output_factor * Cls(il,2)
     cl_bb(il)=output_factor * Cls(il,3)
     cl_te(il)=output_factor * Cls(il,4)
 !    write(300,*) il, cl_tt(il) 
 end do
  
end  subroutine compute_cl 


