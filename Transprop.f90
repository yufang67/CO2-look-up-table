MODULE Transprop
!
!
!     USE def_variables     
     USE def_constants
     USE properties, ONLY: satprop,heat_cap_p,heat_cap_v,entropy,&
&                          inter_energy, dedr_T
!     USE interp_functions
!     USE non_linear_solvers
!      
     IMPLICIT NONE
!
     PRIVATE
     PUBLIC  CO2visco, CO2conduc, cpCO2, cvCO2, entropyCO2,CO2conduc2phase,dedrCO2
!
!
     CONTAINS
!
!===============================================================================
!
SUBROUTINE CO2visco(vis_out, v_in, T_in, x_in, p_in, flag_loca) 
!
!===============================================================================
!
!
! Input: T_in      (the temperature)
!        v_in      (the specific volume v_in) 
!        flag_loca (the location flag) 
! Output: 
!         vis_out (viscosity)
!
!===============================================================================
REAL(pr),INTENT(OUT)  :: vis_out
INTEGER,INTENT(IN)   :: flag_loca
REAL(pr),INTENT(IN)   :: v_in, T_in, x_in,p_in
!
!Local Variables
INTEGER :: i, j,flag
REAL(pr) :: eta0, eta1, delta_etar, delta_etac
REAL(pr) :: b_etas, T_s, sigma, NA, eta_tL, Tr, rhor,rho_tL
REAL(pr) :: delta,qual,press,temp, rho
!REAL(pr) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
REAL(pr) :: rhol, rhov,rhorl,rhorv,delta_etarl,delta_etarv
REAL(pr) :: dummy1,dummy2,dummy3, vvsat, vlsat
!
flag = flag_loca
!
!   IF ( (T_in <= 304.99) .AND. (rho>=300.0) .AND. (rho<=700.0) ) THEN
!      temp = 304.99
!   ENDIF
rho  = 1.0/v_in
temp = T_in
!
!
eta0 = 1.0055*sqrt(temp) / ( av_00 + av_01*temp**(1.0/6.0) + av_02*exp(av_03*temp**(1.0/3.0)) + &
&                            (av_04+av_05*temp**(1.0/3.0))/exp(temp**(1.0/3.0)) + av_06*sqrt(temp) )
!
!
T_s    = temp/200.76
sigma  = 3.78421e-10_pr 
NA     = 6.022140857e23_pr
b_etas = bv_00 +&
&        bv_01/(T_s**(tv_01))+&
&        bv_02/(T_s**(tv_02))+&
&        bv_03/(T_s**(tv_03))+&
&        bv_04/(T_s**(tv_04))+&
&        bv_05/(T_s**(tv_05))+&
&        bv_06/(T_s**(tv_06))+&
&        bv_07/(T_s**(tv_07))+&
&        bv_08/(T_s**(tv_08))
!
!
eta1 = eta0 * b_etas * sigma**3.0 * NA / W 
!
!
rho_tL = 1178.56 ! 1178.46
Tr     = T_in / 216.592
IF (flag /= 5) THEN
   rhor   = rho  / rho_tL
   eta_tL = rho_tL**(2.0/3.0) * sqrt(8.3144598*T_tri) / ( W**(1.0/6.0) * NA**(1.0/3.0) )
!print*,'eta_tL', eta_tL
!
!
   delta_etar = eta_tL * ( cv_01*Tr*rhor**(3.0) + (rhor**(2.0)+rhor**(gamm_v)) / (Tr-cv_02) ) * 1e3_pr 
!
!
   delta_etac = 0.0
!
!
   vis_out = (eta0 +  rho*eta1 + delta_etar + delta_etac) * 1e-3_pr
!
!print*, eta0, rho*eta1, delta_etar
!
ELSE 
   CALL satprop(3, p_in, dummy1, vvsat, vlsat, dummy2, dummy3) 
   rhov  = 1.0 /vvsat
   rhol  = 1.0 /vlsat
   temp = T_in
!
   rhorl   = rhol  / rho_tL
   rhorv   = rhov  / rho_tL
!
   eta_tL = rho_tL**(2.0/3.0) * sqrt(8.3144598*T_tri) / ( W**(1.0/6.0) * NA**(1.0/3.0) )
!print*,'eta_tL', eta_tL
!
!
   delta_etarl = eta_tL * ( cv_01*Tr*rhorl**(3.0) + (rhorl**(2.0)+rhorl**(gamm_v)) / (Tr-cv_02) ) * 1e3_pr
   delta_etarv = eta_tL * ( cv_01*Tr*rhorv**(3.0) + (rhorv**(2.0)+rhorv**(gamm_v)) / (Tr-cv_02) ) * 1e3_pr
!
!
   delta_etac = 0.0
!
!
   vis_out = (eta0 +  rhol*eta1 + delta_etarl + delta_etac) * 1e-3_pr * (1.0-x_in) + &
&            (eta0 +  rhov*eta1 + delta_etarv + delta_etac) * 1e-3_pr * x_in
!
!
ENDIF
!######### CHECK #########
! viscosity
  IF ( (vis_out /= vis_out) .OR. (vis_out < 0.0) .OR. (vis_out >1.0)) THEN
     print*, 'viscosity error',vis_out
     print*, 'CO2 TRANSPROP, Tin= ', T_in, 'vin= ',v_in, 'flag',flag
     vis_out = 0.0
!     STOP
  ENDIF
END SUBROUTINE CO2visco
!
!
!
!===============================================================================
!
SUBROUTINE CO2conduc(T_in, v_in, conduc_out)
!
!===============================================================================
!
!
! Input: T_in      (the temperature)
!        v_in      (the specific volume v_in) 
!        flag_loca (the location flag) 
!        x_in      (the quality for two-phase) 
! Output: 
!         conduc_out (the thermal conductivity)
!
!===============================================================================
REAL(pr),INTENT(OUT)  :: conduc_out
!INTEGER,INTENT(IN)   :: flag_loca
REAL(pr),INTENT(IN)   :: v_in, T_in!, x_in
!
!Local Variables
INTEGER :: i, j,flag
REAL(pr) :: lamb0, delta_lamb, delta_lambc
REAL(pr) :: b_etas, T_s, sigma, NA, eta_tL, Tr, rhor,rho_tL
REAL(pr) :: delta,qual,press,temp, rho
!REAL(pr) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
REAL(pr) :: rhol, rhov,rhorl,rhorv,delta_etarl,delta_etarv

!
!flag = flag_loca
rho  = 1.0/v_in
temp = T_in
!
!
Tr   = temp/T_cr
!
lamb0 = sqrt(Tr) / ( lc_00 + lc_01/Tr + lc_02/Tr**2.0 + lc_03/Tr**3.0) * 1e-3   !W/(m*K)
!
!
rhor = rho/rho_cr
!
delta_lamb = (b1c_01 + b2c_01*Tr) * rhor      + &                !W/(m*K)
&            (b1c_02 + b2c_02*Tr) * rhor**2.0 + &
&            (b1c_03 + b2c_03*Tr) * rhor**3.0 + &
&            (b1c_04 + b2c_04*Tr) * rhor**4.0 + &
&            (b1c_05 + b2c_05*Tr) * rhor**5.0 + &
&            (b1c_06 + b2c_06*Tr) * rhor**6.0
!
!
! delta_lambc has empirical formulation
delta_lambc = ( -17.47 - 44.88*(Tr-1.0) ) / ( 0.8563 - &
&             exp(8.865*(Tr-1.0) + 4.16*(rhor-1.0)**2.0 + 2.302*(Tr-1.0)*(rhor-1.0) - (rhor-1.0)**3.0) -&
&             0.4503*(rhor-1.0) - 7.197*(Tr -1.0)) * 1e-3    ! W/(mK)
!
!
conduc_out = lamb0 + delta_lamb + delta_lambc
!
!
!
END SUBROUTINE CO2conduc
!
!
!
!===============================================================================
!
SUBROUTINE CO2cp(cp_out, v_in, u_in, x_in, flag_loca)
!
!===============================================================================
!
USE def_variables
USE interp_functions
!
REAL(pr), INTENT(OUT) :: cp_out
REAL(pr), INTENT(IN)  :: u_in, v_in,x_in
INTEGER, INTENT(IN)   :: flag_loca
!
! Local Variables
INTEGER :: i, j
REAL(pr) :: x_mesh_test, delta_u, delta_v, v_min_li, v_max_li,v_max_li_log, &
!&           T_down_left, T_down_right, T_up_left, T_up_right, &
!&           p_down_left, p_down_right, p_up_left, p_up_right, &
&           cp_down_left, cp_down_right, cp_up_left, cp_up_right, &
&           DEN, A, B, D, E
REAL(pr) :: delta
!
!
IF (flag_loca == 1) THEN ! ### LL ###
  delta_u = y_mesh_LL(2) - y_mesh_LL(1)
  i = INT((u_in - y_mesh_LL(1))/delta_u) + 1
!
v_min_li = 0_pr
v_max_li = 0_pr
v_max_li_log = 0_pr
!
!
  DO j = 1, ord_spline + 1
     v_min_li = v_min_li + spline_pmax_LL(ord_spline+2-j,i) * u_in**(j-1)
!Evaluate Vmin e Vmax, on the row i(high already identified)
     v_max_li = v_max_li + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!     v_max_li_log = v_max_li_log + spline_Lsat_LL(ord_spline+2-j,i) * u_in**(j-1)
!     v_max_li     = 10_pr ** v_max_li_log
!Approximated through the spline 
  ENDDO
!
!
  x_mesh_test = ((x_mesh_LL(MMM_LL) - x_mesh_LL(1)) / &
& (v_max_li - v_min_li)) * (v_in - v_min_li) + x_mesh_LL(1)
!  
  delta_v = x_mesh_LL(2) - x_mesh_LL(1)
  j = INT((x_mesh_test - x_mesh_LL(1))/delta_v) + 1
!
!
  cp_down_left  = cpcp_LL(i,j)
  cp_down_right = cpcp_LL(i,j+1)
  cp_up_left    = cpcp_LL(i+1,j)
  cp_up_right   = cpcp_LL(i+1,j+1)
!
! Bilinear Interpolation:
!
  DEN = (x_mesh_LL(j+1)- x_mesh_LL(j))*(y_mesh_LL(i)- y_mesh_LL(i+1))
  A   = (x_mesh_LL(j+1)- x_mesh_test) *(y_mesh_LL(i)- u_in)
  B   = (x_mesh_test   - x_mesh_LL(j))*(y_mesh_LL(i)- u_in)
  D   = (x_mesh_LL(j+1)- x_mesh_test) *(u_in - y_mesh_LL(i+1))
  E   = (x_mesh_test   - x_mesh_LL(j))*(u_in - y_mesh_LL(i+1))
!
  cp_out = (A/DEN) * cp_up_left + (B/DEN) * cp_up_right + (D/DEN) * cp_down_left + &
&    (E/DEN) * cp_down_right
!
!
!
ELSEIF (flag_loca == 2) THEN ! ### LH ###
  delta_u = y_mesh_LH(2)     - y_mesh_LH(1)
  i       = INT((u_in - y_mesh_LH(1))/delta_u) + 1
!
  v_min_li     = 0_pr
  v_max_li     = 0_pr
  v_max_li_log = 0_pr
!
  DO j = 1, ord_spline + 1
    v_min_li     = v_min_li     + spline_pmax_LH(ord_spline+2-j,i) * u_in**(j-1)
    v_max_li     = v_max_li     + spline_Lsat_LH(ord_spline+2-j,i) * u_in**(j-1)
!   v_max_li_log = v_max_li_log + spline_Lsat_LH(ord_spline+2-j,i) * u_in**(j-1)
!   v_max_li     = 10_pr ** v_max_li_log
  ENDDO
!
  x_mesh_test = ((x_mesh_LH(MMM_LH) - x_mesh_LH(1)) / (LOG10(v_max_li) - LOG10(v_min_li)))&
& * (LOG10(v_in) - LOG10(v_min_li)) +  x_mesh_LH(1)
!x_mesh_test = ((x_mesh_LH(MMM_LH) - x_mesh_LH(1)) / &
!& (v_max_li - v_min_li)) * (v_in - v_min_li) + x_mesh_LH(1)
!
  delta_v = x_mesh_LH(2)     - x_mesh_LH(1)
  j       = INT((x_mesh_test - x_mesh_LH(1))/delta_v) + 1
!
  cp_down_left  = cpcp_LH(i,j)
  cp_down_right = cpcp_LH(i,j+1)
  cp_up_left    = cpcp_LH(i+1,j)
  cp_up_right   = cpcp_LH(i+1,j+1)
!
! Bilinear Interpolation:
!
  DEN = (x_mesh_LH(j+1)- x_mesh_LH(j))* (y_mesh_LH(i)- y_mesh_LH(i+1))
  A   = (x_mesh_LH(j+1)- x_mesh_test) * (y_mesh_LH(i)- u_in)
  B   = (x_mesh_test   - x_mesh_LH(j))* (y_mesh_LH(i)- u_in)
  D   = (x_mesh_LH(j+1)- x_mesh_test) * (u_in - y_mesh_LH(i+1))
  E   = (x_mesh_test   - x_mesh_LH(j))* (u_in - y_mesh_LH(i+1))
!
  cp_out = (A/DEN) * cp_up_left + (B/DEN) * cp_up_right + (D/DEN) * cp_down_left + &
&    (E/DEN) * cp_down_right
!
!
!
!
ELSEIF (flag_loca == 3) THEN !### R ###
  delta        = (x_mesh_max - x_mesh_min)/(MMM_R-1)
  x_mesh_R     = x_mesh_min + (/(i*delta, i=0,MMM_R-1)/)
     
  CALL Lin_int(cp_out, cp_out, cp_out, u_in, v_in, NNN_R, MMM_R, x_mesh_R, y_mesh_R, &
&         spline_pmin, spline_Vsat, cpcp_R, cpcp_R, cpcp_R)
!
!
!
!
ELSEIF (flag_loca ==4) THEN ! ### HT ###
     delta       = (x_mesh_max - x_mesh_min)/(MMM_HT-1)
     x_mesh_HT   =  x_mesh_min + (/(i*delta, i=0,MMM_HT-1)/)
     CALL Lin_int_Log10(cp_out, cp_out, cp_out, u_in, v_in, NNN_HT, MMM_HT, x_mesh_HT, y_mesh_HT, &
&         spline_right_HT, spline_left_HT, cpcp_HT, cpcp_HT, cpcp_HT)
!
!
!
!
ELSE
  print*, 'two-phase cp'
  STOP
!
!
!
!
ENDIF 
END SUBROUTINE CO2cp
!
!===============================================================================
!
SUBROUTINE cpCO2(cp_out, v_in, vp_out, x_in, T_in, p_in, flag_loca)
!
!===============================================================================
!
!   USE def_variables
!   USE interp_functions

REAL(pr), INTENT(OUT) :: cp_out,vp_out
REAL(pr), INTENT(IN)  :: v_in,x_in, T_in, p_in
INTEGER, INTENT(IN)   :: flag_loca
!
REAL(pr):: heatcp,cpv,cpl
REAL(pr):: vvsat,vlsat
REAL(pr):: dummy1,dummy2,dummy3
REAL(pr)::vl
!
!
IF (flag_loca /= 5) THEN 

   CALL heat_cap_p(T_in,v_in,heatcp)
   cp_out = heatcp

ELSE 
   CALL satprop(3, p_in, dummy1, vvsat, vlsat, dummy2, dummy3)
!   vl   = (v_in - x_in*vp_in) / (1.0-x_in)

   CALL heat_cap_p(T_in,vvsat,cpv)
   CALL heat_cap_p(T_in,vlsat   ,cpl)

   cp_out = x_in*cpv + (1.0-x_in)*cpl
   vp_out  = vvsat

ENDIF

!######### CHECK #########

IF ( (cp_out /= cp_out) .OR. (cp_out < 0.0) .OR. (cp_out >1e6)) THEN
   print*, 'Cp error',cp_out
   print*, 'CO2 TRANSPROP, Tin= ', T_in,'p_in= ',p_in, 'vin= ',v_in, 'flag',flag_loca
   cp_out = 0.0
!     STOP
ENDIF

END SUBROUTINE cpCO2
!
!===============================================================================
!
SUBROUTINE cvCO2(cv_out, v_in, vp_out, x_in, T_in, p_in, flag_loca)
!
!===============================================================================
!
!   USE def_variables
!   USE interp_functions

REAL(pr), INTENT(OUT) :: cv_out,vp_out
REAL(pr), INTENT(IN)  :: v_in,x_in, T_in, p_in
INTEGER, INTENT(IN)   :: flag_loca
!
REAL(pr):: heatcv,cvv,cvl
REAL(pr):: vvsat,vlsat
REAL(pr):: dummy1,dummy2,dummy3
!
REAL(pr)::vl
!
IF (flag_loca /= 5) THEN 

   CALL heat_cap_v(T_in,v_in,heatcv)
   cv_out = heatcv

ELSE 
   CALL satprop(3, p_in, dummy1, vvsat, vlsat, dummy2, dummy3)
!   vl   = (v_in - x_in*vp_in) / (1.0-x_in)

   CALL heat_cap_v(T_in,vvsat,cvv)
   CALL heat_cap_v(T_in,vlsat,cvl)

   cv_out = x_in*cvv + (1.0-x_in)*cvl
   vp_out  = vvsat

ENDIF

!######### CHECK #########

IF ( (cv_out /= cv_out) .OR. (cv_out < 0.0) .OR. (cv_out >1e6)) THEN
   print*, 'Cv error',cv_out
   print*, 'CO2 TRANSPROP, Tin= ', p_in, 'vin= ',v_in, 'flag',flag_loca
   cv_out = 0.0
!     STOP
ENDIF

END SUBROUTINE cvCO2
!
!
!===============================================================================
!
SUBROUTINE entropyCO2(s_out, v_in, vp_out, x_in, T_in, p_in, flag_loca)
!
!===============================================================================
!
!   USE def_variables
!   USE interp_functions

REAL(pr), INTENT(OUT) :: s_out,vp_out
REAL(pr), INTENT(IN)  :: v_in,x_in, T_in, p_in
INTEGER, INTENT(IN)   :: flag_loca
!
REAL(pr):: s,sv,sl
REAL(pr):: vvsat,vlsat
REAL(pr):: dummy1,dummy2,dummy3
!
REAL(pr)::vl
!
IF (flag_loca /= 5) THEN 

   CALL entropy(T_in,v_in,s)
   s_out = s

ELSE 
   CALL satprop(3, p_in, dummy1, vvsat, vlsat, dummy2, dummy3)
!   vl   = (v_in - x_in*vp_in) / (1.0-x_in)

   CALL entropy(T_in,vvsat,sv)
   CALL entropy(T_in,vlsat,sl)

   s_out = x_in*sv + (1.0-x_in)*sl
   vp_out  = vvsat

ENDIF

!######### CHECK #########

IF ( (s_out /= s_out) .OR. (s_out >1e6)) THEN
   print*, 'entropy error',s_out
   print*, 'CO2 TRANSPROP, Tin= ', p_in, 'vin= ',v_in, 'flag',flag_loca
   s_out = 0.0
!     STOP
ENDIF

END SUBROUTINE entropyCO2


!===============================================================================
!
SUBROUTINE CO2conduc2phase(lambda_out, v_in, vp_out, x_in, T_in, p_in, flag_loca)
!
!===============================================================================
!
!   USE def_variables
!   USE interp_functions

REAL(pr), INTENT(OUT) :: lambda_out,vp_out
REAL(pr), INTENT(IN)  :: v_in,x_in, T_in, p_in
INTEGER, INTENT(IN)   :: flag_loca
!
REAL(pr):: lam,lamv,laml
REAL(pr):: vvsat,vlsat
REAL(pr):: dummy1,dummy2,dummy3
!
REAL(pr)::vl
!
IF (flag_loca /= 5) THEN 

   CALL CO2conduc(T_in,v_in,lam)
   lambda_out = lam

ELSE 
   CALL satprop(3, p_in, dummy1, vvsat, vlsat, dummy2, dummy3)
!   vl   = (v_in - x_in*vp_in) / (1.0-x_in)

   CALL CO2conduc(T_in,vvsat,lamv)
   CALL CO2conduc(T_in,vlsat,laml)

   lambda_out = x_in*lamv + (1.0-x_in)*laml
   vp_out  = vvsat

ENDIF

!######### CHECK #########

IF ( (lambda_out /= lambda_out) .OR. (lambda_out >1e6)) THEN
   print*, 'lambda error',lambda_out
   print*, 'CO2 TRANSPROP, Tin= ', p_in, 'vin= ',v_in, 'flag',flag_loca
   lambda_out = 0.0
!     STOP
ENDIF

END SUBROUTINE CO2conduc2phase
!
!
!===============================================================================
!
SUBROUTINE dedrCO2(dedr_out, v_in, vp_out, x_in, T_in, p_in, flag_loca)
!
!===============================================================================
!
!   USE def_variables
!   USE interp_functions

REAL(pr), INTENT(OUT) :: dedr_out,vp_out
REAL(pr), INTENT(IN)  :: v_in,x_in, T_in, p_in
INTEGER, INTENT(IN)   :: flag_loca
!
REAL(pr):: dedr
REAL(pr):: vvsat,vlsat,eev,eel
REAL(pr):: dummy1,dummy2,dummy3
!
REAL(pr)::vl
!
IF (flag_loca /= 5) THEN 

   CALL dedr_T(T_in,v_in,dedr)
   dedr_out = dedr

ELSE 
   CALL satprop(3, p_in, dummy1, vvsat, vlsat, dummy2, dummy3)
   CALl inter_energy(T_in,vvsat,eev)
   CALl inter_energy(T_in,vlsat,eel)
!   vl   = (v_in - x_in*vp_in) / (1.0-x_in)


   dedr_out = -v_in*v_in*(eev-eel)/(vvsat-vlsat)
   vp_out  = vvsat

ENDIF

!######### CHECK #########

IF (dedr_out /= dedr_out) THEN
   print*, 'dedr error',dedr_out
   print*, 'CO2 TRANSPROP, Tin= ', p_in, 'vin= ',v_in, 'flag',flag_loca
   dedr_out = 0.0
!     STOP
ENDIF

END SUBROUTINE dedrCO2
!
!
!
END MODULE Transprop
