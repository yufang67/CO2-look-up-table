PROGRAM new
!
!

USE def_constants
USE def_variables
USE non_linear_solvers
USE properties
USE grid_functions
USE Grid
USE interp_functions
USE Interp_table
USE Transprop
!
IMPLICIT NONE
!
INTEGER :: j, i, N, flag, k, p,domain, flag_loca
! Number of nodes for the test
!LL
INTEGER, PARAMETER ::  NNN_test_LL = 5 * NNN_LL , MMM_test_LL = 5 * MMM_LL
!LH
INTEGER, PARAMETER ::  NNN_test_LH = 5 * NNN_LH , MMM_test_LH = 5 * MMM_LH
!R
INTEGER, PARAMETER ::  NNN_test_R  = 5 * NNN_R  , MMM_test_R  = 5 * MMM_R
!HT
INTEGER, PARAMETER ::  NNN_test_HT = 5 * NNN_HT , MMM_test_HT = 5 * MMM_HT
!
REAL(8) :: dtest_u, dtest_v, &
& res, time1, time2, p_guess, T_guess                
!
REAL(8) :: delta, a_out,out3
REAL(8) :: v_max_log!, v_sat_li, v_sat_left_log, v_sat_left, v_sat_log, &
!&          v_sat_right
!
!LL
REAL(8)                         :: d_y_test_LL, d_x_test_LL
REAL(8), DIMENSION(NNN_test_LL) :: y_test_LL
REAL(8), DIMENSION(MMM_test_LL) :: x_test_LL
!LH
REAL(8)                         :: d_y_test_LH, d_x_test_LH
REAL(8), DIMENSION(NNN_test_LH) :: y_test_LH
REAL(8), DIMENSION(MMM_test_LH) :: x_test_LH
!R
REAL(8)                         :: d_y_test_R, d_x_test_R
REAL(8), DIMENSION(NNN_test_R)  :: y_test_R
REAL(8), DIMENSION(MMM_test_R)  :: x_test_R
!HT
REAL(8)                         :: d_y_test_HT, d_x_test_HT
REAL(8), DIMENSION(NNN_test_HT) :: y_test_HT
REAL(8), DIMENSION(MMM_test_HT) :: x_test_HT
!
! Matrixes declaration for each region
!LL
REAL(8), DIMENSION(NNN_test_LL, MMM_test_LL) :: v_test_LL_log,v_test_LL, c_BL_LL,cp_BL_LL, &
&  p_BL_LL,  T_BL_LL, x_out_LL, p_LOIET_LL, T_LOIET_LL, c_LOIET_LL,cp_LOIET_LL, &
&  err_pBL_LL, err_TBL_LL, err_cBL_LL!, c_BC_LL, p_BC_LL,  T_BC_LL,  &
!&  err_pBC_LL, err_TBC_LL, err_cBC_LL
!
!LH
REAL(8), DIMENSION(NNN_test_LH, MMM_test_LH) :: v_test_LH_log,v_test_LH, c_BL_LH, cp_BL_LH,&
&  p_BL_LH,  T_BL_LH, x_out_LH, p_LOIET_LH, T_LOIET_LH, c_LOIET_LH, cp_LOIET_LH,&
&  err_pBL_LH, err_TBL_LH, err_cBL_LH!, c_BC_LH, p_BC_LH,  T_BC_LH,  &
!&  err_pBC_LH, err_TBC_LH, err_cBC_LH, v_test_LH_log
!
!R
REAL(8), DIMENSION(NNN_test_R, MMM_test_R) :: v_test_R_log, v_test_R, c_BL_R, cp_BL_R, &
&  p_BL_R,  T_BL_R, x_out_R, p_LOIET_R, T_LOIET_R, c_LOIET_R, cp_LOIET_R,&
&  err_pBL_R, err_TBL_R, err_cBL_R!, c_BC_R, p_BC_R,  T_BC_R,  &
!&  err_pBC_R, err_TBC_R, err_cBC_R, v_test_R_log
!HT
REAL(8), DIMENSION(NNN_test_HT, MMM_test_HT) :: v_test_HT_log, v_test_HT, c_BL_HT, cp_BL_HT,&
&  p_BL_HT,  T_BL_HT, x_out_HT, p_LOIET_HT, T_LOIET_HT, c_LOIET_HT,cp_LOIET_HT, &
&  err_pBL_HT, err_TBL_HT, err_cBL_HT!, c_BC_HT, p_BC_HT,  T_BC_HT,  &
!&  err_pBC_HT, err_TBC_HT, err_cBC_HT, v_test_HT_log
   
REAL(8) :: v_min_li, v_max_li, out_2, criticalc,energy,energyc,err_cp
!
!
!
!
!====================================================
! CONSTRUCTION OF THE GRID ON THE (u,v) DIAGRAM
!====================================================
!
CALL cpu_time(time1)
CALL MAKE_GRID() 
CALL cpu_time(time2)
!STOP
print*, '===================================================='
print*, 'TIME INTERVAL FOR GRID CONSTRUCTION = ', time2-time1
print*, '===================================================='
!
! Introducing test values:
!
dtest_u      = 5d1
dtest_v      = 5d-4
domain       = 4

SELECT CASE(domain)
CASE(1)
!
!
! Left Low LL
!
y_test_LL(1)           = y_mesh_LL(1) + dtest_u
y_test_LL(NNN_test_LL) = y_mesh_LL(NNN_LL) - dtest_u
d_y_test_LL            = (y_test_LL(NNN_test_LL) - y_test_LL(1))/(NNN_test_LL - 1)
y_test_LL              = y_test_LL(1) + (/(j*d_y_test_LL, j=0,NNN_test_LL-1)/)

x_test_LL(1)           = x_mesh_LL(1) + dtest_v
x_test_LL(MMM_test_LL) = x_mesh_LL(MMM_LL) - dtest_v
d_x_test_LL            = (x_test_LL(MMM_test_LL) - x_test_LL(1))/(MMM_test_LL - 1)      
x_test_LL              = x_test_LL(1) + (/(j*d_x_test_LL, j=0,MMM_test_LL-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'err_pBL_LL.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'err_TBL_LL.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 63,             FILE = 'err_cBL_LL.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'press_tables_LL.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 65,             FILE = 'Temp_tables_LL.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_LL-1  !For each cell of the original grid to test
   delta = y_mesh_LL(2)  - y_mesh_LL(1)
   i = INT((y_test_LL(k) - y_mesh_LL(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li = 0d0
   v_max_li = 0d0
   DO j = 1, ord_spline + 1
      v_min_li = v_min_li + spline_pmax_LL(ord_spline+2-j,i) * y_test_LL(k) ** (j - 1)
!      v_max_log  = v_max_log + spline_Lsat_LL(ord_spline+2-j,i) * y_test_LL(k) ** (j - 1)
!      v_max_li   = 10d0 ** v_max_log
      v_max_li = v_max_li + spline_Lsat_LL(ord_spline+2-j,i) * y_test_LL(k) ** (j - 1)  !I found the min and max volumes correspondent to that energy level
   ENDDO
   DO p = 1, MMM_test_LL-1
   delta = x_mesh_LL(2)  - x_mesh_LL(1)
     j = INT((x_test_LL(p) - x_mesh_LL(1))/delta) + 1   !I found the j-level of the tested point
     v_test_LL(k,p) = ((x_test_LL(p) - x_mesh_LL(1))/(x_mesh_LL(MMM_LL) - x_mesh_LL(1))) * &
&                        (v_max_li - v_min_li) + v_min_li
!     v_test_LL_log(k,p) = ((x_test_LL(p) - x_mesh_LL(1))/(x_mesh_LL(MMM_LL) - x_mesh_LL(1))) * &
!&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
!     v_test_LL(k,p) = 10d0 ** v_test_LL_log(k,p)

!
!  CO2 BI-LINEAR LOOK-UP TABLES

      CALL CO2BLLT_EQUI(p_BL_LL(k,p), T_BL_LL(k,p), c_BL_LL(k,p), &
     &     x_out_LL(k,p), a_out, res, y_test_LL(k), v_test_LL(k,p),flag_loca)
!
       T_guess = T_BL_LL(k,p)
!
!      Values for comparison
      CALL New_Rap1D(1, T_LOIET_LL(k,p), out_2, res, N,&
      & flag, y_test_LL(k), T_guess, v_test_LL(k,p), out3)
      
      CALL New_Rap1D(1, T_LOIET_LL(k,p), out_2, res, N,&
      & flag, y_test_LL(k), T_LOIET_LL(k,p), v_test_LL(k,p), out3)

!                
!
      CALL pressure(T_LOIET_LL(k,p),v_test_LL(k,p),p_LOIET_LL(k,p))
      CALL sound_speed(T_LOIET_LL(k,p),v_test_LL(k,p), c_LOIET_LL(k,p))
!
      CALL heat_cap_p( T_LOIET_LL(k,p),v_test_LL(k,p),cp_LOIET_LL(k,p) ) 
      CALL CO2cp(cp_BL_LL(k,p),v_test_LL(k,p),y_test_LL(k), 0.0_pr, flag_loca)

!
!      print*, cp_LOIET_LL(k,p),cp_BL_LL(k,p),k,p , flag_loca
      IF (res > 1d-10) THEN
         print*, "1D_res_LL", res, "iter", N, k, p
         STOP
      ENDIF
!                
! For Bilinear
         err_pBL_LL(k,p)   = abs(p_BL_LL(k,p) - p_LOIET_LL(k,p))/abs(p_LOIET_LL(k,p)) 
         err_TBL_LL(k,p)   = abs(T_BL_LL(k,p) - T_LOIET_LL(k,p))/abs(T_LOIET_LL(k,p))
         err_cBL_LL(k,p)   = abs(c_BL_LL(k,p) - c_LOIET_LL(k,p))/abs(c_LOIET_LL(k,p))
         err_cp            = abs(cp_BL_LL(k,p) - cp_LOIET_LL(k,p))/abs(cp_LOIET_LL(k,p))
!
!
!       CALL inter_energy(T_cr+0.0001, 1_pr/rho_cr, energyc)
!       print*, 'internal energy critical', energyc, T_cr, 1_pr/rho_cr
!STOP
!
!      IF (err_cBL_LL(k,p) >0.01) THEN
!          print*, 'k',k,'p',p,T_BL_LL(k,p) 
!          print*, 'c table ',c_BL_LL(k,p),'c SW', c_LOIET_LL(k,p),&
!&                 'error', abs(c_BL_LL(k,p) - c_LOIET_LL(k,p))/abs(c_LOIET_LL(k,p))
!      ENDIF
      IF (err_pBL_LL(k,p) > 0.03_pr) THEN
         print*, "press_error" , err_pBL_LL(k,p),"row", k,"col", p,res
         print*, "p_interp",p_BL_LL(k,p),T_BL_LL(k,p),v_test_LL(k,p),y_test_LL(k)
         CALL inter_energy(T_BL_LL(k,p), v_test_LL(k,p), energy)
         print*, 'internal energy BL', energy
         CALL inter_energy(T_LOIET_LL(k,p) , v_test_LL(k,p), energy)
         print*, 'internal energy SW', energy
         CALL inter_energy(T_tri, 1_pr/rho_tri_L, energy)
         print*, 'internal energy triple liquid', energy
         print*, 'p_eos',p_LOIET_LL(k,p),T_LOIET_LL(k,p)
!         STOP
      ENDIF
      IF (err_cp >0.1) THEN
          print*, k,p,'T',T_BL_LL(k,p), 'error', err_cp
          print*, 'Table ',cp_BL_LL(k,p),'SW', cp_LOIET_LL(k,p)
      ENDIF
         WRITE(41,*) y_test_LL(k), v_test_LL(k,p), err_pBL_LL(k,p)
         WRITE(57,*) y_test_LL(k), v_test_LL(k,p), err_TBL_LL(k,p), err_TBL_LL(k,p)*abs(T_LOIET_LL(k,p))
         WRITE(63,*) y_test_LL(k), v_test_LL(k,p), err_cBL_LL(k,p), err_cp, cp_BL_LL(k,p), cp_LOIET_LL(k,p)
!         WRITE(64,*) y_test_LL(k), v_test_LL(k,p), p_BL_LL(k,p) 
!         WRITE(65,*) y_test_LL(k), v_test_LL(k,p), T_BL_LL(k,p) 
!               
!

   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
   CLOSE(UNIT = 63)
!   CLOSE(UNIT = 64)
!   CLOSE(UNIT = 65)
   
!
!----------------------------------------------------------------------
print*, 'LL Test finished'
!----------------------------------------------------------------------
STOP
CASE(2)
!======================================================
! Left High
!======================================================
!
y_test_LH(1)           = y_mesh_LH(1) + dtest_u
y_test_LH(NNN_test_LH) = y_mesh_LH(NNN_LH) - dtest_u
d_y_test_LH            = (y_test_LH(NNN_test_LH) - y_test_LH(1))/(NNN_test_LH - 1)

y_test_LH              = y_test_LH(1) + (/(j*d_y_test_LH, j=0,NNN_test_LH-1)/)

x_test_LH(1)           = x_mesh_LH(1) + dtest_v
x_test_LH(MMM_test_LH) = x_mesh_LH(MMM_LH) - dtest_v
d_x_test_LH            = (x_test_LH(MMM_test_LH) - x_test_LH(1))/(MMM_test_LH - 1)
      
x_test_LH              = x_test_LH(1) + (/(j*d_x_test_LH, j=0,MMM_test_LH-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'err_pBL_LH.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'err_TBL_LH.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'press_tables_LH.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 65,             FILE = 'Temp_tables_LH.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
OPEN (UNIT = 63,             FILE = 'err_cBL_LH.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_LH-1  !For each cell of the original grid to test
   delta = y_mesh_LH(21)  - y_mesh_LH(20) 
   i     = INT((y_test_LH(k) - y_mesh_LH(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li    = 0d0
   v_max_li    = 0d0
   v_max_log   = 0d0

   DO j = 1, ord_spline + 1
      v_min_li   = v_min_li   + spline_pmax_LH(ord_spline+2-j,i) * y_test_LH(k) ** (j - 1)
      v_max_li   = v_max_li   + spline_Lsat_LH(ord_spline+2-j,i) * y_test_LH(k) ** (j - 1)     
!      v_max_log  = v_max_log + spline_Lsat_LH(ord_spline+2-j,i) * y_test_LH(k) ** (j - 1)
!      v_max_li   = 10d0 ** v_max_log
!      print*, spline_pmax_LH(ord_spline+2-j,i)
   ENDDO

!   stop
   DO p = 1, MMM_test_LH-1

   delta = x_mesh_LH(2)  - x_mesh_LH(1)
     j = INT((x_test_LH(p) - x_mesh_LH(1))/delta) + 1   !I found the j-level of the tested point
!     v_test_LH_log(k,p) = ((x_test_LH(p) - x_mesh_LH(1))/(x_mesh_LH(MMM_LH) - x_mesh_LH(1))) * &
!&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
!     v_test_LH(k,p)     = 10d0 ** v_test_LH_log(k,p)
     v_test_LH(k,p) = ((x_test_LH(p) - x_mesh_LH(1))/(x_mesh_LH(MMM_LH) - x_mesh_LH(1))) * &
&                        (v_max_li - v_min_li) + v_min_li
     
     
     IF ( v_test_LH(k,p) .GE. v_max_li) THEN
        print*, 'Error in meshtest', k, p, y_test_LH(k), v_max_li, v_test_LH(k,p)!----->Problem in spline spin_LH
        stop
     ENDIF   
!        print*, ''
!   print*, v_min_li,v_test_LH(k,p), v_max_li
!       
!        print*, k, p, v_test_LH(k,p)
!  CO2 BI-LINEAR LOOK-UP TABLES
!     print*, y_test_LH(k),v_test_LH(k,p),p
        CALL CO2BLLT_EQUI (p_BL_LH(k,p), T_BL_LH(k,p), c_BL_LH(k,p), &
     &     x_out_LH(k,p), a_out, res, y_test_LH(k), v_test_LH(k,p),flag_loca)
!
        T_guess = T_BL_LH(k,p)
        p_guess = p_BL_LH(k,p) 

        IF (res > 1d-8) THEN
            print*, "middle_res_LH_Bi", res, "iter", N, 'row', k, 'column', p
            STOP
        ENDIF        

        CALL New_Rap1D(1, T_LOIET_LH(k,p), out_2, res, N,&
        & flag, y_test_LH(k), T_guess, v_test_LH(k,p), out3)  
 
!        CALL New_Rap1D(1, T_LOIET_LH(k,p), out_2, res, N,&
!        & flag, y_test_LH(k), T_LOIET_LH(k,p), v_test_LH(k,p), out3)



        CALL pressure(T_LOIET_LH(k,p),v_test_LH(k,p),p_LOIET_LH(k,p))
        CALL sound_speed(T_LOIET_LH(k,p),v_test_LH(k,p), c_LOIET_LH(k,p))
!
        CALL heat_cap_p( T_LOIET_LH(k,p),v_test_LH(k,p),cp_LOIET_LH(k,p) )
        CALL CO2cp(cp_BL_LH(k,p),v_test_LH(k,p),y_test_LH(k), 0.0_pr, flag_loca)
        IF (res > 1d-8) THEN
           print*, "middle_res_LH_1D", res, "iter", N, 'row', k, 'column', p
        !   print*,"p",p_LOIET_LH(k,p),"T",T_BL_LH(k,p)
           STOP
        ENDIF
        
!         print*, k,p,'p TB', p_BL_LH(k,p),'p SW', p_LOIET_LH(k,p),abs(p_BL_LH(k,p) - p_LOIET_LH(k,p))/1e6   
!         print*, k,p,'T TB', T_BL_LH(k,p),'T SW', T_LOIET_LH(k,p)        
!         print*, k,p,'c TB', c_BL_LH(k,p),'c SW', c_LOIET_LH(k,p)        
        
! For Bilinear
         err_pBL_LH(k,p)   = abs(p_BL_LH(k,p) - p_LOIET_LH(k,p))/abs(p_LOIET_LH(k,p)) 
         err_TBL_LH(k,p)   = abs(T_BL_LH(k,p) - T_LOIET_LH(k,p))/abs(T_LOIET_LH(k,p))
         err_cBL_LH(k,p)   = abs(c_BL_LH(k,p) - c_LOIET_LH(k,p))/abs(c_LOIET_LH(k,p))
         err_cp            = abs(cp_BL_LH(k,p) - cp_LOIET_LH(k,p))/abs(cp_LOIET_LH(k,p))


       IF (c_BL_LH(k,p) /= c_BL_LH(k,p)) THEN
       print*, 'y=',k,'x=',p
       print*, 'v',v_test_LH(k,p),'e',y_test_LH(k)
       print*, 'c', c_BL_LH(k,p), c_LOIET_LH(k,p)
       print*, 'TB', T_BL_LH(k,p), p_BL_LH(k,p)
       print*, 'SW', p_LOIET_LH(k,p), T_LOIET_LH(k,p)
       print*, ' '
       ENDIF
!
       IF (err_cp >0.1) THEN
          print*, k,p,'T',T_BL_LH(k,p), 'error', err_cp
          print*, 'Table ',cp_BL_LH(k,p),'SW', cp_LOIET_LH(k,p)
       ENDIF

         WRITE(41,*) y_test_LH(k), v_test_LH(k,p), err_pBL_LH(k,p)
         WRITE(57,*) y_test_LH(k), v_test_LH(k,p), err_TBL_LH(k,p),err_TBL_LH(k,p)*abs(T_LOIET_LH(k,p))
         WRITE(63,*) y_test_LH(k), v_test_LH(k,p), err_cBL_LH(k,p), err_cp, cp_BL_LH(k,p), cp_LOIET_LH(k,p)
      !   WRITE(64,*) y_test_LH(k), v_test_LH(k,p), p_BL_LH(k,p) 
       !  WRITE(65,*) y_test_LH(k), v_test_LH(k,p), T_BL_LH(k,p) 
                
!               
!

   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
   CLOSE(UNIT = 63)
!   CLOSE(UNIT = 64)
!   CLOSE(UNIT = 65)
!
!
!----------------------------------------------------------------------
print*, 'LH Test finished'
!----------------------------------------------------------------------
STOP
CASE(3)
!======================================================
! Right
!==========================================
!
y_test_R(1)           = y_mesh_R(1) + dtest_u
y_test_R(NNN_test_R)  = y_mesh_R(NNN_R) - dtest_u
d_y_test_R            = (y_test_R(NNN_test_R) - y_test_R(1))/(NNN_test_R - 1)
y_test_R              = y_test_R(1) + (/(j*d_y_test_R, j=0,NNN_test_R-1)/)
x_test_R(1)           = x_mesh_R(1) + dtest_v
x_test_R(MMM_test_R)  = x_mesh_R(MMM_R) - dtest_v
d_x_test_R            = (x_test_R(MMM_test_R) - x_test_R(1))/(MMM_test_R - 1)      
x_test_R              = x_test_R(1) + (/(j*d_x_test_R, j=0,MMM_test_R-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'err_pBL_R.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'err_TBL_R.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 58,             FILE = 'err_cBL_R.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'press_tables_R.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 65,             FILE = 'Temp_tables_R.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_R-1  !For each cell of the original grid to test
   delta = y_mesh_R(2)  - y_mesh_R(1) 
   i     = INT((y_test_R(k) - y_mesh_R(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li    = 0d0
   v_max_li    = 0d0

   DO j = 1, ord_spline + 1
      v_max_li   = v_max_li + spline_pmin(ord_spline+2-j,i) * y_test_R(k) ** (j - 1)
      v_min_li   = v_min_li + spline_Vsat(ord_spline+2-j,i) * y_test_R(k) ** (j - 1)  !I found the min and max volumes correspondent to that energy level

   ENDDO
!   print*, ''
!   print*, v_min_li, v_max_li
!   stop
   DO p = 1, MMM_test_R-1
   delta = x_mesh_R(2)  - x_mesh_R(1)
     j = INT((x_test_R(p) - x_mesh_R(1))/delta) + 1   !I found the j-level of the tested point
     v_test_R_log(k,p) = ((x_test_R(p) - x_mesh_R(1))/(x_mesh_R(MMM_R) - x_mesh_R(1))) * &
&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
     v_test_R(k,p)     = 10d0 ** v_test_R_log(k,p)
!     v_test_R(k,p) = ((x_test_R(p) - x_mesh_R(1))/(x_mesh_R(MMM_R) - x_mesh_R(1))) * &
!&                        (v_max_li - v_min_li) + v_min_li
     
     
     IF (( v_test_R(k,p) .LT. v_min_li) .OR. ( v_test_R(k,p) .GT. v_max_li))  THEN
        print*, 'Error in meshtest', k, p, y_test_R(k), v_min_li, v_test_R(k,p)!----->Problem in spline spin_R
        stop
     ENDIF       
      
      
     
!
!  CO2 BI-LINEAR LOOK-UP TABLES

        CALL CO2BLLT_EQUI (p_BL_R(k,p), T_BL_R(k,p), c_BL_R(k,p), &
     &     x_out_R(k,p), a_out, res, y_test_R(k), v_test_R(k,p),flag_loca)
!
        T_guess = T_BL_R(k,p)
        p_guess = p_BL_R(k,p)

        CALL New_Rap1D(1, T_LOIET_R(k,p), out_2, res, N,&
        & flag, y_test_R(k), T_guess, v_test_R(k,p), out3)
!                
!
        CALL pressure(T_LOIET_R(k,p),v_test_R(k,p),p_LOIET_R(k,p))
!
        CALL sound_speed(T_LOIET_R(k,p),v_test_R(k,p), c_LOIET_R(k,p))                                                                                     
!
        CALL heat_cap_p( T_LOIET_R(k,p),v_test_R(k,p),cp_LOIET_R(k,p) )
        CALL CO2cp(cp_BL_R(k,p),v_test_R(k,p),y_test_R(k), 0.0_pr, flag_loca)
!
!  
        IF (res > 1d-8) THEN
           print*, "middle_res_R", res, "iter", N
           STOP
        ENDIF
        
        
        
! For Bilinear
         err_pBL_R(k,p)   = abs(p_BL_R(k,p) - p_LOIET_R(k,p))/abs(p_LOIET_R(k,p)) 
         err_TBL_R(k,p)   = abs(T_BL_R(k,p) - T_LOIET_R(k,p))/abs(T_LOIET_R(k,p))
         err_cBL_R(k,p)   = abs(c_BL_R(k,p) - c_LOIET_R(k,p))/abs(c_LOIET_R(k,p))
         err_cp            = abs(cp_BL_R(k,p) - cp_LOIET_R(k,p))/abs(cp_LOIET_R(k,p))
!
         WRITE(41,*) y_test_R(k), v_test_R(k,p), err_pBL_R(k,p)
         WRITE(57,*) y_test_R(k), v_test_R(k,p), err_TBL_R(k,p), err_TBL_R(k,p)*abs(T_LOIET_R(k,p))
         WRITE(58,*) y_test_R(k), v_test_R(k,p), err_cBL_R(k,p), err_cp, cp_BL_R(k,p), cp_LOIET_R(k,p)
    !     WRITE(64,*) y_test_R(k), v_test_R(k,p), p_BL_R(k,p) 
     !    WRITE(65,*) y_test_R(k), v_test_R(k,p), T_BL_R(k,p) 
                
!               
!

   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
   CLOSE(UNIT = 58)
!   CLOSE(UNIT = 64)
!   CLOSE(UNIT = 65)
!
!
!----------------------------------------------------------------------
print*, 'R Test finished'
!----------------------------------------------------------------------
STOP
CASE(4)
!======================================================
! High Temperature
!==========================================
!
y_test_HT(1)           = y_mesh_HT(1) + dtest_u
y_test_HT(NNN_test_HT) = y_mesh_HT(NNN_HT) - dtest_u
d_y_test_HT            = (y_test_HT(NNN_test_HT) - y_test_HT(1))/(NNN_test_HT - 1)
y_test_HT              = y_test_HT(1) + (/(j*d_y_test_HT, j=0,NNN_test_HT-1)/)
x_test_HT(1)           = x_mesh_HT(1) + dtest_v
x_test_HT(MMM_test_HT)  = x_mesh_HT(MMM_HT) - dtest_v
d_x_test_HT            = (x_test_HT(MMM_test_HT) - x_test_HT(1))/(MMM_test_HT - 1)      
x_test_HT              = x_test_HT(1) + (/(j*d_x_test_HT, j=0,MMM_test_HT-1)/)
!
! I need a different mesh to test some points for each cell
!
!
OPEN (UNIT = 41,             FILE = 'err_pBL_HT.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 57,             FILE = 'err_TBL_HT.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
OPEN (UNIT = 58,             FILE = 'err_cBL_HT.txt', &
&     FORM = 'formatted',  ACTION = 'write',         &
&   STATUS = 'replace')
!
!OPEN (UNIT = 64,             FILE = 'press_tables_HT.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!OPEN (UNIT = 65,             FILE = 'Temp_tables_HT.txt', &
!&     FORM = 'formatted',  ACTION = 'write',         &
!&   STATUS = 'replace')
!
!
DO k = 1, NNN_test_HT-1  !For each cell of the original grid to test
   delta = y_mesh_HT(2)  - y_mesh_HT(1) 
   i     = INT((y_test_HT(k) - y_mesh_HT(1))/delta) + 1  !I found the i-level of the tested point
   v_min_li    = 0d0
   v_max_li    = 0d0

   DO j = 1, ord_spline + 1
      v_max_li   = v_max_li + spline_right_HT(ord_spline+2-j,i) * y_test_HT(k)**(j - 1)
      v_min_li   = v_min_li + spline_left_HT(ord_spline+2-j,i) * y_test_HT(k)**(j - 1)  !I found the min and max volumes correspondent to that energy level

   ENDDO
!   print*, ''
!   print*, v_min_li
!   stop
   DO p = 1, MMM_test_HT-1
   delta = x_mesh_HT(2)  - x_mesh_HT(1)
     j = INT((x_test_HT(p) - x_mesh_HT(1))/delta) + 1   !I found the j-level of the tested point
     v_test_HT_log(k,p) = ((x_test_HT(p) - x_mesh_HT(1))/(x_mesh_HT(MMM_HT) - x_mesh_HT(1))) * &
&                        (LOG10(v_max_li) - LOG10(v_min_li)) + LOG10(v_min_li)
     v_test_HT(k,p)     = 10d0 ** v_test_HT_log(k,p)
!     v_test_HT(k,p) = ((x_test_HT(p) - x_mesh_HT(1))/(x_mesh_HT(MMM_HT) - x_mesh_HT(1))) * &
!&                        (v_max_li - v_min_li) + v_min_li
     
     
     IF (( v_test_HT(k,p) .LT. v_min_li) .OR. ( v_test_HT(k,p) .GT. v_max_li))  THEN
        print*, 'Error in meshtest', k, p, y_test_HT(k), v_min_li, v_test_HT(k,p)!----->Problem in spline spin_HT
        stop
     ENDIF       
      
      
     
!
!  CO2 BI-LINEAR LOOK-UP TABLES

        CALL CO2BLLT_EQUI (p_BL_HT(k,p), T_BL_HT(k,p), c_BL_HT(k,p), &
     &     x_out_HT(k,p), a_out, res, y_test_HT(k), v_test_HT(k,p),flag_loca)
!
        T_guess = T_BL_HT(k,p)
        p_guess = p_BL_HT(k,p)
        
        CALL New_Rap1D(1, T_LOIET_HT(k,p), out_2, res, N,&
        & flag, y_test_HT(k), T_guess, v_test_HT(k,p), out3)
!                
!
        CALL pressure(T_LOIET_HT(k,p),v_test_HT(k,p),p_LOIET_HT(k,p))
!
        CALL sound_speed(T_LOIET_HT(k,p),v_test_HT(k,p), c_LOIET_HT(k,p))
                                                                                     
!  
        CALL heat_cap_p( T_LOIET_HT(k,p),v_test_HT(k,p),cp_LOIET_HT(k,p) )
        CALL CO2cp(cp_BL_HT(k,p),v_test_HT(k,p),y_test_HT(k), 0.0_pr, flag_loca)
!       
        IF (res > 1d-8) THEN
           print*, "middle_res_HT", res, "iter", N
           STOP
        ENDIF
        
! For Bilinear
         err_pBL_HT(k,p)   = abs(p_BL_HT(k,p) - p_LOIET_HT(k,p))/abs(p_LOIET_HT(k,p)) 
         err_TBL_HT(k,p)   = abs(T_BL_HT(k,p) - T_LOIET_HT(k,p))/abs(T_LOIET_HT(k,p))
         err_cBL_HT(k,p)   = abs(c_BL_HT(k,p) - c_LOIET_HT(k,p))/abs(c_LOIET_HT(k,p))
         err_cp            = abs(cp_BL_HT(k,p) - cp_LOIET_HT(k,p))/abs(cp_LOIET_HT(k,p))
!
         WRITE(41,*) y_test_HT(k), v_test_HT(k,p), err_pBL_HT(k,p)
         WRITE(57,*) y_test_HT(k), v_test_HT(k,p), err_TBL_HT(k,p), err_TBL_HT(k,p)*abs(T_LOIET_HT(k,p))
         WRITE(58,*) y_test_HT(k), v_test_HT(k,p), err_cBL_HT(k,p), err_cp, cp_BL_HT(k,p),cp_LOIET_HT(k,p)
                
!               
!

   ENDDO
ENDDO   !For all the tested values NNN_test x MMM_test in the Region
   CLOSE(UNIT = 41)
   CLOSE(UNIT = 57)
   CLOSE(UNIT = 58)
!   CLOSE(UNIT = 64)
!   CLOSE(UNIT = 65)
!
!
!
!----------------------------------------------------------------------
print*, 'HT Test finished'
!----------------------------------------------------------------------
STOP
END SELECT


END PROGRAM new
