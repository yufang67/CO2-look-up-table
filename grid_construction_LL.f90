!============================================================================
      SUBROUTINE grid_construction_left_low()
!============================================================================
! Region Left Low (LL)
! 
! The specific internal energy along from Triple point to Critical
! point. The liquid region is bounded by the saturation line and the maximum
! pressure line, the minimum internal energy and the critical internal
! energy. 
!
!  ** left isobar boundary (e,P)-->(T,v) 2D New_Raph
!     right saturation boundary (e)-->(vl,vv,T)-->(P) 3D New_Raph
!     middle (v,e)-->(T)--(P) 1D New_Raph
! 
! The bilinear mapping is used for the interpolation from physical (v,e)
! doamin to the transformed domaine (X,Y).(Bicubic can be
! considered later by using also the derivatives)
!========================================================================
!
     USE def_constants
     USE def_variables
     USE non_linear_solvers
     USE properties
     USE grid_functions
!        
      IMPLICIT NONE
!
      INTEGER  ::  i, j, N, exitflag, i1, i2, Niter
!      INTEGER, DIMENSION (NNN_LL,MMM_LL) ::  flag_pT_LL, N_pT_LL
!
      REAL(pr)  ::  u_min, res, delta, p_guess, extern, utest_max,&
&                   T_guess_sat,v_guess_sat, vap_guess_sat
!      REAL(pr)  :: p0, v_inters, T_inters, u_inters
!
      REAL(pr), DIMENSION (NNN_LL)  ::  vap_guess, vap_LL
!      
      REAL(pr), DIMENSION (NNN_LL,MMM_LL)  ::  res_pT_LL, res_LL,&
      &                                        T_guess, v_guess,vvv_LL_log
      REAL(pr), DIMENSION (NNN_sat_LL) :: res_sat_LL,vap_sat_LL, T_sat_sat,&
&                                         T_sat_pmax,res_pmax_LL,T_liq_meta
      REAL(pr)  ::  energy,press, temp, out_2, sound,in_2,guess3,out3
!
! data used for bicubic interpolation
!
!REAL(8) :: v1, v2, v3, v4, u1, u2, u3, u4
!REAL(8) :: F1, F2, F3, F4, G1, G2, G3, G4
!REAL(8) :: a1, a2, a3, a4, b1, b2, b3, b4
!REAL(8) :: dF1v, dF2v, dF3v, dF4v, dF1u, dF2u, dF3u, dF4u, dF1vu, dF2vu,
!dF3vu, dF4vu
!REAL(8) :: dF1vv, dF2vv, dF3vv, dF4vv, dF1uu, dF2uu, dF3uu, dF4uu
!REAL(8) :: dF1x, dF2x, dF3x, dF4x, dF1y, dF2y, dF3y, dF4y, dF1xy, dF2xy,
!dF3xy, dF4xy
!
!REAL(8), DIMENSION (4, 4) :: A_bicub, C_bicub
!
!
! Initial guess for point (1,1) and (1,MMM_LL)
!
        T_guess(1,1) = T_tri 
        v_guess(1,1) = 1_pr / rho_tri_L
!        
        TTT_LL(1,MMM_LL) = T_tri
        T_guess(1,MMM_LL) = T_tri
!
        vvv_LL(1,MMM_LL) = 1_pr / rho_tri_L
        v_guess(1,MMM_LL) = 1_pr / rho_tri_L
        ppp_LL(1,MMM_LL) = P_tri
        vap_guess(1) = 1_pr / rho_tri_R
        vap_LL(1) = 1_pr / rho_tri_R
        ccc_LL(1,MMM_LL) = c_tri_L
!--------------------------------------------------------------------------
!
! Construct vector y_mesh_LL(physical domain), x_mesh_LL(tansformed
! domain)
!
!-------------------------------------------------------------------------
        utest_max = 1.0_pr*e_cr ! interal energy smaller than cirtical internal energy
        u_min = e_tri_L       
!
        delta = (utest_max - u_min) / (NNN_LL - 1)             
        y_mesh_LL = u_min + (/(i*delta,i=0,NNN_LL-1)/)
!        y_mesh_LL(NNN_LL) = utest_max
! 
        delta = (x_mesh_max - x_mesh_min) / (MMM_LL -1)
        x_mesh_LL = x_mesh_min + (/(i*delta, i=0,MMM_LL-1)/)
!        x_mesh_LL(NNN_LL) = x_mesh_max
!
! array y_mesh_sat_LL
!
        delta             = (utest_max - u_min)/(NNN_sat_LL-1)
        y_mesh_sat_LL     = u_min + (/(i*delta, i=0,NNN_sat_LL-1)/)
!
!
!--------------------------------------------------------------------------
!       
! Left boundary of LL domain. isobar Pmax,  and u_min-->u_max
!
!--------------------------------------------------------------------------
!
!P (iso) at each point of  left boundary
!
        ppp_LL(:,1) = p_max
!
        DO i=1,NNN_LL
!        
              CALL New_Rap2D(2, TTT_LL(i,1), vvv_LL(i,1), &
              & res_LL(i,1), Niter, exitflag, ppp_LL(i,1), y_mesh_LL(i),&
              & T_guess(i,1), v_guess(i,1))
!
             CALL sound_speed(TTT_LL(i,1),vvv_LL(i,1), sound)
             ccc_LL(i,1) = sound
!
!             CALL inter_energy(TTT_LL(i,1),vvv_LL(i,1),energy)
!        print*, i, abs(y_mesh_LL(i) - energy ) / y_mesh_LL(i)
              IF (res_LL(i,1) > res_ref) THEN
                print*, "left_res IN LL", res_LL(i,1), "iter", Niter,"flag",exitflag
                STOP
              ENDIF
              T_guess(i+1,1) = TTT_LL(i,1)
              v_guess(i+1,1) = vvv_LL(i,1)
!                
         END DO

!---------------------------------------------------------------------------
!
! Right boundary of LL domain. Saturation line u_min-->u_max
!
!--------------------------------------------------------------------------
        DO i=2, NNN_LL-1
              T_guess(i,MMM_LL) = TTT_LL(i-1,MMM_LL)
              v_guess(i,MMM_LL) = vvv_LL(i-1,MMM_LL)
              vap_guess(i) = vap_LL(i-1)
!         
              CALL New_Rap3D(1,vvv_LL(i,MMM_LL),vap_LL(i),TTT_LL(i,MMM_LL),&
              & res_LL(i,MMM_LL), Niter, exitflag, y_mesh_LL(i),in_2,&
              & v_guess(i,MMM_LL), vap_guess(i),T_guess(i,MMM_LL))
!
              CALL pressure(TTT_LL(i,MMM_LL),vvv_LL(i,MMM_LL),press)
              ppp_LL(i,MMM_LL) = press
!
              CALL sound_speed(TTT_LL(i,MMM_LL),vvv_LL(i,MMM_LL), sound)
              ccc_LL(i,MMM_LL) = sound
!
!        print*, "vap_v", vap_LL(i)
              IF (res_LL(i,MMM_LL) > res_ref) THEN
                print*, "right_res IN LL", res_LL(i,MMM_LL), "iter", Niter
                print*, 'i',i,'v',vvv_LL(i,MMM_LL),'vap',vap_LL(i),'e',y_mesh_LL(i)
                STOP
              ENDIF
!
!        print*, "v_right", vvv_LL(i,MMM_LL), 't_right',TTT_LL(i,MMM_LL)
!        print*, 'e', y_mesh_LL()                        
        END DO
! Last point is fixed to critical point values
        vvv_LL(NNN_LL,MMM_LL) = 1_pr/rho_cr
        TTT_LL(NNN_LL,MMM_LL) = T_cr        
        ppp_LL(NNN_LL,MMM_LL) = P_cr
        ccc_LL(NNN_LL,MMM_LL) = c_cr
!        print*, 'LL','n-1','e', y_mesh_LL(NNN_LL-1), vvv_LL(NNN_LL-1,MMM_LL),&
!&                               TTT_LL(NNN_LL-1,MMM_LL),ppp_LL(NNN_LL-1,MMM_LL),ccc_LL(NNN_LL-1,MMM_LL)                         
!        print*, 'LL','n',  'e', y_mesh_LL(NNN_LL), vvv_LL(NNN_LL,MMM_LL),&
!&                               TTT_LL(NNN_LL,MMM_LL),ppp_LL(NNN_LL,MMM_LL),ccc_LL(NNN_LL,MMM_LL)        

!
!              T_guess(NNN_LL-2,MMM_LL) = TTT_LL(NNN_LL-3,MMM_LL)
!              v_guess(NNN_LL-2,MMM_LL) = 1_pr/rho_cr
!              vap_guess(NNN_LL-2) = vap_LL(NNN_LL-3)
!        DO i=NNN_LL-2, NNN_LL
!         
!              CALL New_Rap3D(1,vvv_LL(i,MMM_LL),vap_LL(i),TTT_LL(i,MMM_LL),&
!              & res_LL(i,MMM_LL), Niter, exitflag, y_mesh_LL(i),&
!              & v_guess(i,MMM_LL), vap_guess(i),T_guess(i,MMM_LL))
!
!              CALL pressure(TTT_LL(i,MMM_LL),vvv_LL(i,MMM_LL),press)
!              ppp_LL(i,MMM_LL) = press
!
!             CALL sound_speed(TTT_LL(i,MMM_LL),vvv_LL(i,MMM_LL), sound)
!             ccc_LL(i,MMM_LL) = sound
!              T_guess(i+1,MMM_LL) = TTT_LL(i,MMM_LL)
!              v_guess(i+1,MMM_LL) = vvv_LL(i,MMM_LL)
!              vap_guess(i+1) = vap_LL(i)
!
!       print*, Niter, res_LL(i,MMM_LL)
!              IF (res_LL(i,MMM_LL) > res_ref) THEN
!                print*, "right_res", res_LL(i,MMM_LL), "iter", Niter
!                STOP
!              ENDIF
!
!        print*, "v_right", vvv_LL(i,MMM_LL), 't_right',TTT_LL(i,MMM_LL)
!        print*, 'e', y_mesh_LL(i)                        
!        END DO

!---------------------------------------------------------------------------
!
! Middle domain of the two boundaries. v_min-->v_max, u_min-->u_max
! vvv_LL in the physical space is built by the X in the transformed
! space with a LINEAR SCALING TRANSFORMATION 
!
!--------------------------------------------------------------------------         
        DO j=2,MMM_LL-1
           DO i=1, NNN_LL
!
                vvv_LL(i,j) = ((x_mesh_LL(j) - x_mesh_min) /&
                & (x_mesh_max - x_mesh_min)) *(vvv_LL(i,MMM_LL)-vvv_LL(i,1))+&
                & vvv_LL(i,1)
!                vvv_LL_log(i,j) = ((x_mesh_LL(j) - x_mesh_min) /&
!                & (x_mesh_max - x_mesh_min)) *(log10(vvv_LL(i,MMM_LL))-log10(vvv_LL(i,1)))+&
!                & log10(vvv_LL(i,1))
!               vvv_LL(i,j) = 10_pr**vvv_LL_log(i,j)
!                
                CALL New_Rap1D(1, temp, out_2, res_LL(i,j), Niter,&
               & exitflag, y_mesh_LL(i), TTT_LL(i,j-1),&
               & vvv_LL(i,j), out3)
!                
                CALL New_Rap1D(1, temp, out_2, res_LL(i,j), Niter,&
               & exitflag, y_mesh_LL(i), temp,&
               & vvv_LL(i,j), out3)

!                
                TTT_LL(i,j) = temp
!                
!
              CALL pressure(TTT_LL(i,j),vvv_LL(i,j),press)
              ppp_LL(i,j) = press
!
             CALL sound_speed(TTT_LL(i,j),vvv_LL(i,j), sound)
             ccc_LL(i,j) = sound
!
                IF (res_LL(i,j) > res_ref) THEN
                  print*, "middle_res_LL", res_LL(i,j), "iter", Niter,i,j
                  STOP
                ENDIF
!           
!                CALL inter_energy(temp, vvv_LL(i,j),energy)
!                print*, "error e", (y_mesh_LL(i) - energy)/y_mesh_LL(i)
!                print*, "residu", res_LL(i,j)   
           END DO
        END DO
!        print*, 'LL ','n-1 ','e ', y_mesh_LL(NNN_LL), vvv_LL(NNN_LL,MMM_LL-1),&
!&                                  TTT_LL(NNN_LL,MMM_LL-1),ppp_LL(NNN_LL,MMM_LL-1),ccc_LL(NNN_LL,MMM_LL-1)  
!
!
!print*, 'grid Cp'
! grid for Cp
!
!        DO j = 1,MMM_LL
!           DO  i = 1, NNN_LL
!               CALL heat_cap_p(TTT_LL(i,j),vvv_LL(i,j),cpcp_LL(i,j))
!               IF ( (cpcp_LL(i,j) /= cpcp_LL(i,j)) .OR. (cpcp_LL(i,j) <= 0.0) .OR.&
!&                   (cpcp_LL(i,j) >  cp_cr  ) ) THEN
!!                  print*,'cp in LL kJ/(kgK)', cpcp_LL(i,j)/1000.0,i,j
!                  cpcp_LL(i,j) = cp_cr
!                ENDIF
!             ENDDO
!        ENDDO        
!
!------------------------------------------------------------
! 
!          MESH CONSTRUCTION for spline coefficients
!
!------------------------------------------------------------
!print*, 'start spline'
!
! SATURATION
!
        T_sat_sat(1) = T_tri
        v_Lsat_LL(1) = 1_pr/rho_tri_L
        vap_sat_LL(1)= 1_pr/rho_tri_R
! 21/08/2017 fix last point to critical point
        T_sat_sat(NNN_sat_LL) = T_cr
        v_Lsat_LL(NNN_sat_LL) = 1_pr/rho_cr
        vap_sat_LL(NNN_sat_LL)= 1_pr/rho_cr
!
        T_sat_sat(NNN_sat_LL-1) = 304.1281
        v_Lsat_LL(NNN_sat_LL-1) = 1_pr/rho_cr - 2d-5 
        vap_sat_LL(NNN_sat_LL-1)= 1_pr/rho_cr + 2d-5
! 20/05/2018 fix (NNN_sat_LL-1) point value
        DO i=2, NNN_sat_LL-2
              T_guess_sat = T_sat_sat(i-1)+4d-6
              v_guess_sat = v_Lsat_LL(i-1)!+1d-6
              vap_guess_sat = vap_sat_LL(i-1)!-1d-6
!print*, i, T_guess_sat, v_guess_sat,vap_guess_sat
!         
              CALL New_Rap3D(1,v_Lsat_LL(i),vap_sat_LL(i),T_sat_sat(i),&
              & res_sat_LL(i), Niter, exitflag, y_mesh_sat_LL(i),in_2,&
              & v_guess_sat, vap_guess_sat,T_guess_sat)

!              CALL New_Rap3D(1,v_Lsat_LL(i),vap_sat_LL(i),T_sat_sat(i),&
!              & res_sat_LL(i), Niter, exitflag, y_mesh_sat_LL(i),in_2,&
!              & v_Lsat_LL(i), vap_sat_LL(i), T_sat_sat(i))
!
!print*, i, T_sat_sat(i), v_Lsat_LL(i), vap_sat_LL(i)
!
              IF (res_sat_LL(i) > res_ref) THEN
                print*, "right_sat_res", res_sat_LL(i), "iter", Niter,&
              &         'T=',T_sat_sat(i),i
                print*, 'Tguess',T_guess_sat,'vlquess',v_guess_sat,vap_guess_sat
                print*,'IN LL CONSTRUCTION'
                STOP
              ENDIF
!
        END DO
!
!
! META: compute spinodal liquid between  -290;004 - -190.3 (kj/kg)
!
!e = -427.18 - -290.004  (kj/kg) p = 0
!      T_guess_sat = 263.0
!      v_guess_sat = 1.1e-3_pr
! DO i = 1, 172
!    CALL New_Rap2D(4, T_liq_meta(i), v_liq_meta(i), &
!              & res, N, exitflag, y_mesh_sat_LL(i),in_2,&
!              & T_guess_sat,v_guess_sat)
!!
!   IF (res > 1e-5) THEN
!       print*,'LL_sat_spin-Convergence failed calculating T_sat_sat at spinodal(',i,')'
!       print*,'resnorm     = ',res ,'N_iter = ',N
!       STOP
!   ENDIF
!      T_guess_sat = T_liq_meta(i)
!      v_guess_sat = v_liq_meta(i)+0.1e-4_pr
!      CALL pressure(T_liq_meta(i),v_liq_meta(i), press)
!   print*, 'meta T_liqmeta', T_liq_meta(i), 'vliqmeta', v_liq_meta(i), y_mesh_sat_LL(i),i
!   print*,' press', press
!    T_liq_meta(i) = T_sat_sat(i)
!    v_liq_meta(i) = v_Lsat_LL(i)
! ENDDO
!STOP
!
!
!
!
!e = -290.004 - -201.476 (kj/kg) dp/dv_T = 0
!      T_guess_sat = 275.0
!      v_guess_sat = v_tri_L +1e-4_pr
! DO i = 173, 284
!!
!        CALL New_Rap2D(3, T_liq_meta(i), v_liq_meta(i), &
!              & res, N, exitflag, y_mesh_sat_LL(i),in_2,&
!              & T_guess_sat,v_guess_sat)
!!
!   IF (res > 1e-5) THEN
!       print*,'LL_sat_spin-Convergence failed calculating T_sat_sat at spinodal(',i,')'
!       print*,'resnorm     = ',res ,'N_iter = ',N
!       STOP
!   ENDIF
!   IF (i<200) THEN
!      T_guess_sat = T_liq_meta(i)
!      v_guess_sat = v_liq_meta(i)+1e-4_pr
!!   ELSEIF (i>=284) THEN
!!      T_guess_sat = T_liq_meta(i)+0.1_pr
!!      v_guess_sat = v_liq_meta(i)
!   ELSE 
!      T_guess_sat = T_liq_meta(i)
!      v_guess_sat = v_liq_meta(i)-1e-4_pr
   
!   ENDIF
!      CALL pressure(T_liq_meta(i),v_liq_meta(i), press)
!   print*, 'meta T_liqmeta', T_liq_meta(i), 'vliqmeta', v_liq_meta(i), y_mesh_sat_LL(i),i
!   print*,' press', press
!   print*, 'sat', v_Lsat_LL(i), T_sat_sat(i)
! ENDDO
!!e = -201.476 - -190.3 (kj/kg) saturation line
!   DO i = 285, NNN_sat_LL
!      T_liq_meta(i) = T_sat_sat(i)
!      v_liq_meta(i) = v_Lsat_LL(i)
! ENDDO
!!   v_liq_meta(NNN_sat_LL)  = 1_pr/rho_cr
!
!
!print*, 'start p_max spline'
! P_max
!
!
        T_guess(1,1) = T_tri
        v_guess(1,1) = 1_pr/rho_tri_L
!        
        DO i=1,NNN_sat_LL
!        
              CALL New_Rap2D(2, T_sat_pmax(i), v_Lpmax_LL(i), &
              & res_pmax_LL(i), Niter, exitflag, p_max, y_mesh_sat_LL(i),&
              & T_guess(i,1), v_guess(i,1))
!
              IF (res_pmax_LL(i) > res_ref) THEN
                print*, "left_pmax_res IN LL ", res_pmax_LL(i), "iter", Niter
                STOP
              ENDIF
              T_guess(i+1,1) = T_sat_pmax(i)
              v_guess(i+1,1) = v_Lpmax_LL(i)
!                
         END DO
         
!
!
! construction of spline coefficients associated to saturation
! curve nodes(What I need), meta curve nodes and pmax curve nodes(What I  need)
!
DO i = 1, NNN_LL-1
   j = ord_spline * (i-1) + 1
   CALL polyfit(spline_Lsat_LL(:,i), y_mesh_sat_LL(j:j+ord_spline),v_Lsat_LL(j:j+ord_spline), ord_spline)
!   CALL polyfit(spline_Lsat_LL(:,i), y_mesh_sat_LL(j:j+3),LOG10(v_Lsat_LL(j:j+3)), ord_spline)
!   CALL polyfit(spline_meta_LL(:,i), y_mesh_sat_LL(j:j+ord_spline),v_liq_meta(j:j+ord_spline), ord_spline)
   CALL polyfit(spline_pmax_LL(:,i), y_mesh_sat_LL(j:j+ord_spline),v_Lpmax_LL(j:j+ord_spline), ord_spline)
ENDDO
!print*,'------------------------------------------------------'
print*,           'LL GRID FINISH'
!print*,'------------------------------------------------------'
!
!
!
      END SUBROUTINE grid_construction_left_low

