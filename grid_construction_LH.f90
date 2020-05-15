!
! ===================================================================
!
!                      Grid_construction_Left_High
!
! ===================================================================
!
! Region Left High (LH): 
!
! From the Critical Point specific internal energy value up to the 
! specific internal energy value referred to the maximum point on
! the saturation curve. e_cr --> e_umax
! Bounded by the saturated vapour curve and the maximum pressure line
!
! ===================================================================
SUBROUTINE grid_construction_left_high()
!
USE def_constants
USE def_variables
USE non_linear_solvers
USE properties
USE grid_functions
!
INTEGER  ::  i, j, Niter, exitflag
!INTEGER, DIMENSION (NNN_LH,MMM_LH)  ::  
!
REAL(pr)   ::  utest_max, u_min, res, delta, vliq_guess_sat,&
 & T_guess_sat, p_guess, vliq_guess, out_2,v_guess_sat,in_2,&
 & press, sound, temp, v_guess2, T_guess2, delta_LH, u_six, v_six
!
!REAL(pr), DIMENSION (NNN_LH)     :: T_sat, p_sat, T_pmax, v_max, v_min, c_vmax, c_vmin
REAL(pr), DIMENSION (NNN_LH)     :: vliq_LH
REAL(pr), DIMENSION (NNN_sat_LH) :: T_sat_sat,T_sat_pmax,res_sat_LH,res_pmax_LH,vliq_sat_LH
REAL(pr), DIMENSION (NNN_LH,MMM_LH) :: res_LH,T_guess,v_guess, vvv_log
!
! data used for bicubic interpolation
!
!REAL(pr) :: v1, v2, v3, v4, u1, u2, u3, u4
!REAL(pr) :: F1, F2, F3, F4, G1, G2, G3, G4
!REAL(pr) :: a1, a2, a3, a4, b1, b2, b3, b4
!REAL(pr) :: dF1v, dF2v, dF3v, dF4v, dF1u, dF2u, dF3u, dF4u, dF1vu, dF2vu, dF3vu, dF4vu
!REAL(pr) :: dF1vv, dF2vv, dF3vv, dF4vv, dF1uu, dF2uu, dF3uu, dF4uu
!REAL(pr) :: dF1x, dF2x, dF3x, dF4x, dF1y, dF2y, dF3y, dF4y, dF1xy, dF2xy, dF3xy, dF4xy
!REAL(pr), DIMENSION (4, 4) :: A_bicub, C_bicub
!
!--------------------------------------------------
!
        v_six        = 2.2530239912767237d-3
!
!        T_guess(:,1) = 317.221 
!        v_guess(:,1) = 2.3348d-3
!guess for 50MPa too sensible for c++ not for fortran.
        T_guess(:,1) = 375.87 
        v_guess(:,1) = 1.23186d-3
!        
        TTT_LH(1,MMM_LH) = T_cr
        T_guess(1,MMM_LH) = T_cr+2_pr
        vvv_LH(1,MMM_LH) = 1_pr / rho_cr
        v_guess(1,MMM_LH) = 1_pr / rho_cr + 1d-4    !Values of the critical point fixed
        ppp_LH(1,MMM_LH) = P_cr
        vliq_guess = 1_pr / rho_cr - 1.5d-4
        vliq_LH(1) = 1_pr / rho_cr 
        ccc_LH(1,MMM_LH) = c_cr
!
        Niter = 0
        exitflag=0
!--------------------------------------------------------------------------
!
! Construct vector y_mesh_LH(physical domain), x_mesh_LH(tansformed
! domain)
!
!-------------------------------------------------------------------------
        utest_max = e_umax 
!        u_min     = e_cr+ (e_umax-e_cr)*0.05_pr
        
        delta_LH  = (e_umax - e_cr) / (NNN_LH - 1)
!        u_six     = 6d0 * delta_LH + e_cr
!
!        delta     = (utest_max - u_min) / (NNN_LH - 1)
        y_mesh_LH = e_cr + (/(i*delta_LH,i=0,NNN_LH-1)/)
!        y_mesh_LH(NNN_LH) = utest_max
! 
        delta     = (x_mesh_max - x_mesh_min) / (MMM_LH -1)
        x_mesh_LH = x_mesh_min + (/(i*delta, i=0,MMM_LH-1)/)
!        x_mesh_LH(NNN_LH) = x_mesh_max
!
! array y_mesh_sat_LH
!
        delta_LH          = (e_umax - e_cr)/(NNN_sat_LH-1)
        y_mesh_sat_LH     = e_cr + (/(i*delta_LH, i=0,NNN_sat_LH-1)/)
!
!
!--------------------------------------------------------------------------
!       
! Left boundary of LH domain. isobar Pmax,  and u_min-->u_max
!
!--------------------------------------------------------------------------
!
!P (iso) at each point of  left boundary
!
        ppp_LH(:,1) = p_max
!
!print*, 'start left BC'
        DO i=1,NNN_LH
!print*,i, y_mesh_LH(i),ppp_LH(i,1),T_guess(i,1), v_guess(i,1)
!
!        
              CALL New_Rap2D(2, TTT_LH(i,1), vvv_LH(i,1), &
              & res_LH(i,1), Niter, exitflag, ppp_LH(i,1), y_mesh_LH(i),&
              & T_guess(i,1), v_guess(i,1))
!
             CALL sound_speed(TTT_LH(i,1),vvv_LH(i,1), sound)
             ccc_LH(i,1) = sound
!
!print*,'c',ccc_LH(i,1)
              IF (res_LH(i,1) > res_ref) THEN
                print*, "left_res IN LH", res_LH(i,1), "iter", Niter,"flag",exitflag
                STOP
              ENDIF
              T_guess(i+1,1) = TTT_LH(i,1)
              v_guess(i+1,1) = vvv_LH(i,1)
!                
        END DO
!---------------------------------------------------------------------------
!
! Right boundary of LH domain. Saturation line u_min-->u_max **************
!
!--------------------------------------------------------------------------
!print*, 'start right BC'
        DO i = 2, 6
!        ! Calculation of the volumes in the linear part
           vvv_LH(i,MMM_LH) = (y_mesh_LH(i) - e_cr) / ( y_mesh_LH(7) - e_cr) * &
          & (v_six - 1d0 / rho_cr) + 1d0 / rho_cr
!       
           CALL New_Rap1D(1, TTT_LH(i,MMM_LH), out_2, res_LH(i,MMM_LH), Niter, &
           & exitflag, y_mesh_LH(i), TTT_LH(i-1,MMM_LH),&
           & vvv_LH(i,MMM_LH), out_2)
!           
           CALL pressure(TTT_LH(i,MMM_LH),vvv_LH(i,MMM_LH),press)
              ppp_LH(i,MMM_LH) = press
!!
           CALL sound_speed(TTT_LH(i,MMM_LH),vvv_LH(i,MMM_LH), sound)
              ccc_LH(i,MMM_LH) = sound
!        
           IF (res_LH(i,MMM_LH) > res_ref) THEN
                print*, "right_res IN LH", res_LH(i,MMM_LH),i
                STOP
           ENDIF
!           print*,i, 'press',press, 'c',sound,'v', vvv_LH(i,MMM_LH)     
!           print*, ' '
!           IF ( i == 5) THEN
!              T_guess(i,MMM_LH) = T_cr-0.5_pr
!              v_guess(i,MMM_LH) = 1_pr / rho_cr + 1e-4_pr
!              vliq_guess = 1_pr / rho_cr - 1e-4_pr
!           
!           ENDIF
!       
        ENDDO 
!
!        DO i=6, NNN_LH-1
!print*, 'start right BC 2'
        DO i=2, NNN_LH-1
!        
           IF (i<=6) THEN
              T_guess(i,MMM_LH) = TTT_LH(i,MMM_LH)
              v_guess(i,MMM_LH) = vvv_LH(i,MMM_LH) 
           ELSE
              T_guess(i,MMM_LH) = TTT_LH(i-1,MMM_LH)
              v_guess(i,MMM_LH) = vvv_LH(i-1,MMM_LH)
           ENDIF
!
              CALL New_Rap3D(4,vliq_LH(i),vvv_LH(i,MMM_LH),TTT_LH(i,MMM_LH),&
              & res_LH(i,MMM_LH), Niter, exitflag, y_mesh_LH(i),in_2,&
              & vliq_guess,v_guess(i,MMM_LH),T_guess(i,MMM_LH))
!!
              CALL pressure(TTT_LH(i,MMM_LH),vvv_LH(i,MMM_LH),press)
              ppp_LH(i,MMM_LH) = press
!!
              CALL sound_speed(TTT_LH(i,MMM_LH),vvv_LH(i,MMM_LH), sound)
              ccc_LH(i,MMM_LH) = sound
!!
              IF (res_LH(i,MMM_LH) > 1d-7) THEN
                print*, "right_res IN LH", res_LH(i,MMM_LH),i
                STOP
              ENDIF
!              T_guess(i,MMM_LH) = TTT_LH(i,MMM_LH)
!              v_guess(i,MMM_LH) = vvv_LH(i,MMM_LH)
              vliq_guess = vliq_LH(i) - 2.0d-4
!        print*, i,'press',press,'c',sound,'v', vvv_LH(i,MMM_LH),y_mesh_LH(i)
        END DO
!        print*, vvv_LH(NNN_LH-1,MMM_LH) 
        vvv_LH(NNN_LH,MMM_LH) = v_umax
!        vvv_LH(NNN_LH,MMM_LH) = 1.5e-2_pr
         
!        CALL New_Rap1D(1, temp, out_2, res_LH(NNN_LH,MMM_LH), Niter,&
!&                      exitflag, y_mesh_LH(NNN_LH), T_umax,&
!&                      vvv_LH(NNN_LH,MMM_LH), out_2)
!       print*, vvv_LH(NNN_LH,MMM_LH), vvv_LH(NNN_LH-1,MMM_LH)
        TTT_LH(NNN_LH,MMM_LH) = T_umax
!        print*, res_LH(NNN_LH,MMM_LH), temp
!        TTT_LH(NNN_LH,MMM_LH) = temp
        CALL pressure(TTT_LH(NNN_LH,MMM_LH),vvv_LH(NNN_LH,MMM_LH),press)
              ppp_LH(NNN_LH,MMM_LH) = press
!!
        CALL sound_speed(TTT_LH(NNN_LH,MMM_LH),vvv_LH(NNN_LH,MMM_LH), sound)
              ccc_LH(NNN_LH,MMM_LH) = sound
              
!        print*, 'Umax ',' press',press,'sound',sound,'v', vvv_LH(NNN_LH,MMM_LH),

        
!        DO i = NNN_LH,1,-1   ! 1, NNN_LH-1
!             
!              IF (i == 1) THEN
!                 v_guess2   = 1_pr / rho_cr + 1e-4_pr
!                 T_guess2   = T_cr - 1_pr
!                 vliq_guess = 1_pr / rho_cr - 1e-4_pr
!              ELSEIF(i == NNN_LH) THEN
!                 v_guess2   = v_umax
!                 T_guess2   = T_umax
!                 vliq_guess = 9.647d-4 !vvv_LL(INT(NNN_LL/2_pr),MMM_LL)
!              ELSEIF(i .GT. NNN_LH - 20) THEN
!                 v_guess2   = vvv_LH(i+1,MMM_LH) !/ 3d0
!                 T_guess2   = TTT_LH(i+1,MMM_LH) !+ 1d0
!                 vliq_guess = vliq_LH(i+1)  
!              ELSE
!                 v_guess2   = vvv_LH(i+1,MMM_LH) !/ 1.2d0
!                 T_guess2   = TTT_LH(i+1,MMM_LH) !+ 1d0
!                 vliq_guess = vliq_LH(i+1)
!              ENDIF
              
!              IF (i == 1) THEN
!                 v_guess2   = 1_pr / rho_cr + 1e-4_pr
!                 T_guess2   = T_cr - 1_pr
!                 vliq_guess = 1_pr / rho_cr - 1e-4_pr 
!              ELSE
!                 v_guess2   = vvv_LH(i-1,MMM_LH) !* 1.2d0
!                 T_guess2   = TTT_LH(i-1,MMM_LH) - 1d0
!                 vliq_guess = vliq_LH(i-1)
!              ENDIF
               
!               print*, 'guess i', i, vliq_guess,v_guess2,T_guess2
!                

!             CALL New_Rap3D(4,vliq_LH(i),vvv_LH(i,MMM_LH),TTT_LH(i,MMM_LH),&
!              & res_LH(i,MMM_LH), Niter, exitflag, y_mesh_LH(i),in_2,&
!              & vliq_guess,v_guess2,T_guess2)
              
!             print*, i, 'vliq',vliq_LH(i),'vvap',vvv_LH(i,MMM_LH),'T',TTT_LH(i,MMM_LH)
!
!             IF ((res_LH(i,MMM_LH) > 1e-3) .OR. (vvv_LH(i,MMM_LH) <0.0)) THEN
!                print*, "right_boundary IN LH", res_LH(i,MMM_LH), "iter", Niter
!                STOP
!              ENDIF
!
!              CALL pressure(TTT_LH(i,MMM_LH),vvv_LH(i,MMM_LH),press)
!              ppp_LH(i,MMM_LH) = press
!!!
!              CALL sound_speed(TTT_LH(i,MMM_LH),vvv_LH(i,MMM_LH), sound)
!              ccc_LH(i,MMM_LH) = sound
!
!             ENDDO              
!
         
!              IF (i .LT. 76) THEN
!                 T_guess(i,MMM_LH) = TTT_LH(i,MMM_LH)
!                 v_guess(i,MMM_LH) = vvv_LH(i,MMM_LH)
!                 vliq_guess = vliq_LH(i)
!              ELSE
!                 T_guess(i,MMM_LH) = 0.5d0*(TTT_LH(i,MMM_LH) + T_umax)
!                 v_guess(i,MMM_LH) = vvv_LH(i,MMM_LH)
!                 v_guess(i,MMM_LH) = 0.5d0*(vvv_LH(i,MMM_LH) + v_umax)
!                 vliq_guess = 0.5d0*(vliq_LH(i) + 1d0 / rho_cr) 
!              ENDIF
!
!
!print*, vvv_LH(1,MMM_LH), y_mesh_LH(1) , 'LH'
!---------------------------------------------------------------------------
!
! Middle domain of the two boundaries. v_min-->v_max, u_min-->u_max
! vvv_LH in the physical space is built by the X in the transformed
! space with a LINEAR SCALING TRANSFORMATION 
!
!--------------------------------------------------------------------------         
!print*, 'start middle'
        DO j=2,MMM_LH-1
           DO i=1, NNN_LH
!

                vvv_log(i,j) = ((x_mesh_LH(j) - x_mesh_min) /&
                & (x_mesh_max - x_mesh_min)) *(LOG10(vvv_LH(i,MMM_LH))-LOG10(vvv_LH(i,1)))+&
                & LOG10(vvv_LH(i,1))
                vvv_LH(i,j) = 10d0 ** vvv_log(i,j)
!                vvv_LH(i,j) = ((x_mesh_LH(j) - x_mesh_min)/(x_mesh_max - x_mesh_min)) *&
!                               (vvv_LH(i,MMM_LH)-vvv_LH(i,1)) + vvv_LH(i,1)
!                
                CALL New_Rap1D(1, temp, out_2, res_LH(i,j), Niter,&
               & exitflag, y_mesh_LH(i), TTT_LH(i,j-1),&
               & vvv_LH(i,j), out_2)
                CALL New_Rap1D(1, temp, out_2, res_LH(i,j), Niter,&
               & exitflag, y_mesh_LH(i), temp,&
               & vvv_LH(i,j), out_2)
!                 
                TTT_LH(i,j) = temp
!                
!
              CALL pressure(TTT_LH(i,j),vvv_LH(i,j),press)
              ppp_LH(i,j) = press
!
             CALL sound_speed(TTT_LH(i,j),vvv_LH(i,j), sound)
             ccc_LH(i,j) = sound
!
                IF (res_LH(i,j) > res_ref) THEN
                  print*, "middle_res_LH", res_LH(i,j), "iter", Niter, 'row', i, 'column', j
                  STOP
                ENDIF
!           
           END DO
        END DO
!        print*, 'LH ','2,max ',  'e ', y_mesh_LH(2), vvv_LH(2,MMM_LH),&
!&                                 TTT_LH(2,MMM_LH),ppp_LH(2,MMM_LH),ccc_LH(2,MMM_LH)  
!        print*, 'LH ','1,max-1 ','e ', y_mesh_LH(1), vvv_LH(1,MMM_LH-1),&
!&                                 TTT_LH(1,MMM_LH-1),ppp_LH(1,MMM_LH-1),ccc_LH(1,MMM_LH-1)  
!        print*, 'LH ','1,max ','e ', y_mesh_LH(1), vvv_LH(1,MMM_LH),&
!&                                 TTT_LH(1,MMM_LH),ppp_LH(1,MMM_LH),ccc_LH(1,MMM_LH)  
!        print*, 'LH ','2,max-1 ','e ', y_mesh_LH(2), vvv_LH(2,MMM_LH-1),&
!&                                 TTT_LH(2,MMM_LH-1),ppp_LH(2,MMM_LH-1),ccc_LH(2,MMM_LH-1)  
!        print*, vvv_LH(NNN_LH,MMM_LH), vvv_LH(NNN_LH,MMM_LH-1),vvv_LH(NNN_LH,MMM_LH-14)
!
!
! grid for Cp
!
!        DO j = 1,MMM_LH
!           DO  i = 1, NNN_LH
!               CALL heat_cap_p(TTT_LH(i,j),vvv_LH(i,j),cpcp_LH(i,j))
!!               
!               IF ( (cpcp_LH(i,j) /= cpcp_LH(i,j)) .OR. (cpcp_LH(i,j) <= 0.0) .OR.&
!&                   (cpcp_LH(i,j) >  3000e3_pr  ) ) THEN
!!                  print*,'cp in LH (kJ/(kgK))', cpcp_LH(i,j)/1000.0,i,j
!                  cpcp_LH(i,j) = 3000e3_pr
!               ENDIF
!!
!             ENDDO
!        ENDDO
!!
!
!------------------------------------------------------------
! 
!          MESH CONSTRUCTION for spline coefficients
!
!------------------------------------------------------------
!print*, 'start spline BC'
!
! SATURATION
!
        T_sat_sat(1) = T_cr
        v_Lsat_LH(1) = 1_pr/rho_cr
        vliq_sat_LH(1)= 1_pr/rho_cr
!        T_guess_sat = T_cr - 0.1_pr
!        v_guess_sat = v_Lsat_LH(1) + 1d-5
!        vliq_guess_sat = vliq_sat_LH(1) - 1d-5
        T_guess_sat = 304.12813
        v_guess_sat = 2.144818d-3
        vliq_guess_sat = 2.13175969d-3
!
        DO i=2, NNN_sat_LH-1
!              T_guess_sat = T_sat_sat(i-1)
!              v_guess_sat = v_Lsat_LH(i-1)
!              vliq_guess_sat = vliq_sat_LH(i-1)
!           IF( i == 2) THEN
!print*, 'in',y_mesh_sat_LH(i),in_2,vliq_guess_sat,v_guess_sat, T_guess_sat

                CALL New_Rap3D(4,vliq_sat_LH(i),v_Lsat_LH(i),T_sat_sat(i),&
              & res_sat_LH(i), Niter, exitflag, y_mesh_sat_LH(i),in_2,&
              & vliq_guess_sat,v_guess_sat, T_guess_sat)

!print*, 'out',vliq_sat_LH(i),v_Lsat_LH(i),T_sat_sat(i),res_sat_LH(i), Niter, exitflag
!              
!                CALL New_Rap3D(4,vliq_sat_LH(i),v_Lsat_LH(i),T_sat_sat(i),&
!              & res_sat_LH(i), Niter, exitflag, y_mesh_sat_LH(i),in_2,&
!              & vliq_sat_LH(i),v_Lsat_LH(i),T_sat_sat(i))
!
!           ELSE
! 
!               CALL New_Rap3D(4,vliq_sat_LH(i),v_Lsat_LH(i),T_sat_sat(i),&
!              & res_sat_LH(i), Niter, exitflag, y_mesh_sat_LH(i),in_2,&
!              & vliq_guess_sat,v_guess_sat, T_guess_sat)
!           ENDIF

              T_guess_sat = T_sat_sat(i)
              v_guess_sat = v_Lsat_LH(i)
              vliq_guess_sat = vliq_sat_LH(i)
!
!
              IF (res_sat_LH(i) > res_ref) THEN
                print*, "right_sat_res", res_sat_LH(i), "i", i,&
               &        'T=',T_sat_sat(i)
                print*, 'IN LH CONSTRUCTION'
                STOP
              ENDIF
!
        END DO
!        v_Lsat_LH(NNN_sat_LH) = 1.5e-2_pr
        v_Lsat_LH(NNN_sat_LH) = v_umax
        
!
!---------------------------------------------------------------------------

!
!print*, 'start spline2 BC'
! p_max
!
!
!        T_guess(1,1) = 316.44
!        v_guess(1,1) = 1.8266d-3
        T_guess_sat = 375.87
        v_guess_sat = 1.23186d-3
!        
        DO i=1,NNN_sat_LH
!        
              CALL New_Rap2D(2, T_sat_pmax(i), v_Lpmax_LH(i), &
              & res_pmax_LH(i), Niter, exitflag, p_max, y_mesh_sat_LH(i),&
              & T_guess_sat, v_guess_sat)
!
              IF (res_pmax_LH(i) > res_ref) THEN
                print*, "left_pmax_res IN LH", res_pmax_LH(i), "iter", Niter
                STOP
              ENDIF
              T_guess_sat = T_sat_pmax(i)
              v_guess_sat = v_Lpmax_LH(i)
!                
         END DO

! construction of spline coefficients associated to saturation
! curve nodes(What I need), meta curve nodes and pmax curve nodes(What I
! need)
!
DO i = 1, NNN_LH-1
   j = 3 * (i-1) + 1
!   CALL polyfit(spline_Lsat_LH(:,i), y_mesh_sat_LH(j:j+3),LOG10(v_Lsat_LH(j:j+3)),ord_spline)
   CALL polyfit(spline_Lsat_LH(:,i), y_mesh_sat_LH(j:j+3),v_Lsat_LH(j:j+3),ord_spline)
   CALL polyfit(spline_pmax_LH(:,i), y_mesh_sat_LH(j:j+3),v_Lpmax_LH(j:j+3),ord_spline)
ENDDO
!print*,'------------------------------------------------------'
print*,           'LH GRID FINISH'
!print*,'------------------------------------------------------'
!
!       
!
END SUBROUTINE grid_construction_left_high
