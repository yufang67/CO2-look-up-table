PROGRAM test_grid
!
      USE def_constants
      USE def_variables 
      USE Grid
      USE Interp_table 
!      USE non_linear_solvers
!      USE properties
      IMPLICIT NONE
!
      INTEGER :: flag!, Niter, exitflag,i,j
!      REAL(pr):: T, v, p, e, cv, cp, s, c, v_l, v_v, v_l_Span,v_v_Span,&
!&           resnorm, guess_1, guess_2, T_span, v_span, e_span,helmho1
!      REAL(pr):: p2, e2, cv2, cp2, s2, c2, helmho2,Helmholtz
      REAL(pr):: p_out, T_out, c_out, x_out, a_out,res, u_in,v_in
!
!
!
!
        CALL MAKE_GRID()
!STOP
        u_in = -178470
        v_in= 0.00213858
!        u_in = -207343.3160311
!        u_in = -108327.291813182
!        v_in = 1.8626809245968707e-003
!        v_in = 1.582810157285395E-002
        
        CALL CO2BLLT_EQUI(p_out, T_out, c_out, x_out, a_out, res, u_in,v_in,flag) 

        print*, 'p=', p_out, 'T=',T_out, 'c=',c_out,'x=',x_out,'a=', a_out,flag

!
!================================================================================
!
!               Test creation of grid in LL
!
!================================================================================
!
!print*, "In Physical domaine"
!  Do j = 1,MMM_LL
!    Do i = 1, NNN_LL
!        print*, "Number of point", j,i
!        print*, "specific volume", vvv_LL(i,j)
!        print*, "pressure", ppp_LL(i,j)        
!        print*, "temperature", TTT_LL(i,j)
!        print*, "sound speed", ccc_LL(i,j)
!        print*, " "
!    ENDDO
!  ENDDO
!print*, "======================================================================"
!print*, "Transform domain"
!  DO i = 1,NNN_LL
!        print*, "X", x_mesh_LL(i)
!        print*, "Y", y_mesh_LL(i)
!  ENDDO
!print*,"========================================================================"
!print*, "Spline line"
!  DO i = 1,NNN_sat_LL
!        print*, "Y", y_mesh_sat_LL(i)
!        print*, "v liquid", v_Lsat_LL(i)
!        print*, "v pmax", v_Lpmax_LL(i)
!  ENDDO
!print*,"========================================================================"
!  DO i = 1,NNN_sat_LL
!        print*, "v meta", v_liq_meta(i),i 
!  ENDDO
!print*, "====================================================================="        
!print*, "spline coeff"
!DO i=1, ord_spline+1
!  DO j=1, NNN_LL-1
!       print*, "LL saturation line", spline_Lsat_LL(i,j)
!       print*, "LL Pressure line" , spline_pmax_LL(i,j)
!       print*, "LL meta" ,spline_meta_LL(i,j)
!       print*, " "
!  ENDDO
!ENDDO


!print*, "====================================================================="        
!
!================================================================================
!
!               Test creation of grid in LH
!
!================================================================================
!
!print*, "In Physical domaine"
!  Do j = 1,MMM_LH
!    Do i = 1, NNN_LH
!        print*, "Number of point", j,i
!        print*, "specific volume", vvv_LH(i,j)
!        print*, "pressure", ppp_LH(i,j)        
!        print*, "temperature", TTT_LH(i,j)
!        print*, "sound speed", ccc_LH(i,j)
!        print*, " "
!    ENDDO
!  ENDDO
!print*,"======================================================================"
!print*, "Transform domain"
!  DO i = 1,NNN_LH
!        print*, "X", x_mesh_LH(i)
!        print*, "Y", y_mesh_LH(i)
!  ENDDO
!print*,"========================================================================"
!print*, "Spline line"
!  DO i = 1,NNN_sat_LH
!        print*, i
!        print*, y_mesh_sat_LH(i), v_Lsat_LH(i),i
!        print*, "v vapor", v_Lsat_LH(i),i
!        print*, "v pmax", v_Lpmax_LH(i)
!        print*, " " 
! ENDDO
!DO i = 1,NNN_LH
!        print*, y_mesh_LH(i),vvv_LH(i,MMM_LH),i
!ENDDO


END PROGRAM test_grid
