!============================================================================
      SUBROUTINE grid_construction_TPM()
!============================================================================
! Region two-phase (TP)
! 
! The specific internal energy along from critical point to triple
! point. The region is bounded by the saturation line LH and the triple
! pressure line, the critical internal, point E energy. 
! 
! 
! The bilinear mapping is used for the interpolation from physical (v,e)
! doamin to the transformed domaine (X,Y).(Bicubic can be
! considered later by using also the derivatives)
!========================================================================
!
     USE def_constants, ONLY: pr,ord_spline,e_tri_R,e_cr,x_mesh_max,x_mesh_min,&
&                             e_tri_R,e_tri_L,v_tri_R,v_tri_L,P_tri,T_tri,c_cr,MMM_TPM,NNN_TPM,NNN_sat_TPM,&
&                             P_cr,T_cr,rho_cr
     USE def_variables, ONLY: x_mesh_TPM,y_mesh_TPM,y_mesh_sat_TPM,v_right_TPM,v_left_TPM,&
&        spline_right_TPM,spline_left_TPM,saturP,ppp_TPM,vvv_TPM,xxx_TPM,ccc_TPM,TTT_TPM
!&        vL_psat_spline, vV_psat_spline,uL_psat_spline, uV_psat_spline, Tsat_psat_spline
     USE non_linear_solvers,ONLY: New_Rap3D, BrentRoots
     USE properties, ONLY: pressure,sound_speed,satderiv,satprop
     USE grid_functions,ONLY: polyfit
!        
      IMPLICIT NONE
!     
      REAL(pr) :: u_max, u_min , delta, deltap, vvv_TPM_log
      REAL(pr) :: ev,el,vv,vl,qual,sound,pp,ratio,press,temp,Tsat
      INTEGER  :: i,k,j_tpl, Niter, flag,j
      REAL(pr) :: duL_dp, duV_dp, dvL_dp, dvV_dp,du_dp_x,dv_dp_x
      REAL(pr) :: pguess,res, T_guess,v_six,in_2,out_2,v_guess,vliq,vliq_guess,pguess1,&
&                 Tguess,vvguess,vlguess, res1,dummy,upper, lower




!--------------------------------------------------------------------------
!
! Construct vector y_mesh_TPM(physical domain), x_mesh_TPM(tansformed
! domain)
!
!-------------------------------------------------------------------------
        u_max = e_tri_R
        u_min = e_cr
!
        delta = (u_max - u_min) / (NNN_TPM - 1)
        y_mesh_TPM = u_min + (/(i*delta,i=0,NNN_TPM-1)/)
! 
        delta = (x_mesh_max - x_mesh_min) / (MMM_TPM -1)
        x_mesh_TPM = x_mesh_min + (/(i*delta, i=0,MMM_TPM-1)/)
!
        delta             = (u_max - u_min)/(NNN_sat_TPM-1)
        y_mesh_sat_TPM   = u_min + (/(i*delta, i=0,NNN_sat_TPM-1)/)
!
       
! Right boundary = iso triple pressure line
!print*, 'stat right'
        ev = e_tri_R
        el = e_tri_L
        vv = v_tri_R
        vl = v_tri_L
        ppp_TPM(:,MMM_TPM) = P_tri
        TTT_TPM(:,MMM_TPM) = T_tri
!
        DO k=1,NNN_TPM
!        
           xxx_TPM(k,MMM_TPM) = (y_mesh_TPM(k) - el) / (ev - el)
           qual = (y_mesh_TPM(k) - el) / (ev - el)
           vvv_TPM(k,MMM_TPM) = qual*vv + (1_pr-qual)*vl 
! Speed of sound  
!         
           CALL satderiv(3, ppp_TPM(k,MMM_TPM), duL_dp, duV_dp, dvL_dp, dvV_dp)
!
           ratio = ((ev - el)/(vv - vl))   ! (J/kg)/(m3/kg)
           du_dp_x = qual * duV_dp + (1_pr - qual) * duL_dp
           dv_dp_x = qual * dvV_dp + (1_pr - qual) * dvL_dp
!
           sound = SQRT((P_tri + ratio)/(du_dp_x - ratio * dv_dp_x)) * vvv_TPM(k,MMM_TPM)
!
           ccc_TPM(k,MMM_TPM) = sound
!        print*,k,'v', vvv_TPM(k,MMM_TPM),'x',qual,'sound',ccc_TPM(k,MMM_TPM),y_mesh_TPM(k) 
        END DO
!stop
!        vvv_TPM(1,MMM_TPM) = v_D
!        xxx_TPM(1,MMM_TPM) = (e_cr - el) / (ev - el)
!
! Left boundary == a part of right boundary of LH
!print*, 'stat left'
        vvv_TPM(1,1) = 1_pr/rho_cr
        ppp_TPM(1,1) = P_cr
        TTT_TPM(1,1) = T_cr
        ccc_TPM(1,1) = c_cr
        xxx_TPM(:,1) = 1_pr
!                 
!        v_six = 2.2530239912767237e-3_pr
!        DO i = 2, 18
!        ! Calculation of the volumes in the linear part for 100*100 mesh
!          vvv_TPM(i,1) = (y_mesh_TPM(i) - e_cr) / ( y_mesh_TPM(19) - e_cr) * &
!     &    (v_six - 1d0 / rho_cr) + 1d0 / rho_cr
!
!          CALL New_Rap1D(1, TTT_TPM(i,1), out_2, res, Niter, &
!     &         flag, y_mesh_TPM(i), T_cr-0.01,&
!     &         vvv_TPM(i,1), out_2)
!
!          CALL pressure(TTT_TPM(i,1),vvv_TPM(i,1),press)
!          ppp_TPM(i,1) = press
!!
!          CALL sound_speed(TTT_TPM(i,1),vvv_TPM(i,1), sound)
!          ccc_TPM(i,1) = sound
!
!          IF (res > 1e-12) THEN
!             print*, "Left in TPM", res, "iter", Niter,i
!             STOP
!          ENDIF
!
!  print*, i, vvv_TPM(i,1), ppp_TPM(i,1),TTT_TPM(i,1),ccc_TPM(i,1)
!        ENDDO
             T_guess = T_cr-0.001_pr
             v_guess = 1_pr / rho_cr + 1e-4_pr
             vliq_guess = 1_pr / rho_cr - 1e-4_pr
!             T_guess = TTT_TPM(18,1)
!             v_guess = vvv_TPM(18,1)
!             vliq_guess = 1_pr / rho_cr - 2e-4_pr
        DO i=2, NNN_TPM
!         
         CALL New_Rap3D(4,vliq,vvv_TPM(i,1),TTT_TPM(i,1),&
     &        res, Niter, flag, y_mesh_TPM(i),in_2,&
     &        vliq_guess,v_guess,T_guess)
!
         v_guess = vvv_TPM(i,1)
         vliq_guess = vliq
         T_guess = TTT_TPM(i,1)
!
         CALL pressure(TTT_TPM(i,1),vvv_TPM(i,1),press)
         ppp_TPM(i,1) = press
!
         CALL sound_speed(TTT_TPM(i,1),vvv_TPM(i,1), sound)
         ccc_TPM(i,1) = sound
!
         IF (res > 1e-7) THEN
            print*, "Left in TPM", res, "iter", Niter
            STOP
         ENDIF
!      print*, i, vvv_TPM(i,1), ppp_TPM(i,1),TTT_TPM(i,1),ccc_TPM(i,1)
        END DO
!STOP
!
!              
! Middle points
!
!
!print*, 'stat middle'
!         DO j=2,MMM_TPM-1
          DO i=1, NNN_TPM
!             pguess = ppp_TPM(i,1)
!          DO i=1, NNN_TPM
             DO j=2,MMM_TPM-1
!
!           vvv_TPM(i,j) = ((x_mesh_TPM(j) - x_mesh_min) /&
!     &    (x_mesh_max - x_mesh_min)) *(vvv_TPM(i,MMM_TPM)-vvv_TPM(i,1))&
!     &     +vvv_TPM(i,1)
           vvv_TPM_log = ((x_mesh_TPM(j) - x_mesh_min) /&
     &     (x_mesh_max - x_mesh_min)) *(log10(vvv_TPM(i,MMM_TPM))-log10(vvv_TPM(i,1)))+&
     &     log10(vvv_TPM(i,1))
           vvv_TPM(i,j) = 10_pr**vvv_TPM_log
!
!      print*, vvv_TPL(i,MMM_TPL), vvv_TPL(i,1)
!      print*, y_mesh_TPL(i),vvv_TPL(i,j),pguess            
!           CALL New_Rap1D(2, press, qual, res1, Niter, flag,&
!     &           y_mesh_TPM(i), pguess, vvv_TPM(i,j), temp)
!
           lower = P_tri
           upper = P_cr - 0.01e6_pr
!           upper = pguess+2.0e6_pr
!print*, i,j,'lower',lower,'upper',upper
           CALL BrentRoots(1, press, temp, res, Niter, y_mesh_TPM(i), lower, upper, vvv_TPM(i,j), dummy)           
!
!print*, 'v',vvv_TPM(i,j),'e',y_mesh_TPM(i), 'p,T',press,temp
!print*, ' '
           ppp_TPM(i,j) = press
           TTT_TPM(i,j) = temp
!           pguess       = press !- 1e4_pr 
!
           IF (res  > 1e-10_pr) THEN
             print*, "middle points TPM", res, "iter", Niter,i,j
             STOP
           ENDIF
!          
! Speed of sound  
!         
           CALL satprop(3, press, dummy, vv, vl, ev, el)
           CALL satderiv(3, press, duL_dp, duV_dp, dvL_dp, dvV_dp)
!
           qual = (vvv_TPM(i,j)-vl) / (vv-vl)
!
           ratio = ((ev - el)/(vv - vl))   ! (J/kg)/(m3/kg)
           du_dp_x = qual * duV_dp + (1_pr - qual) * duL_dp
           dv_dp_x = qual * dvV_dp + (1_pr - qual) * dvL_dp
!
           sound = SQRT((press + ratio)/(du_dp_x - ratio * dv_dp_x)) * vvv_TPM(i,j)
!
           ccc_TPM(i,j) = sound
           xxx_TPM(i,j) = qual
!
!      print*,j,i, ppp_TPM(i,j), TTT_TPM(i,j), xxx_TPM(i,j)
          END DO
         END DO
!STOP
!  print*, 'TPM','1,1', y_mesh_TPM(1),  vvv_TPM(1,1),&
!&                          TTT_TPM   (1,1),ppp_TPM(1,1),ccc_TPM(1,1)
!  print*, 'TPM','2,1   ', y_mesh_TPM(2),  vvv_TPM(2,1),&
!&                          TTT_TPM   (2,1),ppp_TPM(2,1),ccc_TPM(2,1)
!  print*, 'TPM','1,2   ', y_mesh_TPM(1),  vvv_TPM(1,2),&
!&                          TTT_TPM   (1,2),ppp_TPM(1,2),ccc_TPM(1,2)
!
! Spline line coefficient right boundary
!
!print*, 'stat spline right'
        ev = e_tri_R
        el = e_tri_L
        vv = v_tri_R
        vl = v_tri_L 
       DO i=1,NNN_sat_TPM
           qual = (y_mesh_sat_TPM(i) - el) / (ev - el)
           v_right_TPM(i) = qual*vv + (1_pr-qual)*vl
       ENDDO
       DO i = 1, NNN_TPM-1
          j = ord_spline * (i-1) + 1
          CALL polyfit(spline_right_TPM(:,i), &
     &               y_mesh_sat_TPM(j:j+ord_spline),&
     &               v_right_TPM(j:j+ord_spline), ord_spline)
       ENDDO
! 
! Spline line coefficient left boundary
!
!print*, 'stat spline left'
        v_left_TPM(1) = 1_pr/rho_cr

!        T_guess = T_cr - 0.0001_pr
!        v_guess =  1_pr/rho_cr + 1e-5_pr
!        vliq_guess = 1_pr/rho_cr - 1e-5_pr
!
! USE test_grid to compute T_guess, v_guess, vliq_guess
        T_guess = 304.1281 
        v_guess =  2.1443612972628459E-003 
        vliq_guess = 2.1322415401290890E-003 

        DO i=2, NNN_sat_TPM
!              T_guess_sat = T_sat_sat(i-1)
!              v_guess_sat = v_Lsat_LH(i-1)
!              vliq_guess_sat = vliq_sat_LH(i-1)
!           IF( i == 2) THEN
!print*,i, 'in', y_mesh_sat_TPM(i),in_2,vliq_guess,v_guess, T_guess

          CALL New_Rap3D(4,vliq,v_left_TPM(i),Tsat,&
      &        res, Niter, flag, y_mesh_sat_TPM(i),in_2,&
      &        vliq_guess,v_guess, T_guess)
!
!print*, 'out', vliq,v_left_TPM(i),Tsat,res,Niter,flag              

          CALL New_Rap3D(4,vliq,v_left_TPM(i),Tsat,&
      &        res,Niter,flag, y_mesh_sat_TPM(i),in_2,&
      &        vliq,v_left_TPM(i),Tsat)
!
!           ELSE
! 
!               CALL New_Rap3D(4,vliq_sat_LH(i),v_Lsat_LH(i),T_sat_sat(i),&
!              & res_sat_LH(i), Niter, exitflag, y_mesh_sat_LH(i),in_2,&
!              & vliq_guess_sat,v_guess_sat, T_guess_sat)
!           ENDIF

              T_guess    = Tsat
              v_guess    = v_left_TPM(i)
              vliq_guess = vliq
!
              IF (res > 1e-8) THEN
               print*, "Left spline TPM", res, "i", i,&
               &        'T=',Tsat
               STOP
              ENDIF
!      print*, i, v_left_TPM(i)
        END DO
!
!
      DO i = 1, NNN_TPM-1
        j = 3 * (i-1) + 1
        CALL polyfit(spline_right_TPM(:,i), y_mesh_sat_TPM(j:j+3),v_right_TPM(j:j+3),ord_spline)
       CALL polyfit(spline_left_TPM(:,i), y_mesh_sat_TPM(j:j+3),v_left_TPM(j:j+3),ord_spline)
!        CALL polyfit(spline_left_TPM(:,i), y_mesh_sat_TPM(j:j+3),LOG10(v_left_TPM(j:j+3)),ord_spline)
      ENDDO 
!
!     print*,'------------------------------------------------------'
     print*,           'TPM GRID FINISH'
!     print*,'------------------------------------------------------'
!
!
      END SUBROUTINE grid_construction_TPM
