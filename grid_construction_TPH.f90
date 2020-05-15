!============================================================================
      SUBROUTINE grid_construction_TPH()
!============================================================================
! Region two-phase (TP)
! 
! The specific internal energy along from triple to umax. 
! The region is bounded by the saturation line LH and the left of R
! the triple point internal, point E energy. 
! 
!NNN_TPH = NNN_R 
! The bilinear mapping is used for the interpolation from physical (v,e)
! doamin to the transformed domaine (X,Y).(Bicubic can be
! considered later by using also the derivatives)
!========================================================================
!
     USE def_constants
     USE def_variables, ONLY: x_mesh_TPH,y_mesh_TPH,y_mesh_sat_TPH,v_left_TPH,&
&        spline_left_TPH,saturP,ppp_TPH,vvv_TPH,xxx_TPH,ccc_TPH,TTT_TPH,&
&        ppp_R,  TTT_R,  vvv_R,  ccc_R
!&        vL_psat_spline, vV_psat_spline,uL_psat_spline, uV_psat_spline, Tsat_psat_spline
     USE non_linear_solvers
     USE properties
     USE grid_functions
!        
      IMPLICIT NONE
!     
      REAL(pr) :: u_max, u_min , delta, deltap, vvv_TPH_log
      REAL(pr) :: ev,el,vv,vl,qual,sound,pp,ratio,press,temp,Tsat
      INTEGER  :: i,k,j_tpl, Niter, flag,j
      REAL(pr) :: duL_dp, duV_dp, dvL_dp, dvV_dp,du_dp_x,dv_dp_x
      REAL(pr) :: pguess,res, T_guess,v_six,in_2,out_2,v_guess,vliq,vliq_guess,pguess1i,dummy,pguess1




!--------------------------------------------------------------------------
!
! Construct vector y_mesh_TPL(physical domain), x_mesh_TPL(tansformed
! domain)
!
!-------------------------------------------------------------------------
        u_max = e_umax
        u_min = e_tri_R
!
        delta = (u_max - u_min) / (NNN_TPH - 1)
        y_mesh_TPH = u_min + (/(i*delta,i=0,NNN_TPH-1)/)
! 
        delta = (x_mesh_max - x_mesh_min) / (MMM_TPH -1)
        x_mesh_TPH = x_mesh_min + (/(i*delta, i=0,MMM_TPH-1)/)
!
        delta             = (u_max - u_min)/(NNN_sat_TPH-1)
        y_mesh_sat_TPH   = u_min + (/(i*delta, i=0,NNN_sat_TPH-1)/)
!
       
! Right boundary = saturation line of R
!
        vvv_TPH(:,MMM_TPH) = vvv_R(:,1)
        ppp_TPH(:,MMM_TPH) = ppp_R(:,1)
        TTT_TPH(:,MMM_TPH) = TTT_R(:,1)
        ccc_TPH(:,MMM_TPH) = ccc_R(:,1)
        xxx_TPH(:,MMM_TPH) = 1_pr
!DO i=1,NNN_TPH
!print*, i, 'eint', y_mesh_TPH(i)
!print*, vvv_TPH(i,MMM_TPH), TTT_TPH(i,MMM_TPH)
!ENDDO
!
! Left boundary == a part of right boundary of LH
        vvv_TPH(NNN_TPH,1) = v_umax
        ppp_TPH(NNN_TPH,1) = P_umax
        TTT_TPH(NNN_TPH,1) = T_umax
        ccc_TPH(NNN_TPH,1) = c_umax
!
        T_guess    = T_E 
        v_guess    = v_E
        vliq_guess = 1.117e-3_pr
        DO i=1, NNN_TPH-1
!         
         CALL New_Rap3D(4,vliq,vvv_TPH(i,1),TTT_TPH(i,1),&
     &        res, Niter, flag, y_mesh_TPH(i),in_2,&
     &        vliq_guess,v_guess,T_guess)
!
         v_guess = vvv_TPH(i,1)
         vliq_guess = vliq
         T_guess = TTT_TPH(i,1)
!
         CALL pressure(TTT_TPH(i,1),vvv_TPH(i,1),press)
         ppp_TPH(i,1) = press
!
         CALL sound_speed(TTT_TPH(i,1),vvv_TPH(i,1), sound)
         ccc_TPH(i,1) = sound
!
         IF (res > 1e-7) THEN
            print*, "Left in TPH", res, "iter", Niter
            STOP
         ENDIF
!       print*, 'e_mesh', y_mesh_TPH(i)
!      print*, i, vvv_TPH(i,1), ppp_TPH(i,1),TTT_TPH(i,1),ccc_TPH(i,1)
        END DO
        xxx_TPH(:,1) = 1_pr
!      STOP
!
! Middle points
         pguess1 = 3.8e6_pr
         DO j=2,MMM_TPH-1
          pguess = pguess1
          DO i=1, NNN_TPH-1
!
!           vvv_TPH(i,j) = ((x_mesh_TPH(j) - x_mesh_min) /&
!     &    (x_mesh_max - x_mesh_min)) *(vvv_TPH(i,MMM_TPH)-vvv_TPH(i,1))&
!     &     +vvv_TPH(i,1)
           vvv_TPH_log = ((x_mesh_TPH(j) - x_mesh_min) /&
                & (x_mesh_max - x_mesh_min)) *(log10(vvv_TPH(i,MMM_TPH))-log10(vvv_TPH(i,1)))+&
                & log10(vvv_TPH(i,1))
           vvv_TPH(i,j) = 10_pr**vvv_TPH_log
!
!      print*, vvv_TPL(i,MMM_TPL), vvv_TPL(i,1)
!      print*, y_mesh_TPL(i),vvv_TPL(i,j),pguess            
           CALL New_Rap1D(2, press, qual, res, Niter, flag,&
     &           y_mesh_TPH(i), pguess, vvv_TPH(i,j), temp)
!
           ppp_TPH(i,j) = press
           TTT_TPH(i,j) = temp
           xxx_TPH(i,j) = qual
           pguess       = press
           IF (i==2) THEN
            pguess1 = press
           ENDIF 
!
           IF (res  > 1e-7_pr) THEN
             print*, "middle points TPH", res, "iter", Niter,i,j
             STOP
           ENDIF
!          
! Speed of sound  
!         
           CALL satprop(3, press, dummy, vv, vl, ev, el)
           CALL satderiv(3, press, duL_dp, duV_dp, dvL_dp, dvV_dp)
         
!
!           qual = (vvv_TPH(i,j)-vl) / (vv-vl)
!
           ratio = ((ev - el)/(vv - vl))   ! (J/kg)/(m3/kg)
           du_dp_x = qual * duV_dp + (1_pr - qual) * duL_dp
           dv_dp_x = qual * dvV_dp + (1_pr - qual) * dvL_dp
!
           sound = SQRT((press + ratio)/(du_dp_x - ratio * dv_dp_x)) * vvv_TPH(i,j)
!
           ccc_TPH(i,j) = sound
           xxx_TPH(i,j) = qual

!
!      print*,j,i, ppp_TPH(i,j), TTT_TPH(i,j), xxx_TPH(i,j)
         END DO
        END DO
!STOP
!
! Spline line coefficient left boundary
!
        v_left_TPH(1) = v_E
!
        T_guess = T_E
        v_guess = v_E
        vliq_guess = 1.117e-3_pr
!
        DO i=2, NNN_sat_TPH
          CALL New_Rap3D(4,vliq,v_left_TPH(i),Tsat,&
      &        res, Niter, flag, y_mesh_sat_TPH(i),in_2,&
      &        vliq_guess,v_guess, T_guess)
!              
          CALL New_Rap3D(4,vliq,v_left_TPH(i),Tsat,&
      &        res,Niter,flag, y_mesh_sat_TPH(i),in_2,&
      &        vliq,v_left_TPH(i),Tsat)
!
              T_guess    = Tsat
              v_guess    = v_left_TPH(i)
              vliq_guess = vliq
!!
              IF (res > 1e-7) THEN
               print*, "Left spline TPH", res, "i", i,&
               &        'T=',Tsat
               STOP
              ENDIF
!      print*, i, v_left_TPH(i)
        END DO
!
!
      DO i = 1, NNN_TPH-1
        j = 3 * (i-1) + 1
       CALL polyfit(spline_left_TPH(:,i), y_mesh_sat_TPH(j:j+3),v_left_TPH(j:j+3),ord_spline)
!        CALL polyfit(spline_left_TPH(:,i), y_mesh_sat_TPH(j:j+3),LOG10(v_left_TPH(j:j+3)),ord_spline)
      ENDDO 
!
!     print*,'------------------------------------------------------'
     print*,           'TPH GRID FINISH'
!     print*,'------------------------------------------------------'
!
!
      END SUBROUTINE grid_construction_TPH
