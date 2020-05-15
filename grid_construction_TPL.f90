!============================================================================
      SUBROUTINE grid_construction_TPL()
!============================================================================
! Region two-phase (TP)
! 
! The specific internal energy along from Triple point to Critical
! point. The region is bounded by the saturation line and the triple
! pressure line, the minimum internal energy and the critical internal
! energy. 
! NNN_LL = NNN_TPL
! 
! The bilinear mapping is used for the interpolation from physical (v,e)
! doamin to the transformed domaine (X,Y).(Bicubic can be
! considered later by using also the derivatives)
!========================================================================
!
     USE def_constants
     USE def_variables, ONLY: x_mesh_TPL,y_mesh_TPL,y_mesh_Ptri_TPL, v_right_TPL,&
&        spline_right_TPL,ppp_TPL, vvv_TPL, xxx_TPL, ccc_TPL, TTT_TPL,&
&        vvv_LL, ppp_LL, ccc_LL, TTT_LL
!     &   vL_psat_spline, vV_psat_spline,uL_psat_spline, uV_psat_spline,&
!     &   Tsat_psat_spline, SaturP
     USE non_linear_solvers
     USE properties
     USE grid_functions
!        
      IMPLICIT NONE
!     
      REAL(pr) :: u_max, u_min , delta, deltap,vvv_TPL_log
      REAL(pr) :: ev,el,vv,vl,qual,sound,pp,ratio,press,temp
      INTEGER  :: i,k,j_tpl, Niter, flag,j
      REAL(pr) :: duL_dp, duV_dp, dvL_dp, dvV_dp,du_dp_x,dv_dp_x
      REAL(pr) :: pguess,res,vvguess,vlguess,Tguess,dummy
      REAL(pr) :: T_guess,v_guess,vap_guess,in_2,lower,upper
      REAL(pr), DIMENSION(NNN_TPL) :: vap_TPL



!--------------------------------------------------------------------------
!
! Construct vector y_mesh_TPL(physical domain), x_mesh_TPL(tansformed
! domain)
!
!-------------------------------------------------------------------------
        u_max = e_cr
        u_min = e_tri_L
!
        delta = (u_max - u_min) / (NNN_TPL - 1)
        y_mesh_TPL = u_min + (/(i*delta,i=0,NNN_TPL-1)/)
! 
        delta = (x_mesh_max - x_mesh_min) / (MMM_TPL -1)
        x_mesh_TPL = x_mesh_min + (/(i*delta, i=0,MMM_TPL-1)/)
!
        delta             = (u_max - u_min)/(NNN_isop_TPL-1)
        y_mesh_Ptri_TPL   = u_min + (/(i*delta, i=0,NNN_isop_TPL-1)/)
!
!
! Bottom point
         ccc_TPL(1,:) = c_tri_L
         ppp_TPL(1,:) = p_tri
         vvv_TPL(1,:) = v_tri_L
         TTT_TPL(1,:) = T_tri
         xxx_TPL(1,:) = 0_pr
!
! Right boundary = iso triple pressure line
        ev = e_tri_R
        el = e_tri_L
        vv = v_tri_R
        vl = v_tri_L
        ppp_TPL(:,MMM_TPL) = P_tri
        TTT_TPL(:,MMM_TPL) = T_tri
        DO k=2,NNN_TPL
!        
           xxx_TPL(k,MMM_TPL) = (y_mesh_TPL(k) - el) / (ev - el)
           qual = (y_mesh_TPL(k) - el) / (ev - el)
           vvv_TPL(k,MMM_TPL) = qual*vv + (1_pr-qual)*vl 
! Speed of sound  
!
           CALL satderiv(3, ppp_TPL(k,MMM_TPL), duL_dp, duV_dp, dvL_dp, dvV_dp)
!
           ratio = ((ev - el)/(vv - vl))   ! (J/kg)/(m3/kg)
           du_dp_x = qual * duV_dp + (1_pr - qual) * duL_dp
           dv_dp_x = qual * dvV_dp + (1_pr - qual) * dvL_dp
!
           sound = SQRT((P_tri + ratio)/(du_dp_x - ratio * dv_dp_x)) * vvv_TPL(k,MMM_TPL)
!
           ccc_TPL(k,MMM_TPL) = sound
!      print*,k,'v', vvv_TPL(k,MMM_TPL),'x',qual,'sound',ccc_TPL(k,MMM_TPL),y_mesh_TPL(k)
        END DO
!      STOP
!
!
! Left boundary
!
        vvv_TPL(NNN_TPL,1) = 1_pr/rho_cr
        TTT_TPL(NNN_TPL,1) = T_cr
        ppp_TPL(NNN_TPL,1) = P_cr
        ccc_TPL(NNN_TPL,1) = c_cr
        xxx_TPL(:,1)       = 0_pr
!
         T_guess   = T_tri
         v_guess   = 1_pr / rho_tri_L
         vap_guess = 1_pr / rho_tri_R
          DO i=2, NNN_TPL-1
!         
              CALL New_Rap3D(1,vvv_TPL(i,1),vap_TPL(i),TTT_TPL(i,1),&
              & res, Niter, flag, y_mesh_TPL(i),in_2,&
              & v_guess, vap_guess,T_guess)
!
              T_guess   = TTT_TPL(i,1)
              v_guess   = vvv_TPL(i,1)
              vap_guess = vap_TPL(i)

              CALL pressure(TTT_TPL(i,1),vvv_TPL(i,1),press)
              ppp_TPL(i,1) = press
!
              CALL sound_speed(TTT_TPL(i,1),vvv_TPL(i,1), sound)
              ccc_TPL(i,1) = sound
!
!        print*, "vap_v", vap_LL(i)
              IF (res > res_ref) THEN
                print*, "res in TPL left", res, "iter", Niter
                print*, 'i',i,'v',vvv_TPL(i,1),'vap',vap_TPL(i),'e',y_mesh_TPL(i)
                STOP
              ENDIF
!
!        print*, "left boundary", vvv_TPL(i,1), 't_right',TTT_TPL(i,1)
!        print*, 'e', y_mesh_TPL(i),i
        END DO

!STOP
!              
! Middle points

!         DO j=2,MMM_TPL
          DO i=2, NNN_TPL
!              pguess = 0.6e6_pr
              pguess = ppp_TPL(i,1) - 0.0001e6_pr
!              pguess = P_tri + 0.00001e6_pr
!          DO i=2, NNN_TPL
           DO j=2,MMM_TPL-1
!
!           vvv_TPL(i,j) = ((x_mesh_TPL(j) - x_mesh_min) /&
!     &    (x_mesh_max - x_mesh_min)) *(vvv_TPL(i,MMM_TPL)-vvv_TPL(i,1))&
!     &     +vvv_TPL(i,1)
           vvv_TPL_log  = ((x_mesh_TPL(j) - x_mesh_min) /&
     &     (x_mesh_max - x_mesh_min)) *(log10(vvv_TPL(i,MMM_TPL))-log10(vvv_TPL(i,1)))+&
     &      log10(vvv_TPL(i,1))
           vvv_TPL(i,j) = 10_pr**vvv_TPL_log
!
!      print*, vvv_TPL(i,MMM_TPL), vvv_TPL(i,1)
!      print*,i,j, y_mesh_TPL(i),vvv_TPL(i,j),pguess            
!           CALL New_Rap1D(2, press, qual, res, Niter, flag,&
!     &                    y_mesh_TPL(i), pguess, vvv_TPL(i,j), temp)
!
           lower = p_tri 
           upper = p_cr-0.001e6_pr
!           upper = pguess + 1.0e6_pr
!      print*, lower,upper, vvv_TPL(i,j),y_mesh_TPL(i)
           CALL BrentRoots(1, press, temp, res, Niter, y_mesh_TPL(i), lower, upper, vvv_TPL(i,j), dummy)
!
            IF (res  > 1e-10) THEN
             print*, "middle points TPL", res, "iter", Niter,i,j
             print*, 'v',vvv_TPL(i,j),'u',y_mesh_TPL(i), pguess, press
             STOP
           ENDIF
           ppp_TPL(i,j) = press
           TTT_TPL(i,j) = temp
           pguess       = press !- 0.0003d6
!           pguess       = press + 0.00014d6
!          
! Speed of sound  
!       
           CALL satprop(3, press, dummy, vv, vl, ev, el) 
           CALL satderiv(3, press, duL_dp, duV_dp, dvL_dp, dvV_dp)
!
           qual = (vvv_TPL(i,j)-vl) / (vv-vl)
           ratio = ((ev - el)/(vv - vl))   ! (J/kg)/(m3/kg)
           du_dp_x = qual * duV_dp + (1_pr - qual) * duL_dp
           dv_dp_x = qual * dvV_dp + (1_pr - qual) * dvL_dp
!
           sound = SQRT((press + ratio)/(du_dp_x - ratio * dv_dp_x)) * vvv_TPL(i,j)
!
           xxx_TPL(i,j) = qual
           ccc_TPL(i,j) = sound
!           xxx_TPL(i,j) = qual
!
!      print*,j,i, 'p,t,x,c',ppp_TPL(i,j), TTT_TPL(i,j), xxx_TPL(i,j),ccc_TPL(i,j)
!      print*,j,i, press + ratio,du_dp_x - ratio * dv_dp_x
!      print*,j,i, 'ev,el,vv,vl',ev,el, vv, vl
!      STOP 
          END DO
!      STOP
         END DO

!  print*, 'TPL','max-1,1', y_mesh_TPL(NNN_TPL-1),  vvv_TPL(NNN_TPL-1,1),&
!&                          TTT_TPL   (NNN_TPL-1,1),ppp_TPL(NNN_TPL-1,1),ccc_TPL(NNN_TPL-1,1)
!  print*, 'TPL','max-1,2   ', y_mesh_TPL(NNN_TPL-1),  vvv_TPL(NNN_TPL-1,2),&
!&                          TTT_TPL   (NNN_TPL-1,2),ppp_TPL(NNN_TPL-1,2),ccc_TPL(NNN_TPL-1,2)
!  print*, 'TPL','max,2   ', y_mesh_TPL(NNN_TPL),  vvv_TPL(NNN_TPL,2),&
!&                          TTT_TPL   (NNN_TPL,2),ppp_TPL(NNN_TPL,2),ccc_TPL(NNN_TPL,2)


!
! Spline line coefficient right boundary
!
        ev = e_tri_R
        el = e_tri_L
        vv = v_tri_R
        vl = v_tri_L 
       DO i=1,NNN_isop_TPL
           qual = (y_mesh_Ptri_TPL(i) - el) / (ev - el)
           v_right_TPL(i) = qual*vv + (1_pr-qual)*vl
       ENDDO
       DO i = 1, NNN_TPL-1
          j = ord_spline * (i-1) + 1
          CALL polyfit(spline_right_TPL(:,i), &
     &               y_mesh_Ptri_TPL(j:j+ord_spline),&
     &               v_right_TPL(j:j+ord_spline), ord_spline)
       ENDDO
!     
!     print*,'------------------------------------------------------'
     print*,           'TPL GRID FINISH'
!     print*,'------------------------------------------------------'
!
!
      END SUBROUTINE grid_construction_TPL
