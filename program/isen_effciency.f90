PROGRAM isen
!      
!      USE def_constants
!      USE def_variables
      USE axl_solvers
      USE properties
!      USE grid_functions
      USE Grid
!      USE Interp_table
!      USE derivees
!      USE interp_functions
!      USE saturation
      IMPLICIT NONE
!======================================================================
!
!       isentropic expension
!
!======================================================================
!
      REAL(8) :: p_in, T_in, s_in, h_in, u_in,  rho_in, v_in,Ain,e_in,vv_in1
      REAL(8) :: Aout,pout,qualout,uout,hout,vvout
      REAL(8) :: htot
      REAL(8) :: psat_guess, qual_guess,dummy3,dummy4,resnorm
      INTEGER :: Niter,exitflag
      REAL(8) :: tt, vvt, vlt, evt, elt

      
      REAL(8) :: p_in2, T_in2, s_in2, h_in2, u_in2,  rho_in2, v_in2,e_in2,vv_in2
      REAL(8) :: pout2,qualout2,uout2,hout2,vv_out2
      REAL(8) :: htot2
      REAL(8) :: psat_guess2, qual_guess2


      CALL MAKE_GRID()

      
      
!---------------------------- INITIAL ----------------------------------
      p_in    = 8.04d6
      T_in    = 299.65 
      rho_in  = 759.13
      v_in    = 1.0 / rho_in     
      u_in    = 37.735
      vv_in1  = 0 !0.21
      Ain     = 6d-3 
!
      CALL entropy(T_in, v_in, s_in)
      CALL inter_energy(T_in, v_in, e_in)
!
      h_in = e_in + p_in*v_in
      htot = h_in + 0.5*(u_in*u_in+vv_in1*vv_in1)
      s_in = s_in

print*, 'entropy_in', s_in, 'htot',htot
!
!---------------------------- OUTLET ----------------------------------
      Aout    = 3.11d-3
      uout    = 83.02
      vvout    = 0 !0.56
      hout    = htot - 0.5*(uout*uout+vvout*vvout)
!
!       CALL BrentRootsaxl(2, psat, dummy, resnorm, Niter, &
!&           htot, lower, upper, flux_in/Adiv(1), dummy)
!
!       print*, 'throat psat',psat, 'res',resnorm, 'Niter',Niter
!
!       
!       CALL entropy(tt,vlt,st)
!
!       print*, 'throat entropy',st
!       CALL sound_speed(tt, vlt,clt)
!       print*, 'c',clt,'velo',flux_in/Adiv(1)*vlt
!  STOP
!!! 
      psat_guess  = 5.6805d6
      qual_guess  = 0.0694

print*, 'hout',hout  
!         CALL New_Rap2Daxl(4, psat, rhoinv, &
!&             resnorm, Niter, exitflag, s_in, htot,ms,dummy,gus1, gus2) 
!            gus1    = psat 
!            gus2    = rhoinv 
         CALL New_Rap2Daxl(5, pout, qualout, &
&             resnorm, Niter, exitflag, s_in, hout,dummy3,dummy4,psat_guess, qual_guess) 
    

!         psat_guess  = pout
!         qual_guess  = qualout

!         CALL New_Rap2Daxl(5, pout, qualout, &
!&             resnorm, Niter, exitflag, s_in, hout,dummy3,dummy4,psat_guess, qual_guess) 
         

       CALL satprop(3, pout, tt, vvt, vlt, evt, elt)

         print*,'res', resnorm,Niter
         print*, 'psat1', pout, 'Tsat1',tt ,'qual1',qualout 



!---------------------------- INITIAL ----------------------------------
      p_in2    = 5.5532d6
      T_in2    = 293.073  
      rho_in2  = 177.519
      v_in2    = 1.0 / rho_in2     
      u_in2    = 3.5205
      vv_in2   = 0.725

!
      CALL entropy(T_in2, v_in2, s_in2)
      CALL inter_energy(T_in2, v_in2, e_in2)
!
      h_in2 = e_in2 + p_in2*v_in2
      htot2 = h_in2 + 0.5*(u_in2*u_in2+vv_in2*vv_in2)
      s_in2 = s_in2

print*, 'entropy_in', s_in2, 'htot',htot2
!
!---------------------------- OUTLET ----------------------------------
      uout2    = 5.46
      vv_out2  = 4.379
      hout2    = htot2 - 0.5*(uout2*uout2+vv_out2*vv_out2)

      psat_guess2  = 5.53969d6
      qual_guess2  = 0.721

print*, 'hout',hout2  
!         CALL New_Rap2Daxl(4, psat, rhoinv, &
!&             resnorm, Niter, exitflag, s_in, htot,ms,dummy,gus1, gus2) 
!            gus1    = psat 
!            gus2    = rhoinv 
         CALL New_Rap2Daxl(5, pout2, qualout2, &
&             resnorm, Niter, exitflag, s_in2, hout2,dummy3,dummy4,psat_guess2, qual_guess2) 
    
       CALL satprop(3, pout2, tt, vvt, vlt, evt, elt)

         print*,'res', resnorm,Niter
         print*, 'psat2', pout2, 'Tsat2',tt ,'qual2',qualout2 

END PROGRAM isen
