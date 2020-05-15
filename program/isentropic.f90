PROGRAM isen
!      
!      USE def_constants
!      USE def_variables
      USE axl_solvers
      USE properties
!      USE grid_functions
      USE Grid
      USE Interp_table
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
      REAL(8) :: p_in, T_in, s_in, h_in, u_in, flux_in, rho_in, v_in,Ain,e_in
      REAL(8) :: eint, volum,gus1, gus2,h,s,psat,qual 
      INTEGER :: Niter,exitflag,i,Num,iocheck,flag
      REAL(8), DIMENSION(1000) :: Adiv,umoin
      REAL(8) :: out_2,resnorm,v_guess, out3, htot, theta, divL,deltaL
      REAL(8) :: u,rho,dummy,lower, upper, rhoinv
      REAL(8) :: tt, vvt, vlt, evt, elt,st,ms,clt,ut,p1,p2,p3,dA,amoin
      REAL(8), DIMENSION(104) :: x, pmid, rmid, emid, Tmid, cmid, xmid, amid, smid
      
      CALL MAKE_GRID()
!-------test-----
!      psat = 6.235379d6
!      psat = 6.6177756d6
!      CALL satprop(3, psat, tt, vvt, vlt, evt, elt)
!      print*, 'temp', tt

!  open (10, file = 'profile_p.txt', form = 'formatted', status = 'old', iostat = iocheck)
  
 ! i = 1
!  do while (.true.)
!!     if (iocheck .ne. 0) exit
!  DO i = 1,104
!     read (10, *, iostat = iocheck) x(i), pmid(i), rmid(i), emid(i) 
!     print*, x(i), pmid(i), rmid(i), emid(i), i 
!!     if (iocheck .ne. 0) exit
!     i = i + 1
!  end do
!  
!  close(10)
!  
!  DO i = 1,104
!     
!     CALL CO2BLLT_EQUI(pmid(i), Tmid(i), cmid(i), xmid(i), amid(i), smid(i), emid(i),1.0/rmid(i),flag)
!    
!     print*, 'i',i,smid(i), flag
!  ENDDO   
!
!STOP
!--------------------------------
! INLET NOZZLE
      p_in    = 9.1d6
!      T_in    = 310.15   !37
      T_in    = 309.65    !36.5
      u_in    = 0.5
      v_guess = 1.0 / 600.0
      Num     = 1000
! GEO
      Ain     = 10d-3 
      divL    = 56.15d-3
!      theta   = 0.076 / 180.0 * 3.1415926    ! RAD
      theta   = 0.612 / 180.0 * 3.1415926    ! RAD
      deltaL   = divL / int(Num)
      DO i = 1,Num
         Adiv(i) = (0.12d-3 + deltaL*(i-1)*tan(theta) ) * 2.0
!         print*, 'section', Adiv(i)
      ENDDO
!STOP
!        
!     CALL New_Rap1Daxl(3, v_in, out_2, resnorm, Niter,&
!     &              exitflag, p_in,v_guess , T_in, out3)
!
!     print*, "density inlet", 1.0/v_in, resnorm
!stop
!---------------------------- INITIAL ----------------------------------
!      rho_in  = 621.475   !nozz1 
      rho_in  = 634.82     !nozz4
      v_in    = 1.0 / rho_in     
!      flux_in = rho_in*u_in*Ain 
!
      CALL entropy(T_in, v_in, s_in)
      CALL inter_energy(T_in, v_in, e_in)
      h_in = e_in + p_in*v_in
!      print*, 'entropy_in', s_in, 'enthalpy_in',h_in
      htot = h_in !+ 0.5*u_in*u_in
!
!---------------------------- THROAT ----------------------------------
!
       psat = 7.0d6
       CALL satprop(3, psat, tt, vvt, vlt, evt, elt)
       CALL axlpress(tt,vvt,vlt,p1,p2,p3)
!print*, 'psat', psat,p1,p2,p3, (psat+p1+p2+p3)/4.0

       ut   = sqrt(htot - elt - psat*vlt) * 2.0 
       CALL sound_speed(tt, vlt,clt)
       CALL entropy(tt,vlt,s_in)
!       print*, 'c',clt,'velo',ut,'s',s_in
       flux_in = ut*Adiv(1) / vlt
!! ut is supersonic
!STOP
!--------------------------------------------------------------------
!       lower = 6.5d6
!       upper = 7.5d6
!!
!       CALL BrentRootsaxl(2, psat, dummy, resnorm, Niter, &
!&           htot, lower, upper, flux_in/Adiv(1), dummy)
!
!       print*, 'throat psat',psat, 'res',resnorm, 'Niter',Niter
!
!       CALL satprop(3, psat, tt, vvt, vlt, evt, elt)
!       
!       CALL entropy(tt,vlt,st)
!
!       print*, 'throat entropy',st
!       CALL sound_speed(tt, vlt,clt)
!       print*, 'c',clt,'velo',flux_in/Adiv(1)*vlt
!  STOP
!--------------------------- In divergence --------------------------
!!! 
      gus1    = psat - 0.1d6
      gus2    = 1.623d-3
      umoin(1)= ut

      DO i = 2,Num

!        ms = flux_in/Adiv(i)         
!         u = flux_in / (rho*Adiv(i))
!         h = htot - 0.5*u*u
!         s = s_in
  
!         CALL New_Rap2Daxl(4, psat, rhoinv, &
!&             resnorm, Niter, exitflag, s_in, htot,ms,dummy,gus1, gus2) 
!            gus1    = psat 
!            gus2    = rhoinv 
         Amoin = (Adiv(i)+Adiv(i-1)) / 2.0
         dA   = ( Adiv(i) - Adiv(i-1) ) / amoin
         CALL New_Rap2Daxl(4, psat, rhoinv, &
&             resnorm, Niter, exitflag, s_in, htot,dA,umoin(i-1),gus1, gus2) 
         
         umoin(i) = flux_in *rhoinv / Adiv(i)
         gus1 = psat  !- 0.1d6
         gus2 = rhoinv



         print*, 'i',i,'res', resnorm,Niter
         print*, 'psat', psat, 'density', 1.0/rhoinv
!stop
      ENDDO










END PROGRAM isen
