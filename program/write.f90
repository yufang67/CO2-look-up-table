PROGRAM write_solut
!      
      USE def_constants
      USE def_variables
      USE non_linear_solvers
      USE properties
      USE grid_functions
      USE Grid
!      USE Interp_table
      IMPLICIT NONE


      INTEGER :: i,j
      REAL(pr) :: cp


        CALL MAKE_GRID()
!---------------------------------------------------------------------------
!
!                     domain
!
!---------------------------------------------------------------------------
! Save two boundarys
!OPEN (UNIT = 22,             FILE = 'boundary_LH.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!
!DO i = 1, NNN_LH
!     WRITE(22,*), vvv_LH(i,1),vvv_LH(i,MMM_LH), y_mesh_LH(i)
!ENDDO

! CLOSE(UNIT = 22)

!
!
!
! Save spline 
!OPEN (UNIT = 22,             FILE = 'spline_LH.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!
!DO i = 1, NNN_sat_LH
!      WRITE(22,*) v_Lsat_LH(i), v_Lpmax_LH(i),  y_mesh_sat_LH(i)
!ENDDO

! CLOSE(UNIT = 22)
OPEN (UNIT = 22,             FILE = 'spline_LL.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')
!
DO i = 1, NNN_sat_LL
      WRITE(22,*) v_Lsat_LL(i), v_Lpmax_LL(i),  y_mesh_sat_LL(i), v_liq_meta(i)
ENDDO
!
 CLOSE(UNIT = 22)
!OPEN (UNIT = 22,             FILE = 'spline_HT.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!!
!DO i = 1, NNN_sat_HT
!      WRITE(22,*) v_right_HT(i), v_left_HT(i),  y_mesh_sat_HT(i)
!ENDDO
!
! CLOSE(UNIT = 22)
!OPEN (UNIT = 22,             FILE = 'spline_R.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!!
!DO i = 1, NNN_sat_R
!     WRITE(22,*) v_Vsat(i), v_Vpmin(i),  y_mesh_sat_R(i)
!ENDDO
!
! CLOSE(UNIT = 22)
STOP
!
!
! Save e, v, p, t, c
    OPEN (UNIT = 42,             FILE   = 'LH_p_T_c', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_LH
         DO j = 1, MMM_LH
            CALL heat_cap_p(TTT_LH(i,j),vvv_LH(i,j),cp)
            WRITE(42,*) y_mesh_LH(i), vvv_LH(i,j), ppp_LH(i,j), &
     &                                          TTT_LH(i,j), ccc_LH(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)

    OPEN (UNIT = 42,             FILE   = 'LL_p_T_c', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_LL
         DO j = 1, MMM_LL
            CALL heat_cap_p(TTT_LL(i,j),vvv_LL(i,j),cp)
            WRITE(42,*) y_mesh_LL(i), vvv_LL(i,j), ppp_LL(i,j), &
     &                                          TTT_LL(i,j), ccc_LL(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)
!
   OPEN (UNIT = 42,             FILE   = 'HT_p_T_c', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_HT
         DO j = 1, MMM_HT
            CALL heat_cap_p(TTT_HT(i,j),vvv_HT(i,j),cp)
            WRITE(42,*) y_mesh_HT(i), vvv_HT(i,j), ppp_HT(i,j), &
     &                                          TTT_HT(i,j), ccc_HT(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)
!
   OPEN (UNIT = 42,             FILE   = 'R_p_T_c', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_R
         DO j = 1, MMM_R
            CALL heat_cap_p(TTT_R(i,j),vvv_R(i,j),cp)
            WRITE(42,*) y_mesh_R(i), vvv_R(i,j), ppp_R(i,j), &
     &                                          TTT_R(i,j), ccc_R(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)





END PROGRAM write_solut
