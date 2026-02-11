        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:23 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AEROINTERPOLATE__genmod
          INTERFACE 
            SUBROUTINE AEROINTERPOLATE(GREF,GLON,GLAT,PNUM,PREF,LAI,    &
     &FVCOVER,Z0,ZP_DISP,RBC,RDC)
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: GLON
              REAL(KIND=4), INTENT(IN) :: GLAT
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: FVCOVER
              REAL(KIND=8), INTENT(INOUT) :: Z0
              REAL(KIND=8), INTENT(INOUT) :: ZP_DISP
              REAL(KIND=8), INTENT(INOUT) :: RBC
              REAL(KIND=8), INTENT(INOUT) :: RDC
            END SUBROUTINE AEROINTERPOLATE
          END INTERFACE 
        END MODULE AEROINTERPOLATE__genmod
