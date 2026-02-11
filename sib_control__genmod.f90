        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SIB_CONTROL__genmod
          INTERFACE 
            SUBROUTINE SIB_CONTROL(GNUM,GREF,LONSIB,LATSIB,DAYNEW,      &
     &DAYLMAX,DOY,SIBG)
              USE MODULE_SIB, ONLY :                                    &
     &          GRIDCELL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: GNUM
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: LONSIB
              REAL(KIND=4), INTENT(IN) :: LATSIB
              LOGICAL(KIND=4), INTENT(IN) :: DAYNEW
              REAL(KIND=4), INTENT(IN) :: DAYLMAX
              INTEGER(KIND=4), INTENT(IN) :: DOY
              TYPE (GRIDCELL_TYPE), INTENT(INOUT) :: SIBG
            END SUBROUTINE SIB_CONTROL
          END INTERFACE 
        END MODULE SIB_CONTROL__genmod
