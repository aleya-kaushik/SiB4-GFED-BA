        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:50 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TM5_INTERP__genmod
          INTERFACE 
            SUBROUTINE TM5_INTERP(INDX,LON,LAT,GDIAGT,GPROGT)
              USE MODULE_SIB, ONLY :                                    &
     &          GDIAG_TYPE,                                             &
     &          GPROG_TYPE
              INTEGER(KIND=4), INTENT(IN) :: INDX
              REAL(KIND=4), INTENT(IN) :: LON
              REAL(KIND=4), INTENT(IN) :: LAT
              TYPE (GDIAG_TYPE) :: GDIAGT
              TYPE (GPROG_TYPE) :: GPROGT
            END SUBROUTINE TM5_INTERP
          END INTERFACE 
        END MODULE TM5_INTERP__genmod
