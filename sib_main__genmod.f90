        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:31 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SIB_MAIN__genmod
          INTERFACE 
            SUBROUTINE SIB_MAIN(GREF,LONSIB,LATSIB,PREF,PNUM,PHYSCONT,  &
     &GPROGT,GDIAGT,SIBGL)
              USE MODULE_SIB, ONLY :                                    &
     &          GPROG_TYPE,                                             &
     &          GDIAG_TYPE,                                             &
     &          LU_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHYS_PARAM
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: LONSIB
              REAL(KIND=4), INTENT(IN) :: LATSIB
              INTEGER(KIND=4), INTENT(IN) :: PREF
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
              TYPE (GPROG_TYPE), INTENT(INOUT) :: GPROGT
              TYPE (GDIAG_TYPE), INTENT(INOUT) :: GDIAGT
              TYPE (LU_TYPE), INTENT(INOUT) :: SIBGL
            END SUBROUTINE SIB_MAIN
          END INTERFACE 
        END MODULE SIB_MAIN__genmod
