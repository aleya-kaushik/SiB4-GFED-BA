        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:02 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADDINC__genmod
          INTERFACE 
            SUBROUTINE ADDINC(GREF,PFTREF,LONSIB,LATSIB,GDIAGT,GPROGT,  &
     &CAST,FLUXT,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          GDIAG_TYPE,                                             &
     &          GPROG_TYPE,                                             &
     &          CAS_TYPE,                                               &
     &          FLUX_TYPE,                                              &
     &          SSCOL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: GREF
              INTEGER(KIND=4), INTENT(IN) :: PFTREF
              REAL(KIND=4), INTENT(IN) :: LONSIB
              REAL(KIND=4), INTENT(IN) :: LATSIB
              TYPE (GDIAG_TYPE), INTENT(INOUT) :: GDIAGT
              TYPE (GPROG_TYPE), INTENT(INOUT) :: GPROGT
              TYPE (CAS_TYPE), INTENT(INOUT) :: CAST
              TYPE (FLUX_TYPE), INTENT(INOUT) :: FLUXT
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE ADDINC
          END INTERFACE 
        END MODULE ADDINC__genmod
