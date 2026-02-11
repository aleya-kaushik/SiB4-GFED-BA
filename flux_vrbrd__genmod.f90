        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:21 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FLUX_VRBRD__genmod
          INTERFACE 
            SUBROUTINE FLUX_VRBRD(PNUM,Z2,NSL,RST,TD,CAST,VEGT,GPROGT,  &
     &GDIAGT,FLUXT,HYDROST)
              USE MODULE_SIB, ONLY :                                    &
     &          GPROG_TYPE,                                             &
     &          GDIAG_TYPE,                                             &
     &          CAS_TYPE,                                               &
     &          VEG_TYPE,                                               &
     &          FLUX_TYPE,                                              &
     &          HYDROS_TYPE
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              REAL(KIND=4), INTENT(IN) :: Z2
              INTEGER(KIND=1), INTENT(IN) :: NSL
              REAL(KIND=8), INTENT(IN) :: RST
              REAL(KIND=8), INTENT(IN) :: TD
              TYPE (CAS_TYPE), INTENT(IN) :: CAST
              TYPE (VEG_TYPE), INTENT(IN) :: VEGT
              TYPE (GPROG_TYPE), INTENT(IN) :: GPROGT
              TYPE (GDIAG_TYPE), INTENT(IN) :: GDIAGT
              TYPE (FLUX_TYPE), INTENT(INOUT) :: FLUXT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
            END SUBROUTINE FLUX_VRBRD
          END INTERFACE 
        END MODULE FLUX_VRBRD__genmod
