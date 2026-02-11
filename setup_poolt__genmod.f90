        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_POOLT__genmod
          INTERFACE 
            SUBROUTINE SETUP_POOLT(POOLCONT,ROOTF,EQUIBDT,POOLDT,EQUIBLT&
     &,POOLLT,FRACT)
              USE MODULE_SIB, ONLY :                                    &
     &          EQUIBD_TYPE,                                            &
     &          EQUIBL_TYPE,                                            &
     &          POOLD_TYPE,                                             &
     &          POOLL_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: ROOTF(10)
              TYPE (EQUIBD_TYPE), INTENT(INOUT) :: EQUIBDT
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (EQUIBL_TYPE), INTENT(INOUT) :: EQUIBLT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE SETUP_POOLT
          END INTERFACE 
        END MODULE SETUP_POOLT__genmod
