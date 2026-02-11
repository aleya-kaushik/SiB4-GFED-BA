        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_ROOTF__genmod
          INTERFACE 
            SUBROUTINE SETUP_ROOTF(KROOT,ROOTD,DZ,LAYER_Z,NODE_Z,ROOTF)
              REAL(KIND=4), INTENT(IN) :: KROOT
              REAL(KIND=4), INTENT(IN) :: ROOTD
              REAL(KIND=8), INTENT(IN) :: DZ(10)
              REAL(KIND=8), INTENT(IN) :: LAYER_Z(10)
              REAL(KIND=8), INTENT(IN) :: NODE_Z(10)
              REAL(KIND=8), INTENT(INOUT) :: ROOTF(10)
            END SUBROUTINE SETUP_ROOTF
          END INTERFACE 
        END MODULE SETUP_ROOTF__genmod
