        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:47 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIAGNOSTIC_WRITELU__genmod
          INTERFACE 
            SUBROUTINE DIAGNOSTIC_WRITELU(NPTS,NSAVEPEROUT,FILEID,      &
     &NUMVARS,NLEVS,VARIDS,WRITE_STEP,OUTSAVE)
              INTEGER(KIND=4), INTENT(IN) :: NLEVS
              INTEGER(KIND=4), INTENT(IN) :: NUMVARS
              INTEGER(KIND=4), INTENT(IN) :: NSAVEPEROUT
              INTEGER(KIND=4), INTENT(IN) :: NPTS
              INTEGER(KIND=4), INTENT(IN) :: FILEID
              INTEGER(KIND=4), INTENT(IN) :: VARIDS(NUMVARS)
              INTEGER(KIND=4), INTENT(IN) :: WRITE_STEP
              REAL(KIND=8), INTENT(IN) :: OUTSAVE(NPTS,NLEVS,NUMVARS,   &
     &NSAVEPEROUT)
            END SUBROUTINE DIAGNOSTIC_WRITELU
          END INTERFACE 
        END MODULE DIAGNOSTIC_WRITELU__genmod
