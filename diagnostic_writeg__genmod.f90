        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:47 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIAGNOSTIC_WRITEG__genmod
          INTERFACE 
            SUBROUTINE DIAGNOSTIC_WRITEG(NPTS,NTPERSAVE,NSAVEPEROUT,    &
     &NUMVARS,FILEID,WRITE_STEP,VARREFS,VARIDS,SIF_SATCOUNT,OUTSAVE)
              INTEGER(KIND=4), INTENT(IN) :: NUMVARS
              INTEGER(KIND=4), INTENT(IN) :: NSAVEPEROUT
              INTEGER(KIND=4), INTENT(IN) :: NPTS
              INTEGER(KIND=4), INTENT(IN) :: NTPERSAVE
              INTEGER(KIND=4), INTENT(IN) :: FILEID
              INTEGER(KIND=4), INTENT(IN) :: WRITE_STEP
              INTEGER(KIND=4), INTENT(IN) :: VARREFS(NUMVARS)
              INTEGER(KIND=4), INTENT(IN) :: VARIDS(NUMVARS)
              INTEGER(KIND=4), INTENT(IN) :: SIF_SATCOUNT(NPTS,         &
     &NSAVEPEROUT,2)
              REAL(KIND=8), INTENT(IN) :: OUTSAVE(NPTS,NUMVARS,         &
     &NSAVEPEROUT)
            END SUBROUTINE DIAGNOSTIC_WRITEG
          END INTERFACE 
        END MODULE DIAGNOSTIC_WRITEG__genmod
