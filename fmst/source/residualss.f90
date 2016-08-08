!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program calculates the RMS traveltime residual and
! variance of the current model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM residual
IMPLICIT NONE
INTEGER :: i,idm1,idm2,isum
INTEGER :: nr,ns,nt,checkstat
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: otr,mt,mtr,rsum,var
CHARACTER (LEN=27) :: mtimes,rtimes,otimes,sources,receivers
!
! otimes = output file of traveltime residuals
! sources = file containing source information
! receivers = file containing receiver information
! mtimes = file containing model traveltimes
! ns= number of sources
! nr = number of receivers
! nt = number of traveltimes
! otr = observed traveltime residual
! mt = model traveltime
! mtr = set of traveltime residuals
! var = variance
! idm1,idm2 = switches that check ray validity
!
! Read in the input parameters
!
OPEN(UNIT=10,FILE='residualss.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a27)')otimes
READ(10,'(a27)')sources
READ(10,'(a27)')receivers
READ(10,'(a27)')mtimes
CLOSE(10)
!
! Determine the number of traveltimes
!
OPEN(UNIT=10,FILE=sources,STATUS='old')
READ(10,*)ns
CLOSE(10)
OPEN(UNIT=10,FILE=receivers,STATUS='old')
READ(10,*)nr
CLOSE(10)
nt=ns*nr
!
! Read in the model and observed traveltimes and take
! the difference.
!
OPEN(UNIT=10,FILE=mtimes,STATUS='old')
OPEN(UNIT=30,FILE=otimes,STATUS='old')
rsum=0.0
isum=0
DO i=1,nt
   READ(10,*)idm1,mt
   READ(30,*)idm2,otr
   IF(idm1.NE.idm2)THEN
      WRITE(6,*)'ERROR!!! Model and observed times'
      WRITE(6,*)'files are inconsistent!!!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
   IF(idm2.EQ.1)THEN
      mtr=mt-otr
      rsum=rsum+mtr**2
      isum=isum+1
   ENDIF
ENDDO
nt=isum
CLOSE(10)
CLOSE(30)
var=rsum/REAL(nt-1)
rsum=SQRT(rsum/REAL(nt))
rsum=1000.0*rsum
WRITE(6,'(f12.2,f12.5)')rsum,var
END PROGRAM residual
