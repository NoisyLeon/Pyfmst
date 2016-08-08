!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program outputs traveltime residual data in a form
! suitable for input to GMT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM resplot
IMPLICIT NONE
INTEGER :: i,j,idm,idm1,idm2,iorf
INTEGER :: nr,ns,nt,checkstat
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: otr,rt,mt,tres,prs,rd1,rd2
CHARACTER (LEN=40) :: mtimes,rtimes,otimes,sources,receivers,ofile
!
! rtimes = input file of reference times
! otimes = output file of traveltime residuals
! sources = file containing source information
! receivers = file containing receiver information
! mtimes = file containing model traveltimes
! ofile = output file for plotting
! ns= number of sources
! nr = number of receivers
! nt = number of traveltimes
! otr = observed traveltime residual
! rt = reference time
! mt = model traveltime
! iorf = initial (0) or final (1) residuals plotted
! prs = print residuals larger than this size
!
! Read in the input parameters
!
OPEN(UNIT=10,FILE='resplotss.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)iorf
READ(10,'(a32)')ofile
READ(10,*)prs
READ(10,'(a32)')otimes
READ(10,'(a32)')rtimes
READ(10,'(a32)')sources
READ(10,'(a32)')receivers
READ(10,'(a32)')mtimes
CLOSE(10)
!
! Determine the number of traveltimes
!
OPEN(UNIT=10,FILE=sources,STATUS='old')
READ(10,*)ns
CLOSE(10)
OPEN(UNIT=10,FILE=receivers,STATUS='old')
READ(10,*)nr
DO i=1,nr
   read(10,*)rd1,rd2
ENDDO
CLOSE(10)
nt=ns*nr
!
! Read in the model and reference traveltimes and take
! the difference and then subtract the observed residuals.
!
OPEN(UNIT=10,FILE=mtimes,STATUS='old')
OPEN(UNIT=20,FILE=rtimes,STATUS='old')
OPEN(UNIT=30,FILE=otimes,STATUS='old')
OPEN(UNIT=40,FILE=ofile,STATUS='unknown')
DO i=1,ns
   DO j=1,nr
      READ(10,*)idm1,mt
      READ(20,*)idm2,rt
      READ(30,*)idm,otr
      IF(idm.NE.idm1.OR.idm.NE.idm2)THEN
         WRITE(6,*)'ERROR!!! Input time files are'
         WRITE(6,*)'inconsistent!!!'
         WRITE(6,*)'TERMINATING PROGRAM!!!'
         STOP
      ENDIF
      IF(idm.EQ.1)THEN
         IF(iorf.EQ.0)THEN
            tres=otr-rt
         ELSE
            tres=otr-mt
         ENDIF
         WRITE(40,*)tres
         IF(ABS(tres).GT.prs)THEN
            WRITE(6,1)j,i,tres
1           FORMAT('Recorder ',i4,' of source ',i4,' has a &
           &residual size of ',f9.4,' s.')
         ENDIF
      ENDIF
   ENDDO
ENDDO
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
END PROGRAM resplot
