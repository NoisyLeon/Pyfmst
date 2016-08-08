!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program constructs a synthetic set of traveltime
! residuals, with the option of adding noise.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM syntht
IMPLICIT NONE
INTEGER :: i
INTEGER :: agn,nr,ns,nt,checkstat,rseed
INTEGER, DIMENSION(:), ALLOCATABLE :: tswitch
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: sdgn,dsd,ot,rt,rsum,tpert
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: tres
REAL, EXTERNAL :: gasdev
CHARACTER (LEN=27) :: mtimes,otimes,sources,receivers
!
! mtimes = input file of synthetic observed times
! otimes = output file of traveltime residuals
! sources = file containing source information
! receivers = file containing receiver information
! agn = add Gaussian noise? (0=no,1=yes)
! sdgn = standard deviation of Gaussian noise
! dsd = default standard deviation for data covariance
! ns= number of sources
! nr = number of receivers
! nt = number of traveltimes
! ot = observed time
! rt = reference time
! tres = set of traveltime residuals
! tswitch = Switch for turning off rays
! tpert = Random value between -1 and 1
! rseed = Random seed for Gaussian noise generation
! gasdev = Gaussian noise value
!
! Read in the input parameters
!
ALLOCATE(tswitch(nt), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM syntht: REAL tswitch'
ENDIF
tswitch=1
OPEN(UNIT=10,FILE='synthtss.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a27)')mtimes
READ(10,'(a27)')sources
READ(10,'(a27)')receivers
READ(10,'(a27)')otimes
READ(10,*)agn
READ(10,*)sdgn
READ(10,*)rseed
READ(10,*)dsd
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
ALLOCATE(tres(nt), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM syntht: REAL tres'
ENDIF
!
! Read in the traveltimes and take
! the difference.
!
OPEN(UNIT=10,FILE=mtimes,STATUS='old')
rsum=0.0
DO i=1,nt
   READ(10,*) tswitch(i),tres(i)
ENDDO
CLOSE(10)
!
! Add Gaussian 
!
IF(agn.EQ.1)THEN
   dsd=sdgn
   DO i=1,nt
      tres(i)=tres(i)+gasdev(rseed)*sdgn
   ENDDO
ENDIF
!
! Write output to file
!
OPEN(UNIT=10,FILE=otimes,STATUS='unknown')
DO i=1,nt
   WRITE(10,1)tswitch(i),tres(i),dsd
ENDDO
1 FORMAT(i4,f10.5,f8.4)
CLOSE(10)
DEALLOCATE(tres, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM syntht: REAL tres'
ENDIF
END PROGRAM syntht

REAL FUNCTION gasdev(idum)
IMPLICIT NONE
INTEGER :: iset,i
INTEGER, PARAMETER :: imax=100000
INTEGER :: idum
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: fac,rsq,v1,v2
REAL(KIND=i10), SAVE :: gset
REAL, EXTERNAL :: ran1
iset=0
IF(iset.EQ.0)then
   DO i=1,imax
      v1=2.0*ran1(idum)-1.
      v2=2.0*ran1(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.LT.1.AND.rsq.NE.0.)EXIT
   ENDDO
   fac=sqrt(-2.0*LOG(rsq)/rsq)
   gset=v1*fac
   gasdev=v2*fac
   iset=1
ELSE
   gasdev=gset
   iset=0
ENDIF
END FUNCTION gasdev

REAL FUNCTION ran1(idum)
IMPLICIT NONE
INTEGER :: idum
INTEGER, PARAMETER :: ia=16807,im=2147483647,iq=127773
INTEGER, PARAMETER :: ir=2836,ntab=32,ndiv=1+(im-1)/ntab
INTEGER :: j,k
INTEGER, SAVE :: iy
INTEGER, DIMENSION (:), SAVE ::iv(ntab)
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10), PARAMETER :: eps=1.2e-7,rnmx=1.0-eps,am=1./im
iv=ntab*0
iy=0
IF(idum.LE.0.OR.iy.EQ.0)THEN
   DO j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      IF(idum.LT.0)idum=idum+im
      IF(j.LE.ntab)iv(j)=idum
   ENDDO
   iy=iv(1)
ENDIF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF(idum.LT.0)idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=min(am*iy,rnmx)
END FUNCTION ran1
