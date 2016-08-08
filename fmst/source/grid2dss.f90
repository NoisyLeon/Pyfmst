!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to generate a 3-D grid of regularly
! spaced velocity vertices in spherical coordinates. Velocity
! is allowed to vary with depth only. Velocity vertices are
! to be interpolated by cubic B-splines. The will also add, if
! required, random structure and/or a checkerboard.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM gridder
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: i,j,l,ars,checkstat
INTEGER :: nvt,nvp,aapmc,adch,nch,usp
INTEGER :: vusp1,vusp2,vusp1o,vusp2o
INTEGER :: asp,nsp,rseed
INTEGER, DIMENSION (100) :: ispt,ispp
REAL(KIND=i10) :: gotn,gots,gopw,gope,gst,gsp
REAL(KIND=i10) :: vback,vel,rssf,velp
REAL(KIND=i10) :: decm,mpc
REAL(KIND=i10) :: chvp,chvp1,chvp2
REAL(KIND=i10), DIMENSION (100) :: spa,spt,spp
REAL, EXTERNAL :: gasdev
CHARACTER (LEN=20) :: ofile,mfile
!
! nvt = number of vertices in theta (N-S) direction
! nvp = number of vertices in phi (E-W) direction
! gotn,gots = N-S bounds of grid
! gopw,gope = E-W bounds of grid
! gst = grid separation in theta (N-S)
! gsp = grid separation in phi (E-W)
! ofile = Name of output grid file
! mfile = Name of velocity model file
! vel = velocity of grid node
! vback = background velocity value
! ars = Add random structure (0=no,1=yes)
! rssf = Random structure scaling factor
! velp = Velocity perturbation obtained randomly
! aapmc = Add a priori model covariance (0=no, 1=yes)
! decm = Diagonal elements of covariance matrix
! adch = Add checkerboard (0=no,1=yes)
! mpc = Maximum perturbation of checkerboard
! nch = size of checkerboard cell
! chvp = Checkerboard velocity perturbation
! chvp1,chvp2 = Dummy checkerboard variables
! usp = use spacing for checkerboard (0=no,1=yes)
! vusp1,vusp2 = Checkerboard spacing variables
! vusp1o,vusp2o = Previous values of above
! asp = Apply spike (0=no,1=yes)
! nsp = Number of spikes
! spa = Amplitude of spikes
! spr,spp = Coordinates of spikes
! ispt,ispp = Grid locations of spikes
! gasdev = Function for returning random noise
! rseed = Random seed for noise generation
! 
OPEN(UNIT=10,FILE='grid2dss.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ofile
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)nvt
READ(10,*)nvp
READ(10,*)gotn,gots
READ(10,*)gopw,gope
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)vback
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ars
READ(10,*)rssf
READ(10,*)rseed
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)aapmc
READ(10,*)decm
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)adch
READ(10,*)mpc
READ(10,*)nch
READ(10,*)usp
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)asp
READ(10,*)nsp
DO i=1,nsp
   READ(10,*)spa(i)
   READ(10,*)spt(i),spp(i)
ENDDO
CLOSE(10)
!
! Compute gst,gsp
!
gst=(gotn-gots)/REAL(nvt-1)
gsp=(gope-gopw)/REAL(nvp-1)
!
! If spikes are required, compute grid locations
!
IF(asp.EQ.1)THEN
   DO i=1,nsp
      ispt(i)=(gotn-spt(i))/gst+1
      ispp(i)=(spp(i)-gopw)/gsp+1
   ENDDO
ENDIF
1 FORMAT(a20)
!
! Now write output to file, noting that
! we automatically add a cushion of
! boundary nodes to the grid.
!
OPEN(UNIT=20,FILE=ofile,STATUS='unknown')
WRITE(20,*)nvt,nvp
WRITE(20,'(3f14.8)')gotn,gopw
WRITE(20,'(3f14.8)')gst,gsp
WRITE(20,'(1X)')
IF(adch.EQ.1)THEN
   chvp1=mpc
   chvp2=mpc
   chvp=mpc
   vusp1=-1
   vusp2=-1
ENDIF
DO i=0,nvp+1
   IF(adch.EQ.1)THEN
      IF(MOD(i,nch).EQ.0)THEN
         chvp1=-chvp1
         IF(usp.EQ.1)THEN
            IF(vusp1.EQ.0)THEN
               IF(vusp1o.EQ.-1)THEN
                  vusp1=1
               ELSE
                  vusp1=-1
               ENDIF
            ELSE
               vusp1o=vusp1
               vusp1=0
            ENDIF
         ENDIF
      ENDIF
      chvp2=chvp1
      vusp2=1
      vusp2o=1
   ENDIF
   DO j=0,nvt+1
      IF(adch.EQ.1)THEN
         IF(MOD(j,nch).EQ.0)THEN
            chvp2=-chvp2
            IF(usp.EQ.1)THEN
               IF(vusp2.EQ.0)THEN
                  IF(vusp2o.EQ.-1)THEN
                     vusp2=1
                  ELSE
                     vusp2=-1
                  ENDIF
               ELSE
                  vusp2o=vusp2
                  vusp2=0
               ENDIF
            ENDIF
         ENDIF
         chvp=chvp2
      ENDIF
      vel=vback
!
!     Add random structure if required.
!
      IF(ars.EQ.1)THEN
         vel=vel+gasdev(rseed)*rssf
      ENDIF
!
!     Add checkerboard if required
!
      IF(adch.EQ.1)THEN
         vel=vel+vusp1*vusp2*chvp
      ENDIF
!
!     Apply spikes if required
!
      IF(asp.EQ.1)THEN
         DO l=1,nsp
            IF(i.EQ.ispp(l).OR.i.EQ.ispp(l)+1)THEN
               IF(j.EQ.ispt(l).OR.j.EQ.ispt(l)+1)THEN
                  vel=vel+spa(l)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!     Write out specified covariance value
!     if required
!
      IF(aapmc.EQ.1)THEN
         WRITE(20,'(2f12.8)')vel,decm
      ELSE
         WRITE(20,'(f12.8)')vel
      ENDIF
   ENDDO
   WRITE(20,'(1X)')
ENDDO
CLOSE(20)
STOP
END PROGRAM gridder

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
