!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to calculate an estimate of model 
! roughness and variance by dicing up a cubic B-spline grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM misfit
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1
INTEGER :: nvt,nvp,cont,conp
INTEGER :: ndt,ndp,nnt,nnp,stt,stp
INTEGER :: checkstat
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: got,gop,gst,gsp,u,mvar,mrough
REAL(KIND=i10) :: rgst,rgsp,rdm,sumi,sumj,sumk,earth
REAL(KIND=i10) :: rdm1,sumi1,sumj1,sumk1,dp,dt,ri,risti
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: vi,wi
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: velsm,velrm,velr
REAL(KIND=i10), PARAMETER :: pi=3.14159265359
CHARACTER (LEN=25) :: rmfile,smfile
!
! nvt = number of vertices in theta (N-S) direction
! nvp = number of vertices in phi (E-W) direction
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! gst = grid separation in theta (N-S)
! gsp = grid separation in phi (E-W)
! smfile = solution model file 
! rmfile = reference model file
! velr = velocity at refined grid point
! velsm = velocity at spline vertex for solution model
! velrm = velocity at spline vertex for reference model
! ndt,ndp = node dicing level in r,theta,phi
! nnt,nnp = number of diced nodes in r,theta,phi
! u = Cubic spline independent variable
! vi,wi = Cubic spline basis functions
! rgst,rgsp = Refined node spacing in r,theta,phi
! sumi,sumj,sumk = Summation variables for constructing spline
! cont,conp = Counters for refining grid
! stt,stp = Refined grid location in r,theta,phi
! mvar = Model variance
! mrough = Model roughness
! sumi1,sumj1,sumk1= Summation variables for reference spline
! dp,dt = Contributions to roughness Laplacian
! ri,risti = Denomenators for difference operator
! earth = Earth radius
!
OPEN(UNIT=10,FILE='misfitss.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a25)')smfile
READ(10,'(a25)')rmfile
READ(10,*)ndt,ndp
READ(10,*)earth
CLOSE(10)
!
! Read in B-spline grid of solution model
!
OPEN(UNIT=20,FILE=smfile,STATUS='old')
READ(20,*)nvt,nvp
READ(20,*)got,gop
READ(20,*)gst,gsp
ALLOCATE(velsm(0:nvt+1,0:nvp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velv'
ENDIF
DO i=0,nvp+1
   DO j=0,nvt+1
      READ(20,*)velsm(j,i)
   ENDDO
ENDDO
CLOSE(20)
!
! Read in B-spline grid of reference model
!
OPEN(UNIT=20,FILE=rmfile,STATUS='old')
READ(20,*)nvt,nvp
READ(20,*)got,gop
READ(20,*)gst,gsp
got=(90.0-got)*pi/180.0
gop=gop*pi/180.0
gst=gst*pi/180.0
gsp=gsp*pi/180.0
ALLOCATE(velrm(0:nvt+1,0:nvp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velv'
ENDIF
DO i=0,nvp+1
   DO j=0,nvt+1
         READ(20,*)velrm(j,i)
   ENDDO
ENDDO
CLOSE(20)
! Calculate total numer of refined nodes in r,theta,phi
! and the refined grid spacing.
!
nnt=(nvt-1)*ndt+1
nnp=(nvp-1)*ndp+1
rgst=gst/ndt
rgsp=gsp/ndp
ALLOCATE(velr(nnt,nnp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velr'
ENDIF
!
! Calculate the values of the basis functions
!
ALLOCATE(vi(nvt+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL vi'
ENDIF
DO i=1,ndt+1
   u=ndt
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
ALLOCATE(wi(nvp+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL wi'
ENDIF
DO i=1,ndp+1
   u=ndp
   u=(i-1)/u
   wi(i,1)=(1.0-u)**3/6.0
   wi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   wi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   wi(i,4)=u**3/6.0
ENDDO
!
! Calculate velocity values on refined grid
!
mvar=0.0
DO i=1,nvp-1
   conp=ndp
   IF(i==nvp-1)conp=ndp+1
   DO j=1,nvt-1
      cont=ndt
      IF(j==nvt-1)cont=ndt+1
      DO l=1,conp
         stp=ndp*(i-1)+l
         DO m=1,cont
            stt=ndt*(j-1)+m
            sumi=0.0
            sumi1=0.0
            DO i1=1,4
               sumj=0.0
               sumj1=0.0
               DO j1=1,4
                  sumj=sumj+vi(m,j1)*velsm(j-2+j1,i-2+i1)
                  sumj1=sumj1+vi(m,j1)*velrm(j-2+j1,i-2+i1)
               ENDDO
               sumi=sumi+wi(l,i1)*sumj
               sumi1=sumi1+wi(l,i1)*sumj1
            ENDDO
            velr(stt,stp)=sumi-sumi1
            mvar=mvar+(sumi-sumi1)**2
         ENDDO
      ENDDO
   ENDDO
ENDDO
mvar=mvar/(nnp*nnt-1.0)
!
! Now estimate model roughness
!
mrough=0.0
DO i=2,nnp-1
   DO j=2,nnt-1
      ri=earth
      risti=ri*sin(got+(j-1)*rgst)
      dp=velr(j,i+1)-2.0*velr(j,i)+velr(j,i-1)
      dp=dp/((risti*rgsp)**2)
      dt=velr(j+1,i)-2.0*velr(j,i)+velr(j-1,i)
      dt=dt/((ri*rgst)**2)
      mrough=mrough+ABS(dp)+ABS(dt)
   ENDDO
ENDDO
mrough=mrough/((nnp-2)*(nnt-2))
WRITE(6,*)'Model variance in (km/s)**2 is ',mvar
WRITE(6,*)'Model roughness in (kms)**(-1) is ',mrough
DEALLOCATE(vi,wi, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM cbsvel: REAL vi,wi'
ENDIF
DEALLOCATE(velsm,velrm,velr, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM cbsvel: REAL velsm,velrm,velr'
ENDIF
STOP
END PROGRAM misfit
