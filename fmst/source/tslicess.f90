!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to convert velocity, traveltime
! and raypath information from the 2-D spherical shell tomography
! program into a form suitable for input to GMT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM slice
IMPLICIT NONE
INTEGER :: i,j,l,m,i1,j1
INTEGER :: checkstat,isw,conp,cont
INTEGER :: nnt,nnp,pvg,prp,nre,nr,parv
INTEGER :: idm2,idm3,nnx,nnz,stp,stt
INTEGER :: ddt,ddp,ptf,pls,plr,ns,nrc
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i5) :: rlat,rlon
REAL(KIND=i10) :: got,gop,rgst,rgsp
REAL(KIND=i10) :: sumi,sumj
REAL(KIND=i10) :: lft,rgt,btm,top
REAL(KIND=i10) :: rdm,rd1,rd2,rd3,u
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,vela
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ttn,veln
REAL(KIND=i5), PARAMETER :: tol=1.0e-5
CHARACTER (LEN=30) :: ifilet,ifilev,irfile,ifilevr,sep
CHARACTER (LEN=30) :: ofileb,ofilev,ofilet,ofiler
CHARACTER (LEN=30) :: ifiles,ofiles,ifilerc,ofilerc
!
! ifilet = input 2-D traveltime grid file
! ifilev = input velocity grid file
! ifilevr = input reference velocity grid file
! irfile = input ray path file
! ofilet = output file for traveltimes
! ofileb = bounds for slice file
! ofilev = output velocity file
! ofiler = output ray path file
! nnt,nnp = number of diced nodes in r,theta,phi
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! rgst,rgsp = Refined node spacing in r,theta,phi
! ttn = traveltime grid values
! lft,rgt,btm,top = plotting bounds
! pvg = plot velocity grid (0=no,1=yes)
! parv = plot absolute (0) or relative (1) velocity
! ptf = plot traveltime field depth slice? (0=no,1=yes)
! veln = velocity grid values
! prp = plot ray paths? (0=no,1=yes)
! nre = number of ray elements
! rlat,rlon = ray point depth, latitude, longitude
! nr = number of receivers
! sep = character marker for separating rays
! ddt,ddp = dicing of velocity grid
! u = bspline independent coordinate
! ui,vi = bspline basis functions
! vela = diced velocity values
! nnx,nnz = dimensions of vela
! conp,cont = variables for edge of bspline grid
! stp,stt = counters for vela grid points
! pls = plot sources (0=no,1=yes)
! plr = plot receivers (0=no,1=yes)
! ifiles = input source file
! ofiles = output source file
! ifilerc = input receiver file
! ofilerc = output receiver file
! ns = number of sources
! nrc = number of receivers
!
OPEN(UNIT=10,FILE='tslicess.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a30)')ifilev
READ(10,'(a30)')ifilevr
READ(10,'(a30)')ifilet
READ(10,'(a30)')irfile
!
! Bounding box of plot
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a22)')ofileb
!
! Source and receiver information
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)pls
READ(10,'(a30)')ifiles
READ(10,'(a30)')ofiles
READ(10,*)plr
READ(10,'(a30)')ifilerc
READ(10,'(a30)')ofilerc
!
! Output source and receiver information if required
!
IF(pls.EQ.1)THEN
   OPEN(UNIT=20,FILE=ifiles,STATUS='old')
   OPEN(UNIT=30,FILE=ofiles,STATUS='unknown')
   READ(20,*)ns 
   DO i=1,ns
      READ(20,*)rlat,rlon
      WRITE(30,*)rlat,rlon
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
IF(plr.EQ.1)THEN
   OPEN(UNIT=20,FILE=ifilerc,STATUS='old')
   OPEN(UNIT=30,FILE=ofilerc,STATUS='unknown')
   READ(20,*)nrc
   DO i=1,nrc
      READ(20,*)rlat,rlon
      WRITE(30,*)rlat,rlon
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
!
! Calculate GMT bounds files. Start off by reading in
! velocity grid.
!
OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)nnt,nnp
READ(20,*)got,gop
READ(20,*)rgst,rgsp
ALLOCATE(veln(0:nnt+1,0:nnp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL veln'
ENDIF
DO i=0,nnp+1
   DO j=0,nnt+1
      READ(20,*)veln(j,i)
   ENDDO
ENDDO
CLOSE(20)
!
! Now read in velocity grid parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)pvg
READ(10,*)parv
READ(10,*)ddt,ddp
READ(10,'(a22)')ofilev
!
! Calculate GMT bounds file for depth slice if required.
!
lft=gop
rgt=gop+(nnp-1)*rgsp
btm=got-(nnt-1)*rgst
top=got
OPEN(UNIT=50,FILE=ofileb,STATUS='unknown')
WRITE(50,'(f16.10)')lft
WRITE(50,'(f16.10)')rgt
WRITE(50,'(f16.10)')btm
WRITE(50,'(f16.10)')top
rd1=rgsp/ddp
rd2=rgst/ddt
WRITE(50,'(f16.10)')rd1
WRITE(50,'(f16.10)')rd2
!
! Read in reference velocity grid file if required
!
IF(pvg.EQ.1.AND.parv.EQ.1)THEN
   isw=0
   OPEN(UNIT=20,FILE=ifilevr,status='old')
   READ(20,*)idm2,idm3
   IF(idm2.NE.nnt.OR.idm3.NE.nnp)isw=1
   READ(20,*)rd2,rd3
   IF(ABS(rd2-got).GT.tol.OR.ABS(rd3-gop).GT.tol)isw=1
   READ(20,*)rd2,rd3
   IF(ABS(rd2-rgst).GT.tol.OR.ABS(rd3-rgsp).GT.tol)isw=1
   IF(isw.EQ.1)THEN
      WRITE(6,*)'ERROR! Actual velocity grid and reference'
      WRITE(6,*)'velocity grid have different dimensions or'
      WRITE(6,*)'different numbers of grid points!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
   DO i=0,nnp+1
      DO j=0,nnt+1
         READ(20,*)rd1
!
!        This gives the velocity anomaly.
!
         veln(j,i)=veln(j,i)-rd1
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Extract velocity slice if required
!
IF(pvg.EQ.1)THEN
!
! Allocate memory to velocity grid array
!
  nnx=(nnt-1)*ddt+1
  nnz=(nnp-1)*ddp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(ddt+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,ddt+1
     u=ddt
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(ddp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,ddp+1
     u=ddp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
  DO i=1,nnp-1
     conp=ddp
     IF(i==nnp-1)conp=ddp+1
     DO j=1,nnt-1
        cont=ddt
        IF(j==nnt-1)cont=ddt+1
        DO l=1,conp
           stp=ddp*(i-1)+l
           DO m=1,cont
              stt=ddt*(j-1)+m
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumj=sumj+ui(m,j1)*veln(j-2+j1,i-2+i1)
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(stt,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofilev,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
DEALLOCATE(veln, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
ENDIF
!
! Read in input parameters for traveltime grid
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ptf
READ(10,'(a22)')ofilet
!
! Open traveltime field file if required
!
IF(ptf.EQ.1)THEN
   OPEN(UNIT=20,FILE=ifilet,FORM='unformatted',STATUS='old')
   READ(20)got,gop
   READ(20)nnt,nnp
   READ(20)rgst,rgsp
   ALLOCATE(ttn(nnp,nnt), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL ttn'
   ENDIF
   DO i=1,nnp
      DO j=1,nnt
         READ(20)ttn(i,j)
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Finish off GMT boundary files
!
WRITE(50,'(f16.10)')rgsp
WRITE(50,'(f16.10)')rgst
CLOSE(50)
!
! Extract slice if required.
!
IF(ptf.EQ.1)THEN
!
!  Write out slice to GMT xy file
!
   OPEN(UNIT=30,FILE=ofilet,STATUS='unknown')
   DO i=1,nnp
      DO j=nnt,1,-1
         WRITE(30,*)ttn(i,j)
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
!
! Read in ray path parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)prp
READ(10,'(a30)')ofiler
CLOSE(10)
!
! Plot raypaths if required.
!
IF(prp.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofiler,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      IF(nre.NE.0)THEN
         DO j=1,nre
            READ(20)rlat,rlon
            WRITE(30,'(2f10.4)')rlon,rlat
         ENDDO
         WRITE(30,'(a1)')sep
      ENDIF
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
STOP
END PROGRAM slice
