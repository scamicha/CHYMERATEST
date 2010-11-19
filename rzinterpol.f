      program inter

C      implicit none.  This program reads the saved file to make
C      the same file with quadruple azimuthal resolution
c      Compile f77 -o interpol interpol.f
      
      IMPLICIT NONE

      Integer  jmax,jmax1,jmax2,newjmax2
      integer  lmax,newlmax
      integer  kmax,kmax1,kmax2,newkmax2
      parameter (jmax=512, jmax1=jmax+1, jmax2=jmax+2,newjmax2=1026)
      parameter (kmax=64, kmax1=kmax+1, kmax2=kmax+2,newkmax2=130)
      parameter (lmax=2048,newlmax=2048)
      integer jreq, j, k, l, lp, newjreq
      real*8 rof3n,zof3n,newrof3n,newzof3n,area,pi
      real*8 delt,den,time,elost,sound,omcen,tmassini,tmass
      real*8 tmassadd,tmassout,tmassacc,totcool,totheat
      real*8 etot,newetot,rhotot,stot,atot,ttot,etoterr
      REAL*8 etotfl,eflufftot,totirr,totdflux
      real*8 r(jmax2),rhf(jmax2),newr(newjmax2),newrhf(newjmax2)
      real*8 rho(jmax2,kmax2,lmax),s(jmax2,kmax2,lmax),t(jmax2,kmax2
     &     ,lmax),a(jmax2,kmax2,lmax),eps(jmax2,kmax2,lmax)
      real*8 newrrho(newjmax2,kmax2,newlmax),newrs(newjmax2,kmax2
     &     ,newlmax),newrt(newjmax2,kmax2,newlmax),newra(newjmax2
     &     ,kmax2,newlmax),newreps(newjmax2,kmax2,newlmax)
      real*8 newrho(newjmax2,newkmax2,newlmax),news(newjmax2,newkmax2
     &     ,newlmax),newt(newjmax2,newkmax2,newlmax),newa(newjmax2
     &     ,newkmax2,newlmax),neweps(newjmax2,newkmax2,newlmax)

      parameter (pi=3.14159265358979d0)

      open(unit=11,file='savedLMAX2048.00550000',form='unformatted')
      print*,'Opening file'
      read(11)s
      read(11)t
      read(11)a
      read(11)rho
      read(11)eps
      read(11)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMCEN
      read(11)tmassini,tmass,tmassadd,tmassout,tmassacc,totcool
     &     ,totdflux,totheat,totirr,etotfl,eflufftot
      print*,'Check format',rho(208,3,4),rho(235,3,2),time
      close(11)

      newrof3n = rof3n/2.d0
      newzof3n = zof3n/2.d0

      do j=1,jmax2
         r(j) = (j-2)*rof3n
         rhf(j) = r(j)+rof3n/2.d0
      ENDDO

      do j=1,newjmax2
         newr(j) = (j-2)*newrof3n
         newrhf(j) = newr(j)+newrof3n/2.d0
      ENDDO      

      
      do j=1,jmax1
         do k=1,kmax2
            do l=1,lmax

C              For rho
               
               newrrho(2*j-1,k,l) = (3.d0/4.d0)*rho(j,k,l)+(1.d0/4.d0)*
     &            rho(j+1,k,l)
               newrrho(2*j,k,l) = (1.d0/4.d0)*rho(j,k,l)+(3.d0/4.d0)*
     &            rho(j+1,k,l)
 

C              For S

               newrs(2*j-1,k,l) = (3.d0/4.d0)*s(j,k,l)+(1.d0/4.d0)*
     &            s(j+1,k,l)
               
               newrs(2*j,k,l) = (1.0/4.0)*s(j,k,l)+(3.0/4.0)*s(j+1,k,l)
               
               
C              For T
               
               newrt(2*j-1,k,l) = (3.d0/4.d0)*t(j,k,l)+(1.d0/4.d0)*
     &            t(j+1,k,l)
               
               newrt(2*j,k,l) = (1.d0/4.d0)*t(j,k,l)+(3.d0/4.d0)*
     &            t(j+1,k,l)
               
               
C              For A
                
               newra(2*j-1,k,l) = (3.d0/4.d0)*a(j,k,l)+(1.d0/4.d0)*
     &            a(j+1,k,l)
               
               newra(2*j,k,l) = (1.d0/4.d0)*a(j,k,l)+(3.d0/4.d0)*
     &            a(j+1,k,l)
               
               
C              For EPS
               
               newreps(2*j-1,k,l) = (3.d0/4.d0)*eps(j,k,l)+(1.d0/4.d0)*
     &            eps(j+1,k,l)
               
               newreps(2*j,k,l) = (1.d0/4.d0)*eps(j,k,l)+(3.d0/4.d0)*
     &            eps(j+1,k,l)
            enddo
         enddo
      enddo




      do j=1,newjmax2
         do k=1,kmax1
            do l=1,lmax

C              For rho
               
               newrho(j,2*k-1,l) = (3.d0/4.d0)*newrrho(j,k,l)+
     &            (1.d0/4.d0)*newrrho(j,k+1,l)
               newrho(j,2*k,l) = (1.d0/4.d0)*newrrho(j,k,l)+
     &            (3.d0/4.d0)*newrrho(j,k+1,l)
 

C              For S

               news(j,(2*k)-1,l) = (3.d0/4.d0)*newrs(j,k,l)+(1.d0/4.d0)*
     &            newrs(j,k+1,l)
               
               news(j,2*k,l) = (1.0/4.0)*newrs(j,k,l)+(3.0/4.0)*
     &            newrs(j,k+1,l)
               
               
C              For T
               
               newt(j,2*k-1,l) = (3.d0/4.d0)*newrt(j,k,l)+(1.d0/4.d0)*
     &            newrt(j,k+1,l)
               
               newt(j,2*k,l) = (1.d0/4.d0)*newrt(j,k,l)+(3.d0/4.d0)*
     &            newrt(j,k+1,l)
               
               
C              For A
                
               newa(j,2*k-1,l) = (3.d0/4.d0)*newra(j,k,l)+(1.d0/4.d0)*
     &            newra(j,k+1,l)
               
               newa(j,2*k,l) = (1.d0/4.d0)*newra(j,k,l)+(3.d0/4.d0)*
     &            newra(j,k+1,l)
               
               
C              For EPS
               
               neweps(j,2*k-1,l) = (3.d0/4.d0)*newreps(j,k,l)+
     &            (1.d0/4.d0)*newreps(j,k+1,l)
               
               neweps(j,2*k,l) = (1.d0/4.d0)*newreps(j,k,l)+(3.d0/4.d0)*
     &            newreps(j,k+1,l)

            enddo
         enddo
      enddo    


      open(unit=10, file='etot.dat',form='formatted')
      
      etot = 0.d0
      newetot = 0.d0
      etoterr = 0.d0

c      do j=1,jmax
c         do k=1,kmax
c            do l=1,lmax
c               newetot = 0.25*(neweps(2*j,2*k,l)+neweps(2*j+1,2*k,l)+
c     &                   neweps(2*j,2*k+1,l)+neweps(2*j+1,2*k+1,l))
c               etot = eps(j+1,k+1,l)
c               etoterr = etoterr + abs((etot-newetot)/etot)
c               write(10,*)newetot,etot,etoterr
c            ENDDO
c         ENDDO
c      ENDDO

c      close(10)

      do j=2,jmax1
         area=pi*(rhf(j+1)**2-rhf(j)**2)
         do k=2,kmax1
            do l=1,lmax
               etot = etot + rof3n*area*t(j,k,l)
            ENDDO
         ENDDO
      ENDDO

      do j=2,newjmax2-1
         area=pi*(newrhf(j+1)**2-newrhf(j)**2)
         do k=2,newkmax2-1
            do l=1,lmax
               newetot = newetot + newrof3n*area*newt(j,k,l)
            ENDDO
         ENDDO
      ENDDO
      

      etoterr = abs((etot-newetot)/etot)
      print*,etoterr

      do j=1,newjmax2
C         print*,j,newrho(j,5,1),news(j,5,1),newt(j,5,1),newa(j,5,1)
C     &        ,neweps(j,5,1)
      enddo

!      print*,'rho',rho(100,7,lmax),rho(100,7,1)
!      print*,'new',newrho(100,7,509),newrho(100,7,510),newrho(100,7,511)
!     &     ,newrho(100,7,512)
      

!      print*,'a',a(100,7,lmax),a(100,7,1)
!      print*,'new',newa(100,7,509),newa(100,7,510),newa(100,7,511)
!     &     ,newa(100,7,512)

!      print*,'rho',rho(5,7,lmax),rho(5,7,1)
!      print*,'new',newrho(5,7,509),newrho(5,7,510),newrho(5,7,511)
!     &     ,newrho(5,7,512)
      

!      print*,'a',a(5,7,lmax),a(5,7,1)
!      print*,'new',newa(5,7,509),newa(5,7,510),newa(5,7,511)
!     &     ,newa(5,7,512)

C      print*,'etoterr = ', etoterr

      print*,newrof3n,newzof3n,tmassini

      newjreq = 2*(jreq-1)+1

      print*, "newjreq = ", newjreq

      open(unit=8,file='highrzsavedLMAX2048.00550000',
     &     form='unformatted')
      WRITE(8) newS
      WRITE(8) newT
      WRITE(8) newA
      WRITE(8) newRHO
      WRITE(8) newEPS
      Write(8) newROF3N,newZOF3N,DELT,TIME,ELOST,DEN,SOUND,newJREQ,OMCEN
      write(8) tmassini,tmass,tmassadd,tmassout,tmassacc,totcool,
     &     totdflux,totheat,totirr,etotfl,eflufftot
      close(8)


      end
