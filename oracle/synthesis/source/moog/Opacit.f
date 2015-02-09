
      subroutine opacit (modeop,waveop)
c******************************************************************************
c     This routine calculates the continuous opacity *kaplam* or            
c     *kapref* and optical depth *taulam* as needed. the formulas are       
c     ones lifted from the write-up of program *atlas*                      
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      include 'Dummy.com'
      real*8 taucheck(100), tauinteg(100)


c*****initialization                                                        
      freq = 2.997925d18/waveop
      freqlg = dlog(freq)
      evhkt(:) = 0.0
      do i = 1,ntau
         hkt = 6.6256d-27/(1.38054d-16*t(i))
         evhkt(i) = dexp(-freq*hkt)
      enddo

c      do i=1,ntau
       aH1(:) =      1.d-99
       aHminus(:) =  1.d-99
       aHeminus(:) = 1.d-99
       aC1(:) =      1.d-99
       aMg1(:) =     1.d-99
       aMg2(:) =     1.d-99
       aAl1(:) =     1.d-99
       aSi1(:) =     1.d-99
       aSi2(:) =     1.d-99
       aFe1(:) =     1.d-99
c      enddo


c*****compute the opacities
       call opacH1
       call opacHminus
       call opacHscat
       call opacH2scat
       call opacHeminus
       call opacHescat
       call opacescat
       call opacC1
       call opacMg1
       call opacMg2
       call opacAl1
       call opacSi1
       call opacSi2
       call opacFe1


c*****sum up all the opacities
      do i=1,ntau
         kaplamabs(i) = aH1(i) + aHminus(i) + aHeminus(i) + aC1(i) + 
     .                  aMg1(i) + aMg2(i) + aAl1(i) + aSi1(i) +
     .                  aSi2(i) + aFe1(i)
         kaplamsca(i) = sigH(i) + sigH2(i) + sigHe(i) + sigel(i)
         kaplam(i)    = kaplamabs(i) + kaplamsca(i)
c         stop   
      enddo
c      print *, "OPACIT SUM", aH1(1), aHminus(1), aHeminus(1), 
c     .      aC1(1), aMg1(1), aMg2(1), aAl1(1), aSi1(1), aSi2(1),
c     .      aFe1(1), sigH(1), sigH2(1), sigHe(1), sigel(1)
c      print *, "MODEOP", modeop

c      print *, "kaplam", kaplam(60)
      

c*****write out the opacities
      if (debug .gt. 0) then
         write (*,1001) waveop
         do i=1,ntau
            write (*,1002) i, nint(t(i)), dlog10(kaplam(i)),
     .      dlog10(aH1(i)), dlog10(aHminus(i)), dlog10(sigH(i)), 
     .      dlog10(aHeminus(i)), dlog10(sigHe(i)), dlog10(sigel(i)),
     .      dlog10(sigH2(i))
         enddo
         write (*,1003) waveop
         do i=1,ntau
            write (*,1002) i, nint(t(i)), dlog10(kaplam(i)),
     .      dlog10(aC1(i)),  dlog10(aMg1(i)), dlog10(aMg2(i)), 
     .      dlog10(aAl1(i)), dlog10(aSi1(i)), dlog10(aSi2(i)), 
     .      dlog10(aFe1(i))
         enddo
      endif


c*****add in an arbitrary amount of "extra" continuous opacity;
c     this is a pure fudge factor, so experiment with care
      if (fudge .gt. 0.0) then
         do i=1,ntau
            kaplam(i) = kaplam(i)*((fudge*10000)/t(i))
            if (i .eq. 1) then
               write(*,1005)
               write(*,1006) fudge
            endif
         enddo
      endif


c*****compute an optical depth array at this wavelength, and exit
      if (modeop .ne. 1) then
         do  i=1,ntau
            dummy1(i) = tauref(i)*kaplam(i)/(0.4343*kapref(i))
         enddo
         first = tauref(1)*kaplam(1)/kapref(1)
         dummy1(1) = rinteg(xref,dummy1,taulam,ntau,first)
         taulam(1) = first
         do i=2,ntau
            taulam(i) = taulam(i-1) + taulam(i)
         enddo
         if (debug .gt. 0) write(*,1004) (taulam(i),i=1,ntau)
         return
      endif


c*****here is the assignment of opacities kaplam to kapref; exit normally
      do i=1,ntau
         kapref(i) = kaplam(i)
      enddo

c      print *, "KAPREF IS ", kapref
c      stop
      if(modprintopt .lt. 2) return


c*****here is an internal tauref check on the externally read in tau scale;
c     then exit normally
      do i=1,ntau
         dummy1(i) = dlog10(tauref(i))
         tauinteg(i) = tauref(i)*kaplam(i)/(0.4343*kapref(i))
      enddo
      first = tauref(1)
      tauinteg(1) = rinteg(dummy1,tauinteg,taucheck,ntau,first)
      taucheck(1) = first
      do i=2,ntau
         taucheck(i) = taucheck(i-1) + taucheck(i)
      enddo
      if (debug .gt. 0) write (*,1004) (taucheck(i),i=1,ntau)
      return


c*****format statements
1001  format (/'log opacities due to H, He, e- at a wavelength of ',
     .       f10.2,' A'/' i  T(i)  kaplam     aH1     aH-    sigH',
     .       '    aHe-   sigHe   sigel   sigH2')
1002  format (i2,i6,10f8.2)
1003  format (/'log opacities due to "metal" edges at a wavelength',
     .       ' of ',f10.2,' A'/' i  T(i)  kaplam     aC1    aMg1', 
     .       '    aMg2    aAl1    aSi1    aSi2    aFe1')
1004  format(/'COMPUTED TAULAM ARRAY'/(1p7e11.3))
1005  format ('Opacities artificially are now ',
     .        'altered by scaling law of: ')
1006  format ('      kaplam(i) = kaplam(i) * ',
     .        '((', f6.2, ' * 10000)/t(i))')


      end




