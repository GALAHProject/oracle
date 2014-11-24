
      subroutine synspec
c******************************************************************************
c     This routine does synthetic spectra                                
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      real*8 dd(5000)

c      print *, "forcing lineflat = 0"
c*****initialize the synthesis
      if (debug .gt. 0) then 
         write (nf1out,1101)
         write (nf2out,1002) moditle(1:73)
         if (iunits .eq. 1) then
            write (nf2out,1103) start/1.d4, sstop/1.d4,
     .                       step/1.d4, delta/1.d4
         else
            write (nf2out,1102) start,sstop,step,delta
         endif
      endif
c      if (iraf .eq. 1) then
c         npoint = (sstop-start)/step
c         write (nf4out,1104) npoint,wave,wave,step,step
c         write (nf4out,1105)
c         write (nf4out,1106) moditle
c         do j=1,93
c            if (pec(j) .gt. 0 ) then
c               dummy1(j) = dlog10(xabund(j)) + 12.0
c               write (nf4out,1107) names(j),dummy1(j)
c            endif
c         enddo
c         write (nf4out,1108) vturb(1)
c         write (nf4out,1109)
c      endif
      if (mode .ne. 4) then 
         lim1line = 0
         lim2line = 0
         lim1obs = 0
         lim2obs = 0
         lim1 = 0
         lim2 = 0
      endif


c*****now step in wavelength as the spectrum is computed
      kount = nint((sstop - start + (step/4.0) )/step) + 1
      num = 0
      wavl = 0.
      nkount = kount

c*****first calculate or recalculated continuum quantities at the 
c     spectrum wavelength, if needed
c      print *, "NUMBER OF POINTS", kount
      
      do n=1,kount
         num = num + 1
         wave = oldstart + (n-1)*step
         if (dabs(wave-wavl)/wave .ge. 0.001) then
            wavl = wave   
            call opacit (2,wave)    
c            if (debug .ge. 0) 
c     .          write (nf1out,1001) wave,(kaplam(i),i=1,ntau)
            call cdcalc (1)  
            first = 0.4343*cd(1)
            flux = rinteg(xref,cd,dummy1,ntau,first)
c            if (iunits .eq. 1) then
c               write (nf1out,1003) 1.d-4*wave,flux
c            else
c               write (nf1out,1004) wave,flux
c            endif
         endif


c*****find the appropriate set of lines for this wavelength, reading 
c     in a new set if this is the initial depth calculation or if
c     needed because the line list end has been reached
         if (mode .eq. 3) then
20          call linlimit
            if (lim2line .lt. 0) then
               call inlines (2)
               call nearly (1)
               go to 20
            endif
            lim1 = lim1line
            lim2 = lim2line
         endif


c*****compute a spectrum depth at this point; if there are no absorption
c     lines in the interval then just set the depth to zero without
c     extensive line calculations
         lineflag = 0
c         print *, "forcing lineflag = 0"
c         print *, "lineflag", lineflag, n
c         print *, "num", mode, wave, num, lim1, lim2, lim1line, lim2line
         if (lineflag .lt. 0) then
            d(num) = 0.
         else
            call taukap   
            call cdcalc (2)
            first = 0.4343*cd(1)
            d(num) = rinteg(xref,cd,dummy1,ntau,first)
            computed_wls(num) = wave

            if (d(num) .lt. 0.0) then
               computed_fluxes(num) = 1.
            else
               computed_fluxes(num) = 1. - d(num)
            endif
         endif

c        Update the synthesis array
c         print *, "OK", num, d(num-10:num)
c         print *, "COMPUTED", computed_fluxes(num-10:num)

         if (mod(n,10) .eq. 0 .and. debug .gt. 0) then
c            if (iraf .eq. 1) then
c               do j = 1,10
c                  dd(num-10+j) = 1. - d(num-10+j)
c               enddo
c               write (nf4out,1110) (dd(num-10+j),j=1,10)
c            endif
            if (iunits .eq. 1) then
               wave3 = 1.d-4*(wave - 9.0*step)
               write (nf1out,1112) wave3,(d(num-10+j),j=1,10)
            else
               wave3 = wave - 9.0*step
            write (nf1out,1111) wave3,(d(num-10+j),j=1,10)
            endif
            if (nf2out .gt. 0) write (nf2out,1110) (d(num-10+j),j=1,10)
         endif
         if (num .ge. 5000) num = 0
      enddo


c*****finish the synthesis
c      print *, "MODNUM 10", mod(num,10)
      nn = mod(num,10)
      if (nn .ne. 0 .and. debug .gt. 0) then
c         computed_fluxes(num-nn:num) = 1.0 - dd(num-nn:num)
c         print *, "NUM, NN", num, nn, dd(num-nn:num)
         if (iraf .eq. 1) then
            do j=1,nn
               dd(num-nn+j) = 1. - d(num-nn+j)
            enddo
            write (nf4out,1110) (dd(num-nn+j),j=1,nn)
         endif
         if (iunits .eq. 1) then
            wave3 = 1.d-4*(wave - 9.0*step)
            write (nf1out,1112) wave3,(d(num-nn+j),j=1,nn)
         else
            wave3 = wave - 9.0*step
            write (nf1out,1111) wave3,(d(num-nn+j),j=1,nn)
         endif
         if (nf2out .gt. 0) write (nf2out,1110) (d(num-nn+j),j=1,nn)
      endif
      if (debug .gt. 0) then
         if (iunits .eq. 1) then
            write (nf1out,1113) 1.d-4*wave
         else
            write (nf1out,1114) wave
         endif
      endif

c*****exit normally
       return 


c*****format statements
1001  format ('  kaplam from 1 to ntau at wavelength',f10.2/
     .        (6(1pd12.4)))
1002  format ('MODEL: ',a73)
1003  format ('AT WAVELENGTH/FREQUENCY =',f11.7,
     .        '  CONTINUUM FLUX/INTENSITY =',1p,d12.5)
1004  format ('AT WAVELENGTH/FREQUENCY =',f11.3,
     .        '  CONTINUUM FLUX/INTENSITY =',1p,d12.5)
1101  format (/'SPECTRUM DEPTHS')
1102  format (4f11.3)
1103  format (4f10.7)
1104  format ('SIMPLE  =    t'/'NAXIS   =     1'/'NAXIS1  = ',i10,/
     .        'W0      =',f10.4/'CRVAL1  =',f10.4/'WPC     =',f10.4/
     .        'CDELT1  =',f10.4)
1105  format (16HORIGIN  = 'moog'/21HDATA-TYP= 'synthetic'/
     .        18HCTYPE1  = 'lambda'/21HCUNIT1  = 'angstroms')
1106  format (11HTITLE   = ',A65,1H')
1107  format ('ATOM    = ',1H',7x,a2,1H',/,'ABUND   = ',f10.2)
1108  format ('VTURB   = ',d10.4,'     /  cm/sec  ')
1109  format ('END')
1110  format (10f7.4)
1111  format (f10.3,': depths=',10f6.3)
1112  format (f10.7,': depths=',10f6.3)
1113  format ('FINAL WAVELENGTH/FREQUENCY =',f10.7/)
1114  format ('FINAL WAVELENGTH/FREQUENCY =',f10.3/)


      end                                




