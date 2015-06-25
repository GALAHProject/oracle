
      subroutine inlines (num)
c******************************************************************************
c     This subroutine reads in the line data
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Dummy.com'
      include 'Pstuff.com'
      include 'Quants.com'
      include 'Factor.com'
      real*8 swave1(40), satom1(40), se(40),sgf(40),
     .       sdampnum(40),sd0(40),swidth(40), scharge(40)
      integer n2

c      if (num .eq. 2) go to 4
c      if (num .eq. 6) go to 340
      n1marker = 1
      n2 = 0


c*****decide if certain element abundances need to be modified.
c      print *, "inlines", numpecatom
      if (numpecatom .gt. 0) then
         do iatom=3,95
            xabund(iatom) = 10.**pecabund(iatom,isynth)*
     .                      10.**abfactor(isynth)*xabu(iatom)
c            print *, "xabund(iatom)", i, xabund(iatom)
         enddo
      endif
      if (num .ne. 5) then
         if (debug .gt. 0) write (*,1004)
         xmetals = abscale + abfactor(isynth)
         if (ninetynineflag .eq. 1) then
            if (debug .gt. 0) write (*,1005) xmetals
         else
            if (debug .gt. 0) write (*,1006) abscale
         endif
         do j=1,93
            if (pec(j) .gt. 0 ) then
               dummy1(j) = dlog10(xabund(j)) + 12.0
               if (dummy1(j) .le. -10.) then 
                  if (debug .gt. 0) write (*,1008) names(j),dummy1(j)
               else
                  if (debug .gt. 0) write (*,1007) names(j),dummy1(j)
               endif
            endif
         enddo
      endif


c*****output information about the isotopic ratios
      if (numiso .gt. 0) then
         do i=1,numiso
            iiso = isotope(i)
            if (debug .gt. 0) 
     .       write (*,1015) iiso, isotope(i), isoabund(i,isorun)
         enddo
      endif



c*****Inititalize strong line printing
c     if 'printstrong' gt 0 then the strong lines have 
c     been printed
      printstrong = -1

      if (num .ne. 4) then  
c         rewind nflines
         wave = start
c         read (nflines,1001) linitle
      endif


c*****read in the strong lines if needed
c 302   nstrong = 0
c      if (dostrong .gt. 0 ) then
c         rewind nfslines
c         do j=1,41
c            if (linfileopt .eq. 0) then
c               read (nfslines,1002,end=340) swave1(j),satom1(j),se(j),
c     .                             sgf(j),sdampnum(j),sd0(j),swidth(j)
c            else
c               read (nfslines,*,end=340) swave1(j),satom1(j),se(j),
c     .                             sgf(j),sdampnum(j),sd0(j),swidth(j)
c            endif
c            nstrong = nstrong + 1
c            iatom = satom1(j)
c            scharge(j) = 1.0 + dble(int(10.0*(satom1(j) - iatom)
c     .          +0.0001))
c            if (scharge(j) .gt. 3.) then
c               write (*,1003) swave1(i), satom1(i)
c               stop
c            endif
c         enddo
c         if (nstrong .gt. 40) then
c            write(*,*) 'STRONG LINE LIST HAS MORE THAN 40 LINES. THIS'
c            write(*,*) 'IS NOT ALLOWED. I QUIT!'
c            stop
c         endif
c      endif

c      nlines = 0 + nlines_absolute
c      nstrong = 0 + nstrong_absolute

c      print *, "inlines says number of lines is ", nlines, nstrong

      if (nstrong .gt. 40) then
        print *, "more than 40 strong lines. stahp"
        call f2pystop
c        stop
      endif

      if (nstrong .gt. 0) then
        do j=1,nstrong
          swave1(j) = strong_transitions(j, 1)
          satom1(j) = strong_transitions(j, 2)
          se(j) = strong_transitions(j, 3)
          sgf(j) = strong_transitions(j, 4)
          sdampnum(j) = strong_transitions(j, 5)
          sd0(j) = strong_transitions(j, 6)
          swidth(j) = strong_transitions(j, 7)

c         Check the charge
          iatom = satom1(j)
          scharge(j) = 1.0 + dble(int(10.0*(satom1(j) - iatom)
     .      + 0.0001))
          if (scharge(j) .gt. 3.) then
             print *, "strong line charge greater than 3"
             print *, "A", wave1(j), iatom, charge(j), gf(j)
             print *, "stahp"
             call f2pystop
c             stop
          endif
        enddo
      endif

      if (nlines .gt. 0) then
        do j=1,nlines
           wave1(j) = transitions(j, 1)
           atom1(j) = transitions(j, 2)
           e(j, 1) = transitions(j, 3)
           gf(j) = transitions(j, 4)
           dampnum(j) = transitions(j, 5)
           d0(j) = transitions(j, 6)
           width(j) = transitions(j, 7)

           iatom = atom1(j)
           charge(j) = 1.0 + dble(int(10.0*(atom1(j) - iatom)
     .       +0.0001))
           if (charge(j) .gt. 3.) then
              print *, "line charge greater than 3"
              print *, "B", wave1(j), iatom, charge(j), gf(j)
              print *, "stahp"
              call f2pystop
c              stop
           endif
        enddo
      endif

c 340   nlines = 2500 - nstrong
c      j = 1
c 333   if (linfileopt .eq. 0) then
c         read (nflines,1002,end=311) wave1(j),atom1(j),e(j,1),gf(j),
c     .                             dampnum(j),d0(j),width(j)
c      else
c         read (nflines,*,end=311) wave1(j),atom1(j),e(j,1),gf(j),
c     .                             dampnum(j),d0(j),width(j)
c      endif
c      iatom = atom1(j)
c      charge(j) = 1.0 + dble(int(10.0*(atom1(j) - iatom)+0.0001))
c      if (charge(j) .gt. 3.) then
c         write (*,1003) wave1(j), atom1(j)
c         stop
c      endif 
c      if (width(j) .lt. 0.) go to 333
c      if (iunits .eq. 1) wave1(j) = 1.d+4*wave1(j)
c      j = j + 1
c      if (j .le. nlines) go to 333
c 311   nlines = j - 1 


c*****append the strong lines here if necessary
      if (nstrong .gt. 0) then
         do k=1,nstrong
            wave1(nlines+k) = swave1(k)
            atom1(nlines+k) = satom1(k)
            e(nlines+k,1) = se(k)
            gf(nlines+k) = sgf(k)
            dampnum(nlines+k) = sdampnum(k)
            d0(nlines+k) = sd0(k)
            width(nlines+k) = swidth(k)
            charge(nlines+k) = scharge(k)
         enddo
      endif


c*****here groups of lines for blended features are defined
      do j=1,nlines+nstrong
         if (wave1(j) .lt. 0.) then
            group(j) = 1
            wave1(j) = dabs(wave1(j))
            width(j) = width(j-1)
         else
            group(j) = 0
         endif
      enddo     


c*****here excitation potentials are changed from cm^-1 to eV, if needed
      do j=1,nlines+nstrong
         if (e(j,1) .gt. 50.) then
            do jj=1,nlines+nstrong
               e(jj,1) = 1.2389e-4*e(jj,1)
            enddo
            exit
         endif
      enddo
 

c*****here log(gf) values are turned into gf values, if needed
      do j=1,nlines+nstrong
         if (gfstyle.eq.1 .or. gf(j) .lt. 0) then
            do jj=1,nlines+nstrong
               gf(jj) = 10.**gf(jj)
            enddo
            exit
         endif
      enddo         


c*****turn log(RW) values and EW values in mA into EW values in A.  Stuff
c     duplicate EW values of the first line of a blend into all blend members.
      do j=1,nlines+nstrong
         if (width(j) .lt. 0.) then
            width(j) = 10.**width(j)*wave1(j)
         else
            width(j) = width(j)/1000.
         endif
       enddo


c*****here some parameters for the lines are assigned or calculated; 
c     there is a block of statements for moleculer lines, 
c     and a different one for atomic lines
      do j=1,nlines+nstrong
         iatom = atom1(j)
         atom10 = 10.*atom1(j)
         e(j,2) =  e(j,1) + 1.239d+4/wave1(j)


c*****here are the calculations specific to molecular lines
         if (iatom .ge. 100) then
            call sunder (atom1(j),ia,ib)
            if (ia .gt. ib) then
c               write (*,1010) ia,ib
               print *, "stahp ia, ib", ia, ib
               call f2pystop
c               stop
            endif
            if (atom10-int(atom10) .le. 0.0) then
               amass(j) = xam(ia) + xam(ib)    
               mas1 = xam(ia) + 0.0000001
               mas2 = xam(ib) + 0.0000001
            else 
               jat100 = int(100.*(atom10+0.00001))
               mas1 = jat100 - 100*int(atom10)
               jat10000 = int(10000.*(atom10+0.00001))
               mas2 = jat10000 - 100*jat100
               if (mas1.gt.mas2 .or. mas1.le.0.0 .or. 
     .             mas2.le.0.0) then
c                  write (*,1011) mas1, mas2
                  print *, "stahp mas1, mas2", mas1, mas2
                  call f2pystop
c                  stop
               endif
               amass(j) = mas1 + mas2
            endif
c*****use an internal dissociation energy for molecules if the user
c     does not read one in
            if (d0(j) .eq. 0.) then
               do k=1,110
                  if (int(datmol(1,k)+0.01) .eq.
     .                int(atom1(j)+0.01)) then
                     d0(j) = datmol(2,k)
                     go to 390
                  endif
                enddo
                if (debug .gt. 0) write (*,1013) atom1(j)
                call f2pystop
c                stop
            endif
390         rdmass(j) = mas1*mas2/amass(j)
            chi(j,1) = 0.
            chi(j,2) = 0.
            chi(j,3) = 0.


c*****here are the calculations specific to atomic lines
         else
c            print *, "atom10", atom10, int(atom10), xam(iatom)
            if (atom10-int(atom10) .le. 0.0) then
               amass(j) = xam(iatom)
            else 
               atom10 = atom10 + 0.00001
               amass(j) = int(1000*(atom10-int(atom10)))
            endif
            rdmass(j) = 0.
            chi(j,1) = xchi1(iatom)
            chi(j,2) = xchi2(iatom)
            chi(j,3) = xchi3(iatom)
         endif
      enddo


c*****quit the routine normally
      if (nlines+nstrong .lt. 2500) then
         if (sstop .gt. wave1(nlines)+10.) sstop = wave1(nlines)+10.
      endif
      lim1line = 1
      return  


c****prepare to get another chunk of line data 
c 4     n2 = n1marker + lim1line - 1
c      n1marker = n2
c      rewind nflines
c      do j=1,n2
c         read (nflines,1001)
c      enddo
c      start = wave
c      go to 302


c*****format statements
1001  format (a80)
1002  format (7e10.3)
1003  format ('INPUT STRONG LINE: LAMBDA = ', f10.3, ' AND ID = ',
     .        f6.1, ' CANNOT BE DONE!'/
     .        'NO TRIPLE OR GREATER IONS; I QUIT!')
1004  format (/'For these computations, some abundances have ',
     .       ' been altered:')
1005  format ('Changing overall metallicity: ', f6.2, ' dex')
1006  format ('ALL abundances NOT listed below differ ',
     .        'from solar by ', f6.2, ' dex')
1007  format ('element ', a2, ':  abundance = ', f5.2)
1008  format ('element ', a2, ':  abundance = ', f5.1)
1010  format ('ATOMIC NUMBERS IN MOLECULAR NAME (',
     .        2i2, ') ARE IN WRONG ORDER'/'I QUIT!!!')
1011  format ('ISOTOPIC MASS NUMBERS IN MOLECULAR NAME (',
     .        2i3, ') ARE IN WRONG ORDER OR ARE WEIRD;'/'I QUIT!!!')
1013  format (f6.1, ' IS AN UNKOWN MOLECULE; I QUIT!')
1014  format ('Isotopic Ratios given for this synthesis')
1015  format ('Isotopic Ratio: [', i4, '/', f10.5, '] = ', f10.3)
      

      end



