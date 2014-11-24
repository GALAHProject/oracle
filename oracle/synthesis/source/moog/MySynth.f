
      function synthesise(teff_, logg_, mh_, vturb_,
     .   photospheric_structure_, photospheric_abundances,
     .   transitions_, synlimits_, opacity_contribute_, 
     .   npoints_, debug_, ntau_, ncols_, natoms_,
     .   nlines_, wavelengths, fluxes)

      implicit real*8 (a-h,o-z)
      real*8, intent(in) :: teff_, logg_, mh_, vturb_
      real*8, intent(in) :: opacity_contribute_ 
      real*8, dimension(ntau_, ncols_), intent(in) ::
     .   photospheric_structure_
      real*8, dimension(3) :: synlimits_
      real*8, dimension(natoms_, 2), intent(in) ::
     .   photospheric_abundances
      real*8, dimension(nlines_, 7), intent(in) :: transitions_
      integer :: npoints_
      integer, optional :: debug_
      real*8, dimension(npoints_), intent(out) :: wavelengths
      real*8, dimension(npoints_), intent(out) :: fluxes

c******************************************************************************
c     This program synthesizes a section of spectrum and compares it
c     to an observation file.
c******************************************************************************

      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'


cc      i = nint((synlimits_(2) - synlimits_(1))/synlimits_(3) + 1)
c      print *, "SETTING AS ",i
c      allocate (output(2001))

      if (opacity_contribute_ .lt. 1.0) then
         delta = 1.0
         olddelta = 1.0
      else
         delta = 0.0 + opacity_contribute_
         olddelta = 0.0 + opacity_contribute_
      endif

      start = synlimits_(1)
      oldstart = synlimits_(1)

      sstop = synlimits_(2)
      oldstop = synlimits_(2)

      step = synlimits_(3)
      oldstep = synlimits_(3)
      step1000 = 1000. * step

c      print *, "start stop", start, sstop, step
c      print *, "delta ", delta, olddelta

      debug = debug_
      control = 'synth   '
      silent = 'y'
      smterm = 'x11'
      smt1 = 'x11'
      smt2 = 'x11'
   
      nstrong = 0
      transitions(:nlines_, :) = transitions_
c      if (nstrong .gt. 0) then
c         strong_transitions(:nstrong_, :) = strong_transitions_
c      endif

c     MOOG plays with these. So let's keep absolute reference values
      nlines_absolute = 0 + nlines_
      nstrong_absolute = 0 + nstrong

      vturb_absolute = vturb_
c      print *, "input vturb", vturb_, vturb_absolute
      

      logepsilon_abundances(:natoms_, :) = photospheric_abundances
      photospheric_structure(:ntau_, :ncols_) = photospheric_structure_

c     These should not change...
      ntau = ntau_
      modtype = 'WEBMARCS'
      moditle = 'atmosphere comment'
      natoms = natoms_
      abscale = mh_

c*****examine the parameter file
      call params

      
      numpecatom = 0 + natoms_
      numatomsyn = 1

      do i=1,numpecatom
         jatom = photospheric_abundances(i, 1)
         if (jatom .eq. 99) then
            do kk=1,numatomsyn
c               abfactor(kk) = photospheric_abundances(i, 1+kk)
                abfactor(kk) = 0.0
            enddo
         else
            do kk=1,numatomsyn
c               pecabund(jatom, kk) = photospheric_abundances(i, 1+kk)
                pecabund(jatom, kk) = 0.0
            enddo
            pec(jatom) = 1
         endif
      enddo

c do l=1,numpecatom
c            read (nfparam,*) jatom,(deltalogab(kk),kk=1,numatomsyn)
c            linecount = linecount + 1
c            if (jatom .eq. 99) then
c               do kk=1,numatomsyn
c                  abfactor (kk) = deltalogab(kk)
c               enddo
c            else
c               do kk=1,numatomsyn
c                  pecabund(jatom,kk) = deltalogab(kk)
c               enddo
c               pec(jatom) = 1
c            endif
c         enddo


c*****open the files for: standard output, raw spectrum depths, smoothed 
c     spectra, and (if desired) IRAF-style smoothed spectra
c      nf1out = 20     
c      lscreen = 4
c      array = 'STANDARD OUTPUT'
c      nchars = 15
c      call infile ('output ',nf1out,'formatted  ',0,nchars,
c     .             f1out,lscreen)
c      nf2out = 21               
c      lscreen = lscreen + 2
c      array = 'RAW SYNTHESIS OUTPUT'
c      nchars = 20
c      call infile ('output ',nf2out,'formatted  ',0,nchars,
c     .             f2out,lscreen)
c      if (plotopt .ne. 0) then
c         nf3out = 22               
c         lscreen = lscreen + 2
c         array = 'SMOOTHED SYNTHESES OUTPUT'
c         nchars = 25
c         call infile ('output ',nf3out,'formatted  ',0,nchars,
c     .                f3out,lscreen)
c         if (f5out .ne. 'optional_output_file') then
c            nf5out = 26
c            lscreen = lscreen + 2
c            array = 'POSTSCRIPT PLOT OUTPUT'
c            nchars = 22
c            call infile ('output ',nf5out,'formatted  ',0,nchars,
c     .                   f5out,lscreen)
c         endif
c      endif
c      if (iraf .ne. 0) then
c         nf4out = 23               
c         lscreen = lscreen + 2
c         array = 'IRAF ("rtext") OUTPUT'
c         nchars = 24
c         call infile ('output ',nf4out,'formatted  ',0,nchars,
c     .                f4out,lscreen)
c      endif


c*****open and read the model atmosphere file
c      nfmodel = 30
c      lscreen = lscreen + 2
c      array = 'THE MODEL ATMOSPHERE'
c      nchars = 20
c      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
c     .             fmodel,lscreen)
      call inmodel


c*****open the line list file and the strong line list file
c      nflines = 31
c      lscreen = lscreen + 2
c      array = 'THE LINE LIST'
c      nchars = 13
c      call infile ('input  ',nflines,'formatted  ',0,nchars,
c     .              flines,lscreen)
c      if (dostrong .gt. 0) then
c         nfslines = 32
c         lscreen = lscreen + 2
c         array = 'THE STRONG LINE LIST'
c         nchars = 20
c         call infile ('input  ',nfslines,'formatted  ',0,nchars,
c     .                 fslines,lscreen)
c      endif
c      print *, "numpecatom", numpecatom, numatomsyn
      if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
         isynth = 1
         isorun = 1
         nlines = 0
         mode = 3
         call inlines (1)
         call eqlib
         call nearly (1)
         call synspec
      else
         do n=1,numatomsyn
            isynth = n
            isorun = n
            start = oldstart
            sstop = oldstop
            mode = 3
            molopt = 2
            call inlines (1)
            call eqlib
            call nearly (1)
            call synspec
            linprintopt = 0
         enddo
      endif
      
         

c*****do the syntheses
c      choice = '1'
c      do i=1,100
c         if (i .eq. 100) then
c            write (*,1002) 
c            stop
c         endif
c         ncall = 1
c         call getsyns (lscreen,ncall)
         

c*****now either don't make a plot (plotopt = 0) 
c                plot the synthetic spectrum, (plotopt = 1)
c                plot syntheses and observation (plotopt = 2) 
c                or just smooth the syntheses (plotopt = 3)
c         if (choice .eq. 'n') then
c            ncall = 2
c         else
c            ncall = 1
c         endif
c         if     (plotopt .eq. 0) then
c            choice = 'q'
c         elseif (plotopt .eq. 1) then
c            call pltspec (lscreen,ncall)
c         elseif (plotopt .eq. 2) then
c            nfobs = 33               
c            lscreen = lscreen + 2
c            array = 'THE OBSERVED SPECTRUM'
c            nchars = 21
c            if (specfileopt .eq. 1) then
c               call infile ('input  ',nfobs,'unformatted',2880,nchars,
c     .                      fobs,lscreen)
c            else
c               call infile ('input  ',nfobs,'formatted  ',0,nchars,
c     .                      fobs,lscreen)
c            endif
c            call pltspec (lscreen,ncall)
c         elseif (plotopt .eq. 3) then
c            call smooth (-1,ncall)
c            choice = 'q'
c         else
c            write (*,1001)
c            stop
c         endif
c         if (choice .eq. 'q') then
c            call finish (0)
c            exit
c         endif
c      enddo

c      print *, "nkount", nkount
c      print *, "computed_wls", computed_wls(:10)
c      print *, "computed_fluxes", computed_fluxes(:10) 

c      print *, "computed_wls", computed_wls(nkount-10:nkount)
c      print *, "computed_fluxes", computed_fluxes(nkount-10:nkount) 

c      allocate (output(nkount))


      wavelengths = computed_wls(:nkount)
      fluxes = computed_fluxes(:nkount)
      return

c*****format statements
1001  format ('for syntheses, parameter "plot" must be 0, 1, 2, or 3;',
     .        ' I QUIT!')
1002  format ('something wrong:  max number (99) of synthesis ',
     .        'cycles exceeded; I QUIT!')



      end 





