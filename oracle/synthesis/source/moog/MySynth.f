
      function synthesise(in_metallicity, in_xi,
     .   in_photosphere, in_logepsilon_abundances,
     .   in_transitions, in_synlimits, in_opacity_contributes, 
     .   in_modtype, damping, in_npoints, in_debug, wavelengths, fluxes,
     .   data_path, in_ntau, in_ncols, in_natoms, in_nlines)

      implicit real*8 (a-h,o-z)
      real*8, intent(in) :: in_metallicity, in_xi
      real*8, dimension(in_ntau, in_ncols), intent(in) ::
     .   in_photosphere
      real*8, dimension(in_natoms, 2), intent(in) ::
     .   in_logepsilon_abundances
      real*8, dimension(in_nlines, 7), intent(in) :: in_transitions
      real*8, dimension(3) :: in_synlimits
      real*8, intent(in) :: in_opacity_contributes 
      integer :: damping, in_npoints
      character*10, intent(in) :: in_modtype
      integer, optional :: in_debug
      character(300), intent(in) :: data_path

      real*8, dimension(in_npoints), intent(out) :: wavelengths
      real*8, dimension(in_npoints), intent(out) :: fluxes


      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      include 'Dampdat.com'


      nfbarklem = 35
      open (nfbarklem,file=TRIM(data_path) // "/Barklem.dat")

      nfbarklemUV = 36
      open (nfbarklemUV,file=TRIM(data_path) // "/BarklemUV.dat")


      nfmodel =  0 
      nflines =  0
      nfslines = 0
      nfobs =    0
      nftable =  0
      modprintopt  = 2
      molopt       = 2
      linprintopt  = 2
      fluxintopt   = 0
      plotopt      = 0
      dampingopt   = 0 + damping
      specfileopt  = 0
      linfileopt   = 0
      iunits       = 0
      itru         = 0
      iscale       = 0
      iraf         = 0
      histoyes     = 0
      byteswap     = 0
      deviations   = 0
      scatopt      = 0
      gfstyle      = 0
      maxshift     = 0
      dostrong     = 0
      fudge = -1.0

      oldstart = 0.
      start = 0.
      sstop = 0.
      step = 0.
      delta = 0.
      cogatom = 0.
      contnorm = 1.0


c  INITIALIZE SOME VARIABLES: line limit parameters
      ncurve = 0
      lim1line = 0
      lim2line = 0
      lim1obs = 0
      lim2obs = 0
      lim1 = 0
      lim2 = 0


cc      i = nint((in_synlimits(2) - in_synlimits(1))/in_synlimits(3) + 1)
c      print *, "SETTING AS ",i
c      allocate (output(2001))

      modtype = in_modtype

      if (in_opacity_contributes .lt. 1.0) then
         delta = 1.0
         olddelta = 1.0
      else
         delta = 0.0 + in_opacity_contributes
         olddelta = 0.0 + in_opacity_contributes
      endif

      start = in_synlimits(1)
      oldstart = in_synlimits(1)

      sstop = in_synlimits(2)
      oldstop = in_synlimits(2)

      step = in_synlimits(3)
      oldstep = in_synlimits(3)
      step1000 = 1000. * step

c      print *, "start stop", start, sstop, step
c      print *, "delta ", delta, olddelta

      debug = in_debug
      control = 'synth   '
      silent = 'y'
      smterm = 'x11'
      smt1 = 'x11'
      smt2 = 'x11'
   
c       if (nstrong .gt. 0) then
c         strong_transitions(:nstrong_, :) = strong_in_transitions
c      endif

c     MOOG plays with these. So let's keep absolute reference values
      nlines = 0 + in_nlines
      nstrong = 0 
      
      transitions(:nlines, :) = in_transitions

      vturb(1) = in_xi


      logepsilon_abundances(:in_natoms, :) = in_logepsilon_abundances
      photospheric_structure(:in_ntau, :in_ncols) = in_photosphere

c     These should not change...
      ntau = in_ntau
      moditle = 'atmosphere comment'
      natoms = in_natoms
      abscale = in_metallicity

c*****examine the parameter file
c      call params

      
      numpecatom = 0 + in_natoms
      numatomsyn = 1

      do i=1,numpecatom
         jatom = in_logepsilon_abundances(i, 1)
         if (jatom .eq. 99) then
            do kk=1,numatomsyn
               abfactor(kk) = in_logepsilon_abundances(i, 1+kk)
            enddo
         else
            do kk=1,numatomsyn
               pecabund(jatom, kk) = in_logepsilon_abundances(i, 1+kk)
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
      if (debug .gt. 0) print *, "numpecatom", numpecatom, numatomsyn

      if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
         isynth = 1
         isorun = 1
c         nlines = 0
         mode = 3
c        inlines updates nlines and nstrong, so we will have to update them

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

      wavelengths = computed_wls(:in_npoints)
      fluxes = computed_fluxes(:in_npoints)
      return

      call f2pystop

c*****format statements
1001  format ('for syntheses, parameter "plot" must be 0, 1, 2, or 3;',
     .        ' I QUIT!')
1002  format ('something wrong:  max number (99) of synthesis ',
     .        'cycles exceeded; I QUIT!')

   

      end 





