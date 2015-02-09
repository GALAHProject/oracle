
      function abundances(mh_, vturb_,
     .   photospheric_structure_, photospheric_abundances,
     .   transitions_, modtype_, debug_, output, ntau_, 
     .   ncols_, natoms_, nlines_)

      implicit real*8 (a-h,o-z)
      real*8, intent(in) :: mh_, vturb_
      real*8, dimension(ntau_, ncols_), intent(in) ::
     .   photospheric_structure_
      real*8, dimension(natoms_, 2), intent(in) ::
     .   photospheric_abundances
      real*8, dimension(nlines_, 7), intent(in) :: transitions_
      character*10, intent(in) :: modtype_
      integer, optional :: debug_

      real*8, dimension(nlines_), intent(out) :: output


c      if (.not. present(debug_)) then
c         debug_ = 0
c      endif

c      real*8, dimension(:), intent(inout) :: abundances
c      real*8, dimension(nstrong_, 7), intent(in):: strong_transitions_

c      real*8, dimension(nlines_, 1), intent(inout) :: line_abundances


c      real*8, dimension(1, 7) :: strong_transitions_

   
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'






      t(:) = 0.0
      theta(:) = 0.0 
      tkev(:) = 0.0 
      tlog(:) = 0.0 
      pgas(:) = 0.0 
      ne(:) = 0.0 
      nhtot(:) = 0.0 
      numdens(:,:,:) = 0.0 
      molweight(:) = 0.0 
      vturb(:) = 0.0 
      scont(:) = 0.0 
      kapref(:) = 0.0 
      kaplam(:) = 0.0 
      tauref(:) = 0.0 
      taulam(:) = 0.0 
      kaplamabs(:) = 0.0 
      kaplamsca(:) = 0.0 
      rho(:) = 0.0 
      rhox(:) = 0.0 
      xref(:) = 0.0 
      xdepth(:) = 0.0 
      elem(:) = 0.0 
      xabund(:) = 0.0 
      xabu(:) = 0.0 
      u(:,:,:) = 0.0 

      a(:,:) = 0.0
      dopp(:,:) = 0.0
      kapnu0(:,:) = 0.0
      gf(:) = 0.0
      wave1(:) = 0.0
      atom1(:) = 0.0
      e(:,:) = 0.0
      chi(:,:) = 0.0
      amass(:) = 0.0
      charge(:) = 0.0
      d0(:) = 0.0
      dampnum(:) = 0.0
      gf1(:) = 0.0
      width(:) = 0.0
      abundout(:) = 0.0
      widout(:) = 0.0
      strength(:) = 0.0
      rdmass(:) = 0.0
      gambark(:) = 0.0
      alpbark(:) = 0.0
      gamrad(:) = 0.0
      wid1comp(:) = 0.0
      computed_wls(:) = 0.0
      computed_fluxes(:) = 0.0
      kapnu(:) = 0.0
      taunu(:) = 0.0
      cd(:) = 0.0
      sline(:) = 0.0
      d(:) = 0.0
      dellam(:) = 0.0
      w(:) = 0.0
      rwtab(:) = 0.0
      gftab(:) = 0.0


c     fails after this.
      pmol(:) = 0.0
      xmol(:,:) = 0.0
      xamol(:,:) = 0.0
      xatom(:) = 0.0
      patom(:) = 0.0
      amol(:) = 0.0
      print *, "smallmollist", smallmollist(1), largemollist(1)
      print *, "const", datmol(1,1), const(1,1)
c      smallmollist(:) = 0.0
c      largemollist(:) = 0.0
c      datmol(:,:) = 0.0
c      const(:,:) = 0.0


c     Params crap. Can probably be removed eventually.
      nfmodel =  0 
      nflines =  0
      nfslines = 0
      nfobs =    0
      nftable =  0
      modprintopt  = 1
      molopt       = 1
      linprintopt  = 1
      fluxintopt   = 0
      plotopt      = 0
      dampingopt   = 1
      specfileopt  = 0
      linfileopt   = 0
      iunits       = 0
      itru         = 0
      iscale       = 0
      nlines       = 0
      iraf         = 0
      histoyes     = 0
      byteswap     = 0
      deviations   = 0
      scatopt      = 0
      gfstyle      = 0
      maxshift     = 0
      dostrong     = 0
      molset       = 1
      fudge = -1.0


c  INITIALIZE SOME VARIABLES: spectrum run parameters
c      oldstart = 0.
c      start = 0.
c      sstop = 0.
c      step = 0.
c      delta = 0.
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

c      May need to initialise these later:
c      numpecatom = 0
c      numatomsyn = 0
c      newnumpecatom = 0
c      newnumatomsyn = 0
c      ninetynineflag = 0
c      pec(:) = 0
c      newpec(:) = 0
c      abfactor(:) = 0
c      pecabund(:, :) = 0
c      newpecabund(:, :) = 0.
c      numiso = 0
c      numisosyn = 0
c      newnumiso = 0
c      newnumisosyn = 0
c      isotope(:) = 0.0
c      newisotope(:) = 0.0
c      isoabund(:,:) = 0.0
c      newisoabund(:,:) = 0.0













      modtype = modtype_


c      print *, "teff_ etc", teff_, logg_, mh_, vturb_
c      print *, "photospheric_struct", photospheric_structure_
c      print *, "photospheric_abundances", photospheric_abundances_
c      print *, "transitions_", transitions_
c      print *, "ntau_", ntau_, ncols_, natoms_, nlines_
c
c      print *, "ntau_ ncols_ natoms_ nlines_", ntau_, ncols_,
c     .   natoms_, nlines_

c      i = nlines_ + nstrong_
c      allocate(temp_line_abundances(1:i))

      debug = debug_
      control = 'abfind  '
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
      moditle = 'atmosphere comment'
      natoms = natoms_
      abscale = mh_



c*****open the files for standard output and summary abundances
c      nf1out = 20
c      lscreen = 4
c      array = 'STANDARD OUTPUT'
c      nchars = 15
c      call infile ('output ',nf1out,'formatted  ',0,nchars,
c     .             f1out,lscreen)
c      nf2out = 21
c      lscreen = lscreen + 2
c      array = 'SUMMARY ABUNDANCE OUTPUT'
c      nchars = 24
c      call infile ('output ',nf2out,'formatted  ',0,nchars,
c     .             f2out,lscreen)
c      nf5out = 26
c      lscreen = lscreen + 2
c      array = 'POSTSCRIPT PLOT OUTPUT'
c      nchars = 22
c      call infile ('output ',nf5out,'formatted  ',0,nchars,
c     .             f5out,lscreen)


c*****open and read the model atmosphere
c      array = 'THE MODEL ATMOSPHERE'
c      nchars = 20
c 102   nfmodel = 30
c      lscreen = lscreen + 2
c      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
c     .             fmodel,lscreen)
      call inmodel


c*****open and read the line list file; get ready for the line calculations
c      array = 'THE LINE LIST'
c      nchars = 13
c      nflines = 31
c      lscreen = lscreen + 2
c      call infile ('input  ',nflines,'formatted  ',0,nchars,
c     .              flines,lscreen)


c*****get ready for the line calculations: generate a curve-of-growth lookup
c     table, read the linelist, etc.
c      print *, "doing fake line", nlines
c      print *, "vturb", vturb(:ntau)

100   call fakeline

c      print *, "calling inlines(1)"
      call inlines (1)
c      print *, "number of lines now", nlines, nstrong
c      print *, "and the absolute value says", nlines_absolute, nstrong_absolute
c      print *, "and the input value says", nlines_, nstrong_
      call eqlib
      call nearly (1)
 
c*****set some parameters
      ewsynthopt = -1
      mode = 2
      cogatom = 0.
      lim1line = 0


c*****find the range of lines of a species
5     call linlimit
      if (lim1line .lt. 0) then
c         call finish (0)
c         print *, "finishing"
         return
      endif
      lim1obs = lim1line
      lim2obs = lim2line
      lim2line = 0 + nlines_
c      print *, "lim1line, lim2line", lim1line, lim2line


c*****find out whether molecular equilibrium is involved in the species
      call molquery


c*****force each abundance of a species member to predict the 
c     line equivalent width; here is the code for ordinary species
      if (molflag .eq. 0) then
c         print *, "IBATOM", iabatom, dlog10(xabund(iabatom)) + 12.0
         abundin =  dlog10(xabund(iabatom)) + 12.0
c         print *, "ABUNDIN IS GOING TO BE", abundin
         do lim1=lim1line,lim2line
            call lineabund(abundin)
         enddo
c         print *, "CALLING STATS NOW--------------------"
c         print *, "array", array
c         call stats
c         print *, "CALLLLING LINE INFO NOW--------------"
         call lineinfo (3)
      else


c*****and here is the code for species involved in molecular equilibrium;
c     this procedure is iterated until input and output abundances are in 
c     agreement
         iternumber = 1
10       abundin =  dlog10(xabund(iabatom)) + 12.0
         do lim1=lim1line,lim2line
            call lineabund(abundin)
         enddo
c         call stats
         call lineinfo (3)
         if (t(jtau5).lt.3800                 .or. 
     .       atom1(lim1line).gt.100.0         .or.
     .       int(atom1(lim1line)+0.0001).eq.6 .or.
     .       int(atom1(lim1line)+0.0001).eq.8) then
            if (iternumber .lt. 6) then
               if (dabs(average-abundin) .gt. 0.02) then
c                  print *, "UPDATING FROM MOL EXQUILIB"
                  xabund(iabatom) = 10.**(average-12.)
                  iternumber = iternumber + 1
                  call eqlib
                  call nearly (2)
                  go to 10
               else
c                  write (array,1001) iternumber
                   ikount=kount+14
c                  nchars = 53
c                  call putasci (nchars,ikount)
               endif
            else
c               write (array,1003) molecule, dlog10(abundin),
c     .                            dlog10(average)
c               lscreen = lscreen + 2
c               call prinfo (lscreen)
c               print *, "molecule!"
               stop
            endif
         endif
      endif


c*****here a plot may be made on the terminal (and paper) if there 
c     are enough lines; then the user will be prompted on some
c     options concerning what is seen on the plot
c      if (plotopt .ne. 0) then
c         call pltabun
c         if     (choice.eq.'v') then
c            rewind nf1out
c            rewind nf2out
c            write (nf2out,1002) linitle,moditle
c            choice = ' '
c            go to 100
c         elseif (choice .eq. 'm') then
c            close (unit=nfmodel)
c            close (unit=nflines)
c            rewind nf1out
c            rewind nf2out
c            rewind nf5out
c            array = 'THE NEW MODEL ATMOSPHERE'
c            nchars = 24
c            fmodel =  'no_filename_given'
c            lim1line = 0
c            lim2line = 0
c            lim1obs = 0
c            lim2obs = 0
c            go to 102
c         endif
c      endif

 
             
c*****quit, or go on to another species?
c      if (silent .eq. 'y') then
c         choice = 'y'
c         nchars = 0
c      else
c         array = 'DO ANOTHER SPECIES ([y]/n)? '
c         nchars = 28
c         call getasci (nchars,maxline)
c         choice = chinfo(1:1)
c      endif
c      if (choice.eq.'y' .or. nchars.le.0) then
c         if (mode .eq. 2) then
c            go to 5
c         else
c            call finish (0)
c            return
c         endif
c      else
c         call finish (0)
c         return
c      endif
c       j = SUM(photospheric_structure(1,:))
c       call finish(0)
c       print *, "finishing"
c       print *, "abundances", abundout(:nlines)

       output = abundout(1:nlines)
c       abundances = abundout(:nlines)
c       what = abundout(:nlines)

       return

c*****format statements
1001  format ('THIS REQUIRED', i2,' ITERATIONS WITH MOLECULAR ',
     .        'EQUILIBRIUM')
1002  format (a80)
1003  format ('FOR SPECIES ', f10.1,' NO CONVERGENCE: '/
     .        'LAST ITERATIONS YIELD', 2f10.3, '  I QUIT!')


       end
