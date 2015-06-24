
      function abundances(in_metallicity, in_xi,
     .   in_photosphere, in_logepsilon_abundances,
     .   in_transitions, in_modtype, in_debug, output, in_ntau, 
     .   in_ncols, in_natoms, in_nlines)

      implicit real*8 (a-h,o-z)
      real*8, intent(in) :: in_metallicity, in_xi
      real*8, dimension(in_ntau, in_ncols), intent(in) ::
     .   in_photosphere
      real*8, dimension(in_natoms, 2), intent(in) ::
     .   in_logepsilon_abundances
      real*8, dimension(in_nlines, 7), intent(in) :: in_transitions
      character*10, intent(in) :: in_modtype
      integer, optional :: in_debug

      real*8, dimension(in_nlines), intent(out) :: output
  
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      include 'Mol.com'
      include 'Pstuff.com'
      include 'Dampdat.com'

cc      include 'Factor.com' 


c     Params crap. Can probably be removed eventually.
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
      dampingopt   = 1
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


      dummy1(:) = 0.
      dummy2(:) = 0.
      dummy3(:) = 0.
      dummy4(:) = 0.



      t(:) = 0.0
      theta(:) = 0.0 
      tkev(:) = 0.0 
      tlog(:) = 0.0 
      pgas(:) = 0.0 
      ne(:) = 0.0 
      nhtot(:) = 0.0 
      numdens(:,:,:) = 0.0 
      molweight(:) = 0.0 

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
cc      pmol(:) = 0.0
cc      xmol(:,:) = 0.0
cc      xamol(:,:) = 0.0
cc      xatom(:) = 0.0
cc      patom(:) = 0.0
cc      amol(:) = 0.0
cc      smallmollist(:) = 0.0
cc      largemollist(:) = 0.0
cc      datmol(:,:) = 0.0
cc      const(:,:) = 0.0









c  INITIALIZE SOME VARIABLES: spectrum run parameters
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

c*****open data files carried with the source code: Barklem UV damping
c      nfbarklem = 35
c      num = 60
c      call getcount (num,moogpath)
c      if (moogpath(num:num) .ne. '/') then
c         num = num + 1
c         moogpath(num:num) = '/'
c      endif
c      fbarklem(1:num) = moogpath(1:num)
c      fbarklem(num+1:num+11) = 'Barklem.dat'
      nfbarklem = 35
      open (nfbarklem,file="$DATA_DIR/Barklem.dat")

 
c      nfbarklemUV = 36
c      num = 60
c      call getcount (num,moogpath)
c      if (moogpath(num:num) .ne. '/') then
c         num = num + 1
c         moogpath(num:num) = '/'
c      endif
c      fbarklemUV(1:num) = moogpath(1:num)
c      fbarklemUV(num+1:num+13) = 'BarklemUV.dat'
      nfbarklemUV = 36
      open (nfbarklemUV,file="$DATA_DIR/BarklemUV.dat")




c     Pass information to the global variables
      debug = in_debug

      modtype = in_modtype
      ntau = in_ntau
      vturb(1) = in_xi

      moditle = 'atmosphere comment'
      
      natoms = in_natoms
      abscale = in_metallicity
      
      nstrong = 0
      nlines = 0 + int(in_nlines)
      

      control = 'abfind  '
      silent = 'y'
      smterm = 'x11'
      smt1 = 'x11'
      smt2 = 'x11'

      do i=1,natoms
         element(i) = in_logepsilon_abundances(i, 1)
         logepsilon(i) = in_logepsilon_abundances(i, 2)
      enddo

      photospheric_structure(:in_ntau, :in_ncols) = in_photosphere


      call inmodel


      transitions(:in_nlines, :) = in_transitions

c     MOOG plays with these. So let's keep absolute reference values
      
      


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


      mode = 0
100   call fakeline

      nstrong = 0
      nlines = 0 + int(in_nlines)
      call inlines (1) 
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
         if (debug .gt. 0) print *, "finishing"
         return
      endif
      lim1obs = lim1line
      lim2obs = lim2line

c*****find out whether molecular equilibrium is involved in the species
      call molquery


c*****force each abundance of a species member to predict the 
c     line equivalent width; here is the code for ordinary species
      if (molflag .eq. 0) then
         abundin =  dlog10(xabund(iabatom)) + 12.0
         do lim1=lim1line,lim2line
            call lineabund(abundin)
         enddo
c         call stats
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
                  xabund(iabatom) = 10.**(average-12.)
                  iternumber = iternumber + 1
                  call eqlib
                  call nearly (2)
                  go to 10
               else
                  if (debug .gt. 0) write (array,1001) iternumber
                   ikount=kount+14
                  nchars = 53
c                  call putasci (nchars,ikount)
               endif
            else
               if (debug .gt. 0) 
     .             write (array,1003) molecule, dlog10(abundin),
     .                            dlog10(average)
c               lscreen = lscreen + 2
c               call prinfo (lscreen)
               print *, "molecule!"
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
   
      if (debug .gt. 0) print *, "abundances", abundout(1:10)
      output = abundout(1:nlines)
       
      choice = 'y'
      nchars = 0
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
      if (choice.eq.'y' .or. nchars.le.0) then
         if (mode .eq. 2) then
            go to 5
         else
c            call finish (0)
            return
         endif
      else
c         call finish (0)
         return
      endif

      return

c*****format statements
1001  format ('THIS REQUIRED', i2,' ITERATIONS WITH MOLECULAR ',
     .        'EQUILIBRIUM')
1002  format (a80)
1003  format ('FOR SPECIES ', f10.1,' NO CONVERGENCE: '/
     .        'LAST ITERATIONS YIELD', 2f10.3, '  I QUIT!')


       end
