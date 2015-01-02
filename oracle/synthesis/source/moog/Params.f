
      subroutine params 
c****************************************************************************** 
c     This subroutine reads in the commands from the parameter file 
c****************************************************************************** 

      implicit real*8 (a-h,o-z) 
      include 'Atmos.com' 
      include 'Linex.com' 
      include 'Factor.com' 
      include 'Pstuff.com' 
      include 'Obspars.com'
      include 'Mol.com'
      include 'Multistar.com'
      real*8 deltalogab(5)
      character keyword*20
      character arrayz*80
      integer kk
      data newcount, linecount /0, 0/


      if (linecount .eq. 0) oldcount = 0


c  IF DOING MULTIPLE RUNS: if this is not the first reading of the
c  parameter file, then read down to correct place in parameter file, 
c  using "linecount", and then skip the re-initialization of the 
c  various variables
cc      rewind nfparam
cc      read (nfparam,1001,end=100) arrayz
cc      if (linecount .ne. 0) then
cc         do i=1,linecount
cc            read (nfparam,1001,end=100) arrayz
cc         enddo
cc         go to 4
cc      endif


c  INITIALIZE SOME VARIABLES: output file names and output file numbers
c  f1out, nf1out:    verbose standard output
c  f2out, nf2out:    raw synthetic spectra, or summary abundances, or
c                    summary curves-of-growth
c  f3out, nf3out:    smoothed synthetic spectra
c  f4out, nf4out:    IRAF-text-style synthetic spectra output
c  f5out, nf5out:    postscript plot output
c  f6out, nf6out:    synthetic/observed comparison text output (gridsyn)
c  f7out, nf7out:    raw synthetic spectrum of primary of a binary, 
c                    output summary run table, or
c                    to be used for some other purpose as needed; at
c                    present not given keyword input
c  f8out, nf8out:    raw synthetic spectrum of secondary of a binary,
c                    or the kept lines in a line weedout, or for
c                    some other purpose as needed
c  f9out, nf9out:    combined smoothed synthetic spectrum of a binary star,
c                    or the discarded lines in a line weedout, or the
c                    mean raw synthesis for a population of stars, or
c                    some other purpose as needed
c                    to be used for some other purpose as needed`
c  f10out, nf10out:  combined smoothed synthetic spectrum of a binary star,
c                    to be used for some other purpose as needed`
      f1out =   'no_filename_given'
      f2out =   'no_filename_given'
      f3out =   'no_filename_given'
      f4out =   'no_filename_given'
      f5out =   'optional_output_file'
      f6out =   'optional_output_file'
      f7out =   'no_filename_given'
      f8out =   'no_filename_given'
      f9out =   'no_filename_given'
      f10out =  'no_filename_given'
      nf1out =   0
      nf2out =   0
      nf3out =   0
      nf4out =   0
      nf5out =   0
      nf6out =   0
      nf7out =   0
      nf8out =   0
      nf9out =   0
      nf10out =  0
      modelnum = 0


c  INITIALIZE SOME VARIABLES: input file names and input file numbers
      fmodel =  'no_filename_given'
      flines =  'no_filename_given'
      fslines = 'no_filename_given'
      fobs =    'no_filename_given'
      ftable =  'no_filename_given' 
      nfmodel =  0 
      nflines =  0
      nfslines = 0
      nfobs =    0
      nftable =  0


c  INITIALIZE SOME VARIABLES: terminal type;
c  set smterm = ' ' for the plotting package 'sm'
      smterm = 'x11'
      

c  INITIALIZE SOME VARIABLES: 
c  atmosphere data printing:              modprintopt  [old ipr(1)]
c  molecular equilibrium:                 molopt       [old ipr(2)]
c  line data:                             linprintopt  [old ipr(3)] 
c  flux/intensity                         fluxintopt   [old ipr(4)]
c  synthesis/single line computations:    [deleted]    [old ipr(5)] 
c  line formation level                   [deleted]    [old ipr(6)]
c  plotting:                              plotopt      [old ipr(7)] 
c  damping                                dampingopt   [old ipr(8)]
c  observed spectrum file type            specfileopt  [old ipr(9)]
c  format of linelist (formatted/not)     linfileopt
c  source function with scat+abs          scatopt
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
 

c  INITIALIZE SOME VARIABLES:
c   if fudge is less than or equal to 0 then it does not scale it...make it
c   -1.0 just to be sure we dont have floating point problems with 0.0
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


c  INITIALIZE SOME VARIABLES: spectroscopic binary parameters
      deltaradvel = 0.
      lumratio = 1.


c  INITIALIZE SOME VARIABLES: elements with special abundance data; for
c  practical reasons, keyword "abundances" and "isotopes" must be specified
c  in each RUN; they must be reset each time.
4     neq = 0
      numpecatom = 0
      numatomsyn = 0
      newnumpecatom = 0
      newnumatomsyn = 0
      ninetynineflag = 0
      do i=1,95
         pec(i) = 0
         newpec(i) = 0
         do j=1,5
            abfactor(j) = 0.
            pecabund(i,j) = 0.
            newpecabund(i,j) = 0.
         enddo
      enddo
      numiso = 0
      numisosyn = 0
      newnumiso = 0
      newnumisosyn = 0
      do i=1,20
         isotope(i) = 0.
         newisotope(i) = 0.
         do j=1,5
            isoabund(i,j) = 0.
            newisoabund(i,j) = 0.
         enddo
      enddo


c  read a line of the parameter file into "arrayz"; decode it to get the 
c  key word, which is always within the first 20 characters;  dump the 
c  rest of arrayz into "array"
c 5     write (array,1007)
c      read (nfparam,1001,end=98) arrayz
c      linecount = linecount + 1
c      i=index(arrayz,' ')
c      keyword = arrayz(1:i-1)
c      array = arrayz(i:)


c  keyword 'RUN' signals that there are either multiple syntheses being
c  done or multiple comparisons with observed spectra
c      if (keyword .eq. 'RUN') then
c         read (array,*) newcount
c         if (newcount .gt. oldcount+1) then
c            linecount = linecount - 1
c            oldcount = syncount
c            go to 100
c         else
c            syncount = newcount
c            go to 5
c         endif
c      endif


c  keyword 'freeform' indicates whether or not the linelist will be read
c  in under format control (7e10.3) or will be free-form.  If freeform = 0,
c  the default value, then the old-style formatted input will be used;
c  If freeform = 1, unformatted read will be used, BUT the user must then
c  give values for all quantities (that is, explicit zeros will need to
c  be put instead of blank spaces.
c      if     (keyword .eq. 'freeform') then
c         read (array,*) linfileopt
 
 
c  keyword 'standard_out' controls the name of the verbose standard output
c      elseif (keyword .eq. 'standard_out') then
c         read (array,*) f1out


c  keyword 'summary_out' controls the name of either the EW summary or
c  the raw synthesis output
c      elseif (keyword .eq. 'summary_out') then
c         read (array,*) f2out


c  keyword 'hardpost_out' controls the name of a postscript plot output
c      elseif (keyword .eq. 'hardpost_out') then
c         read (array,*) f5out


c  keyword 'speccomp_out' controls the name of a text file containing the
c  comparisons (wavelength shifts, sigmas, etc.) between observed and
c  synthetic spectra
c      elseif (keyword .eq. 'speccomp_out') then
c         read (array,*) f6out


c  keyword 'bin_raw_out' controls the name of a file containing the
c  raw synthesis of a spectroscopic binary, with an appropriate velocity 
c  difference and luminosity ratio dialed in
c      elseif (keyword .eq. 'bin_raw_out') then
c         read (array,*) f9out


c  keyword 'bin_smo_out' controls the name of a file containing the
c  smoothed synthesis of a spectroscopic binary
c      elseif (keyword .eq. 'bin_smo_out') then
c         read (array,*) f10out

 
c  keyword 'summary_in' controls the name of the raw synthesis file,
c  created previously, that will be read in for plotting purposes
c      elseif (keyword .eq. 'summary_in') then
c         read (array,*) f2out


c  keyword 'smoothed_out' controls the name of the smoothed synthesis output
c      elseif (keyword .eq. 'smoothed_out') then
c         read (array,*) f3out


c  keyword 'keeplines_out' controls the name of the list of kept lines
c  for future synthetic spectrum runs
c      elseif (keyword .eq. 'keeplines_out') then
c         read (array,*) f8out


c  keyword 'tosslines_out' controls the name of the list of discarded lines
c  that are too weak to keep in future synthetic spectrum runs
c      elseif (keyword .eq. 'tosslines_out') then
c         read (array,*) f9out


c  keyword 'iraf_out' controls the name of the optional IRAF output
c      elseif (keyword .eq. 'iraf_out') then
c         read (array,*) f4out


c  keyword 'model_in' controls the name of input model atmosphere file
c      elseif (keyword .eq. 'model_in') then
c         read (array,*) fmodel


c  keyword 'lines_in' controls the name of the input line list
c      elseif (keyword .eq. 'lines_in') then
c         read (array,*) flines


c  keyword 'stronglines_in' controls the name of the input strong line list
c      elseif (keyword .eq. 'stronglines_in') then
c         read (array,*) fslines


c  keyword 'observed_in' controls the name of the input observed spectrum
c      elseif (keyword .eq. 'observed_in') then
c         read (array,*) fobs


c  keyword 'table_in' controls the name of the extra input instruction file
c      elseif (keyword .eq. 'table_in   ') then
c         read (array,*) ftable


c  keyword 'table_out' controls the name of the extra input instruction file
c      elseif (keyword .eq. 'table_out  ') then
c         read (array,*) f7out


c  keyword 'popsyn_out' controls the name of the extra input instruction file
c      elseif (keyword .eq. 'popsyn_out ') then
c         read (array,*) f9out


c  keyword 'rawbin_out ' controls the name of the input observed spectrum
c      elseif (keyword .eq. 'rawbin_out ') then
c         read (array,*) f9out


c  keyword 'smoobin_out' controls the name of the input observed spectrum
c      elseif (keyword .eq. 'smoobin_out') then
c         read (array,*) f10out


c  keyword 'atmosphere' controls the output of atmosphere quantities
c           0 = do not print out the atmosphere
c           1 = print out the standard things about an atmsophere
c           2 = print standard things and additional stuff like continuous
c                    opacities, etc.
c      elseif (keyword .eq. 'atmosphere') then
c         read (array,*) modprintopt


c  keyword 'molecules ' controls the molecular equilibrium calculations
c           0 = do not do molecular equilibrium
c           1 = do molecular equilibrium but do not print results
c           2 = do molecular equilibrium and print results
c      elseif (keyword .eq. 'molecules') then
c         read (array,*) molopt
c         if     (molopt .eq. 0) then
c            nchars = 64
c            write (array,1009) 
c            call getasci (nchars,l0)
c            if (chinfo(1:1) .eq. 'n') then
c               stop
c            else
c               molopt = 1
c            endif
c         endif


c  keyword 'molset' controls the choice of which set of molecules will be
c  used in molecular equilibrium calculations.
c          1 = the small set involving H, C, N, O, Mg, Ti (DEFAULT)
c          2 = the large set more useful for very cool stars
c      elseif (keyword .eq. 'molset') then
c         read (array,*) molset


c  keyword 'deviations' controls whether, for synthetic spectrum computations,
c  an 'obs-comp' plot will be made in addition to the normal spectrum plot
c           0 = do not plot the obs-comp plot
c           1 = plot the obs-comp plot
c      elseif (keyword .eq. 'deviations') then
c         read (array,*) deviations


c  keyword 'lines     ' controls the output of line data
c           0 = print out nothing about the input lines
c           1 = print out standard information about the input line list
c           2 = gory line data print (usually for diagnostic purposes)
c      elseif (keyword .eq. 'lines') then
c         read (array,*) linprintopt
c         linprintalt = linprintopt


c  keyword 'gfstyle   ' controls the output of line data
c           0 = base-10 logarithms of the gf values (DEFAULT)
c           1 = straight gf values
c      elseif (keyword .eq. 'gfstyle') then
c         read (array,*) gfstyle


c  keyword 'contnorm  ' allows multiplicative adjustment of the
c           continuum; useful probably only for batch syntheses
c           the numbers employed should be around 1.0;
c           default is 1.000000
c      elseif (keyword .eq. 'contnorm') then
c         read (array,*) contnorm


c  keyword 'plotpars  ' allows you to set all of the plotting 
c     parameters if you know them in advance
c     0 = none set (default); user can change in plotting routine
c     1 = given in following lines as follows 
c xlow         xhi         ylo       yhi
c vshift       lamshift    obsadd    obsmult
c smooth-type  FWHM-Gauss  vsini     L.D.C.    FWHM-Macro     FWHM-Loren
c      elseif (keyword .eq. 'plotpars') then
c         read (array,*) iscale
c         if (iscale .ne. 0) then
c            read (nfparam,*) xlo, xhi, ylo, yhi
c            linecount = linecount + 1
c            read (nfparam,*) veladd, xadd, yadd, ymult
c            if (xadd .ne. 0.) then
c               veladd = 3.0d5*xadd/((xlo+xhi)/2.)
c               xadd = 0.
c            endif
c            linecount = linecount + 1
c            read (nfparam,*) smtype, fwhmgauss, vsini, limbdark, vmac,
c     .                      fwhmloren
c            linecount = linecount + 1
c         endif


c keyword 'trudamp     '  should moog use the detailed line damping for
c                         those transitions that have information stored in
c                         subroutine trudamp? (Default is *no*)
c      elseif (keyword .eq. 'trudamp') then
c         read (array,*) itru


c keyword 'veladjust   '  shoud moog try to do a cross-correlation between
c                         observed and synthetic spectra and use that to
c                         align the spectra better in wavelength
c                         (Default is *no*)
c      elseif (keyword .eq. 'veladjust') then
c         read (array,*) maxshift


c keyword 'units      ' controls the units in which moog 
c          outputs the final spectrum
c            0 = angs
c            1 = microns
c            2 = 1/cm
c      elseif (keyword .eq. 'units') then
c         read (array,*) iunits
c         if (iunits .ne. 0) then
c            write (*,1010)
c            stop
c         endif


c keyword 'iraf       ' allows the user to output a raw spectrum in 
c          a form suitable for IRAF's rtext input command
c            0 = don't do this, make output the normal way.
c            1 = make an IRAF-compatible output
c      elseif (keyword .eq. 'iraf') then
c         read (array,*) iraf


c  keyword 'scat       'allows the user to employ a source function
c          which has both scattering and absorption components     
c          0 = NO scattering
c          1 = scattering
c      elseif (keyword .eq. 'scat') then
c         read (array,*) scatopt


c  keyword 'flux/int  ' choses integrated flux or central intensity
c           0 = integrated flux calculations
c           1 = central intensity calculations
c      elseif (keyword .eq. 'flux/int') then
c         read (array,*) fluxintopt


c*****here are the calculations to set up the damping; for atomic lines
c     there are several options:
c        dampingopt = 0 and dampnum < 0 --->
c                             gammav = 10^{dampnum(i)}*(T/10000K)^0.3*n_HI
c        dampingopt = 0 and dampnum = 0 --->
c                             c6 = Unsold formula
c        dampingopt = 0 and dampnum > 10^(-10) --->
c                             c6 =  (Unsold formula)*dampnum(i)
c        dampingopt = 0 and dampnum(i) < 10^(-10) --->
c                             c6 = dampnum(i)
c        dampingopt = 1 --->
c                             gammav = gamma_Barklem if possible,
c                                        otherwise use dampingopt=0 options
c        dampingopt = 2 --->
c                             c6 = c6_Blackwell-group
c        dampingopt = 3 and dampnum <= 10^(-10) --->
c                             c6 = c6_NEXTGEN for H I, He I, H2
c        dampingopt = 3 and dampnum > 10^(-10) --->
c                             c6 = (c6_NEXTGEN for H I, He I, H2)*dampnum
c     for molecular lines (lacking a better idea) --->
c                                        c6 done as in dampingopt = 0
c      elseif (keyword .eq. 'damping') then
c         read (array,*) dampingopt


c  keyword 'obspectrum' controls the file type of the observed spectrum
c           0 = no observed spectrum is to be input
c           1 = read a true FITS file with internal read statements
c          -1 = as if obspectrum = 1, but on a byte-swapping machine
c           2 = (not implemented yet)
c           3 = read a true Fits file with the FITSIO package
c           4 = (not implemented yet)
c           5 = read a special MONGO style (wavelength, flux pair) file
c      elseif (keyword .eq. 'obspectrum') then
c         read (array,*) specfileopt
c         if (specfileopt .lt. 0) then
c            byteswap = 1
c            specfileopt = iabs(specfileopt)
c         endif


c   keyword 'histogram' makes histogram plots of observed spectra if
c   histoyes = 1
c      elseif (keyword .eq. 'histogram') then
c         read (array,*) histoyes


c  keyword 'terminal  ' gives the sm plotting window type
c           smterm = a character string of the sm window type (see the 
c           appropriate sm manual for a list)
c      elseif (keyword .eq. 'terminal') then
c         read (array,*) smterm


c  keyword 'plot      ' decides whether or not to make a plot of results
c           0 = do not make a plot
c           For syntheses: 1 = plot only synthetic spectra
c                          2 = plot synthetic and observed spectra
c                          3 = smooth the syntheses but don't plot
c           For line analyses: # = the minimum number of lines of a 
c                                  species necessary to trigger a plot
c           For curves-of-growth: 1 = make plots
c           For flux curves: 1 = make plots
c      elseif (keyword .eq. 'plot') then
c         read (array,*) plotopt


c  keyword 'abundances' gives the changes to be applied to the abundances
c           # = the number of different syntheses to run
c               (the next line gives the different abundance factors
c               to use)
c  minimum error check:  numatomsyn must equal numisosyn or code will stop

c      elseif (keyword .eq. 'abundances') then
c         neq = 0
c         numpecatom = 0
c         numatomsyn = 0
c         newnumpecatom = 0
c         newnumatomsyn = 0
c         ninetynineflag = 0
c         do i=1,95
c            pec(i) = 0
c            newpec(i) = 0
c            do j=1,5
c               abfactor(j) = 0.
c               pecabund(i,j) = 0.
c               newpecabund(i,j) = 0.
c            enddo
c         enddo
c         read (array,*) numpecatom,numatomsyn
c         if (numisosyn .ne. 0) then
c            if (numatomsyn .ne. numisosyn) then
cc               write (array,1002) numatomsyn, numisosyn
c               print *, "numatomsyn != numisosyn. stahp"
cc               call putasci (77,6)
c               stop
c            endif
c         endif
c         do l=1,numpecatom
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
c         if (numpecatom.eq.1 .and. jatom.eq.99) ninetynineflag = 1


c keyword 'isotopes   ' gives the isotopes used in the line list and their
c                       abundance relative to the parent spiecies
c  minimum error check:  numatomsyn must equal numisosyn or code will stop
c      elseif (keyword .eq. 'isotopes') then
c         numiso = 0
c         numisosyn = 0
c         newnumiso = 0
c         newnumisosyn = 0
c         do i=1,20
c            isotope(i) = 0.
c            newisotope(i) = 0.
c            do j=1,5
c               isoabund(i,j) = 0.
c               newisoabund(i,j) = 0.
c            enddo
c         enddo
c         read (array,*) numiso,numisosyn
c         if (numatomsyn .ne. 0) then
c            if (numatomsyn .ne. numisosyn) then
c               write (array,1002) numatomsyn, numisosyn
c               call putasci (77,6)
c               stop
c            endif
c         endif
c         do  j=1,numiso
c            read (nfparam,*) isotope(j),(isoabund(j,kk),kk=1,numisosyn)
c            linecount = linecount + 1
c         enddo


c  keyword 'lumratio' gives the ratio of the luminosity of two stars at a
c                     specific wavelength in a binary star system (used 
c                     only with driver "binary")
c      elseif (keyword .eq. 'lumratio') then
c         read (array,*) lumratio


c  keyword 'deltaradvel' gives the velocity difference between the stars
c                        binary star system (used only with driver "binary")
c      elseif (keyword .eq. 'deltaradvel') then
c         read (array,*) deltaradvel


c  keyword 'synlimits ' gives the wavelength parameters for syntheses;
c                       start and sstop are beginning and ending
c                       wavelengths, step is the step size in the
c                       syntheses, and delta is the wavelength range
c                       to either side of a synthesis point to consider
c                       for line opacity calculations
c      elseif (keyword .eq. 'synlimits') then
c         read (nfparam,*) start, sstop, step, delta
c         oldstart = start
c         oldstop  = sstop
c         oldstep  = step
c         olddelta = delta
c         step1000 = 1000.*step
c         if (dble(idnint(step1000))-step1000 .ne. 0.) then
c            write (*,1008) step
c            stop
c         endif
c         linecount = linecount + 1


c  keyword 'fluxlimits' gives the wavelength parameters for flux curves;
c                       start and sstop are beginning and ending
c                       wavelengths, and step is the step size in the
c                        flux curve
c      elseif (keyword .eq. 'fluxlimits') then
c         read (nfparam,*) start, sstop, step
c         linecount = linecount + 1


c  keyword 'blenlimits' gives the parameters for blended line abundance
c                       matches.  delwave is the wavelength offset 
c                       to the blue of first and to the red of the 
c                       last line in the blend to extend the syntheses; 
c                       step is the wavelength step size in the 
c                       computations; cogatom is the name of the
c                       element whose abundance should be varied
c                       to achieve an EW match with observations.
c      elseif (keyword .eq. 'blenlimits') then
c         read (nfparam,*) delwave, step, cogatom
c         linecount = linecount + 1


c  keyword 'coglimits ' gives the log(W/lambda) limits for curves-of-growth
c                       rwlow and rwhigh are the beginning
c                       and ending points of the log(red.width) values,
c                       rwstep is the step in log(red.width),
c                       cogatom is the declaration of which element
c                       will have its abundance varied (necessary only
c                       for spectrum synthesis curves-of-growth,
c                       and wavestep is a forced (if desired) step size
c                       in wavelength along the line (this applies to
c                       single line computations only
c      elseif (keyword .eq. 'coglimits') then
c         read (nfparam,*) rwlow, rwhigh, rwstep, wavestep, cogatom
c         linecount = linecount + 1


c  keyword 'limits    ' old limits format...tell the user to change the
c                       keyword and quit.

c      elseif (keyword .eq. 'limits') then
c      write(*,*) 'Warning: keyword changed to *synlimits*, *coglimits*'
c      write(*,*) 'for Syntesis and COG calculations.'
c      write(*,*) 'Here are the proper formats:'
c      write(*,*) 
c      write(*,*) 'synlimits <start> <stop> <step> <delta>'
c      write(*,*) 'coglimits <rwlow> <rwhigh> <rwstep> <wavestep> ',
c     .           '<cogatom>'
c         stop


c   keyword of strong for lines which are to be considered for all of the 
c   synthesis
c      elseif (keyword .eq. 'strong') then
c         read (array,*) dostrong
     

c  keyword word of opacit which takes the continuus opacity and scales it
c  with the form of kaplam(i)= kaplam(i)*((factor*10000)/t(i))
c  in Opacit.f after it calulates the normal kaplam
c  if value is <= 0 then it does not do it
c      elseif (keyword .eq. 'opacit') then
c         read (array,*) fudge
     
      
c  any other keyword causes great grudge
c      else
c         write (array,1006) keyword
c         call prinfo (5)
c         stop
c
c      endif
 

c  loop back to get another parameter
c      go to 5
 

c  wrap things up with a few assignments
c  98    if (control.eq.'gridsyn' .or. control.eq.'gridplo' .or.
c     .    control.eq.'binary ' .or. control.eq.'abandy ') then
c         control = 'gridend'
c      endif


c  assign plotting window type; if no type has been given in the
c  parameter file, then ask for it
c 100   if (smterm .eq. ' ') then
c         array = 'GIVE THE SM TERMINAL NAME : '
c         nchar = 28
c         call getasci (nchar,12)
c         smterm = chinfo(1:nchar)
c         ivstat = ivcleof(12,1)
c      endif
c      if (smterm.eq.'x11' .or. smterm.eq.'X11') then
c         if    (control .eq. 'synth  ' .or.
c     .          control .eq. 'synpop ' .or.
c     .          control .eq. 'synplot' .or.
c     .          control .eq. 'isoplot' .or.
c     .          control .eq. 'gridsyn' .or.
c     .          control .eq. 'gridplo' .or.
c     .          control .eq. 'doflux ' .or.
c     .          control .eq. 'cogsyn ' .or.
c     .          control .eq. 'cog    ' .or.
c     .          control .eq. 'isotop ' .or.
c     .          control .eq. 'binary ') then
c             smterm = smt1
c         else
c             smterm = smt2
c         endif
c      endif


c  for syntheses, store the plotting parameters
c      if (control.eq.'synth  ' .or. control.eq.'synplot' .or.
c     .    control.eq.'gridsyn' .or. control.eq.'gridplo' .or.
c     .    control.eq.'binary ' .or. control.eq.'synpop ') then
c         if (oldstart .eq. 0) then
c            write (*,1011) 
c            stop
c         endif
c         if (iscale .eq. 0) call plotremember (0)
c         call plotremember (1)
c      endif


c*****exit normally
      return


c*****format statements
c 1001  format (a80)
c 1002  format ('# OF ABUNDANCE (',i1,') AND ISOTOPIC (',i1,')',
c     .        ' SYNTHESES DO NOT AGREE!   I QUIT!       ')
c 1006  format ('THIS OPTION IS UNKNOWN TO MOOG: ', a10, ' I QUIT!')
c 1007  format (79(' '))
c 1008  format ('step =', f10.5, 'A but it cannot be more precise than ',
c     .        'the nearest 0.001A; I QUIT!')
c 1009  format ('WARNING: molecular eq. always done if ',
c     .        'Teff < 8000K;  OK (y/n)???')
c 1010  format ('the units=1 option is fragile; rerun with only ',
c     .       'Angstroms and units=0')
c 1011  format ('SYNTHESIS START NOT SET; I QUIT!')

      end






