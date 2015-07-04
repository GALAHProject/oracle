
      subroutine damping (linnumber)
c******************************************************************************
c     This subroutine computes damping 'gamma' factors 
c     and then Voigt parameters 'a'
c******************************************************************************

      implicit real*8 (a-h,o-z)
      real*8      sigma, alpha, gx, gammaf, vbar
      PARAMETER (PI=3.14159265)
      PARAMETER (A0=5.29177249E-11)

      include 'Atmos.com'
      include 'Linex.com'
      include 'Quants.com'


      j = linnumber
      iwave = int(wave1(j))
      iatom10 = nint(10.*atom1(j))
      if (dampnum(j) .lt. 0.) dampnum(j) = 10.**dampnum(j)


c*****for a few lines, explicit detailed broadening terms have 
c     appeared in the literature, and so do these lines with a 
c     sepaarate subroutine
      if (itru .eq. 0) then
c     Ca II
         if (iatom10 .eq. 201) then
            if (iwave .eq. 8498 .or.
     .          iwave .eq. 8542 .or.
     .          iwave .eq. 8662 .or.
     .          iwave .eq. 3933) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
c     CH
         elseif(iatom10 .eq. 1060) then
            if (iwave .eq. 3693) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
c     Ca I
         elseif (iatom10 .eq. 200) then
            if (iwave.eq.6717 .or. iwave.eq.6318
     .          .or. iwave.eq.6343 .or. iwave.eq.6361) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
c     Ca I autoionization
         elseif (iatom10 .eq. 200) then
            if (iwave.eq.6318 .or. 
     .          iwave.eq.6343 .or.
     .          iwave.eq.6361) then
               call trudamp (j)
               damptype(j) = 'TRUEgam'
               return
            endif
         endif
      endif
  

c*****here are the calculations to set up the damping; for atomic lines 
c     there are several options:
c        dampingopt = 0 and dampnum = 0 ---> 
c                             c6 = Unsold formula
c        dampingopt = 0 and dampnum < 10^(-15) ---> 
c                             c6 = dampnum
c        dampingopt = 0 and 10^(-15) < dampnum < 10^(-5) ---> 
c                             gamma = dampnum
c        dampingopt = 0 and dampnum(i) > 10^(-5) ---> 
c                             c6 =  (Unsold formula)*dampnum
c        dampingopt = 1 --->    
c                             gammav = gamma_Barklem if possible, 
c                                      otherwise use dampingopt=0 options
c        dampingopt = 2 ---> 
c                             c6 = c6_Blackwell-group
c        dampingopt = 3 and dampnum <= 10^(-10) --->             
c                             c6 = c6_NEXTGEN for H I, He I, H2
c        dampingopt = 3 and dampnum > 10^(-10) --->             
c                             c6 = (c6_NEXTGEN for H I, He I, H2)*dampnum
c     for molecular lines (lacking a better idea) --->
c                                        c6 done as in dampingopt = 0
c        dampingopt = 4 and dampnum > 20 
c                             c6 = int(sigma).alpha use ABO theory
c        dampingopt = 4 and dampnum < 20
c                             do as if dampingopt = 2



c*****these damping calculations are done at each atmosphere level
      if (debug .gt. 0) write (nf1out,1001) j, wave1(j)
      do i=1,ntau
         ich = nint(charge(j))
         v1 = dsqrt(2.1175d8*t(i)*(1.0/amass(j)+1.008))


c*****first calculate an Unsold approximation to gamma_VanderWaals
         if (atom1(j) .gt. 100.) then
            ebreakup = 7.0
         else 
            ebreakup = chi(j,ich)
         endif
         if (e(j,1).ge.ebreakup .or. e(j,2).ge.ebreakup) then
            unsold = 1.0e-33
         else
            unsold = dabs(1.61d-33*(13.598*charge(j)/(ebreakup -
     .               e(j,1)))**2 - 1.61d-33*(13.598*charge(j)/
     .               (ebreakup-e(j,2)))**2)
         endif


c*****dampingopt = 0 or 
c*****dampingopt = 1 and no Barklem data
         if     (dampingopt .eq. 0 .or.
     .          (dampingopt.eq.1 .and. gambark(j).lt.0)) then
            if     (dampnum(j) .eq. 0.0) then
               damptype(j) = 'UNSLDc6'
               gammav = 17.0*unsold**0.4*v1**0.6*numdens(1,1,i)
            elseif (dampnum(j) .lt. 1.0d-15) then
               damptype(j) = '   MYc6'
               gammav = 17.0*dampnum(j)**0.4*v1**0.6*numdens(1,1,i)
            elseif (dampnum(j) .lt. 1.0d-04) then
               damptype(j) = 'MYgamma'
               gammav = dampnum(j)*(t(i)/10000.)**0.3*numdens(1,1,i)
            else
               damptype(j) = 'MODUNc6'
               gammav = 
     .            17.0*(unsold*dampnum(j))**0.4*v1**0.6*numdens(1,1,i)
            endif


c*****dampingopt = 1 with extant Barklem data
         elseif ((dampingopt.eq.1 .and. gambark(j).gt.0.) .or. 
     .      (dampingopt .eq. 4 .and. dampnum(j) .lt. 20.)) then
               damptype(j) = 'BKgamma'
               gammav =
     .            gambark(j)*(t(i)/10000.)**alpbark(j)*numdens(1,1,i)


c*****dampingopt = 2
         elseif (dampingopt .eq. 2) then
            damptype(j) = 'BLKWLc6'
            gammav = 17.0*((1.0 + 0.67*e(j,1))*unsold)**0.4*
     .               v1**0.6*numdens(1,1,i)


c*****dampingopt = 3
         elseif (dampingopt .eq. 3) then
            damptype(j) = 'NXTGNc6'
            if (dampnum(j) .le. 1.0d-10) dampnum(j) = 1.0
                 c6h = dabs(1.01d-32*(charge(j)**2)*
     .              (13.598/(ebreakup - e(j,1)))**2 - 1.61d-33*
     .              (13.598/(ebreakup-e(j,2)))**2)
                 c6he = dabs((0.204956/0.666793)*1.01d-32*
     .                  (charge(j)**2)*(13.598/(ebreakup - 
     .                  e(j,1)))**2 - 1.61d-33*(13.598/(ebreakup-
     .                  e(j,2)))**2)
                 c6ht = dabs((0.806/0.666793)*1.01d-32*
     .                  (charge(j)**2)*(13.598/(ebreakup - 
     .                  e(j,1)))**2 - 1.61d-33*(13.598/(ebreakup-
     .                  e(j,2)))**2)
               gammav = 17.0*v1**0.6*(c6h**0.4*numdens(1,1,i) +
     .                  c6he**0.4*numdens(2,1,i) +
     .                  c6ht**0.4*numdens(8,1,i))*dampnum(j)**0.4


c*****dampingopt = 4
         elseif (dampingopt .eq. 4 .and. dampnum(j) .gt. 20) then
               damptype(j) = 'ABO'

               sigma = int(dampnum(j))
               alpha = dampnum(j) - int(dampnum(j))
               
               gx = (2.0 - alpha*0.5) - 1.0
               gammaf = 1 + (-0.5748646 + (0.9512363 + 
     .           (-0.6998588 +(0.4245549 - 0.1010678*gx)
     .           *gx)*gx)*gx)*gx

c    Compute the half-width per unit perturber number density for
c    this temperature.

               gammav = (4./PI)**(alpha*0.5)*gammaf*1.E4*sigma*A0*A0

c    K = 1.380658E-23 and M0 = 1.660540E-27
c    K/M = 8314.51214665109
               vbar = (8.*8314.51214665109*t(i)/PI
     .             *(1./1.008+1./xam(atom1(j))))**0.5
               gammav = gammav*((vbar/1.E4)**(1.-alpha))

c    While nhtot gives the number density of hydrogen, it seems
c    MOOG just approximates the number density of helium with
c    nhe = xabund(2)*nhtot(i) (see Trudamp.f)
               gammav = gammav*(nhtot(i) + 
     .            0.41336*xabund(2)*nhtot(i))*1.E6*2.0

         endif
         
c*****compute radiative broadening either by an approximate formula or 
c*****the value in Barklem.dat) 
         if (gamrad(j).ne.0.0 .and. dampingopt .eq. 1) then
            gammar = gamrad(j)
         else
            gammar = 2.223d15/wave1(j)**2
         endif
        

c*****now Stark broadening (approximate formulae)
         excdiff = chi(j,nint(charge(j))) - e(j,2)
         if (excdiff .gt. 0.0 .and. atom1(j).lt.100.) then
            effn2 = 13.6*charge(j)**2/excdiff
         else
            effn2 = 25.
         endif
         gammas = 1.0e-8*ne(i)*effn2**2.5


c*****now finish by summing the gammas and computing the Voigt *a* values
         gammatot = gammar + gammas + gammav
         a(j,i) = gammatot*wave1(j)*1.0d-8/(12.56636*dopp(j,i))
         if (debug .gt. 0) write (nf1out,1002) i, gammar, 
     .      gammas, gammav, gammatot, a(j,i)
      enddo
      return


c*****format statements
1001  format(//' LINE BROADENING PARAMETERS FOR LINE', i4,
     .       ' AT WAVELENGTH',f8.2/
     .       '  i',4x,'natural',6x,'Stark',4x,'VdWaals',
     .       6x,'total',5x,'a(j,i)')
1002  format (i3,1p5e11.3)


      end




