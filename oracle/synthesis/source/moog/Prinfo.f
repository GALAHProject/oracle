
      subroutine prinfo (line)
c******************************************************************************
c     this routine prints out information on the text terminal
c******************************************************************************
 
      include 'Pstuff.com'
      include 'Atmos.com'
      character errm1*80,array1*80


c*****do not print out information in real time if the code is in 
c     batch mode
      if (silent .eq. 'y') return

      print *, "we shouldn't be here"
      return
       
c      if (line .eq. 1) then
c         istat = ivcleof(4,1)
c      endif
c
c      if (line .eq. maxline-5) then
c         errm1 = errmess
c         array1 = array
c10       array = 'WANT TO SEE MORE ([y]/n)? '
c         nchars = 26
c         call getasci (nchars,4+line)
c         if (chinfo(1:1).eq.'y' .or. nchars.le.0) then
c            istat = ivcleof(4,1)
c            line = 1
c            array = array1
c            errmess = errm1
c         elseif (chinfo(1:1) .eq. 'n') then
c            errmess = 'stopinfo!'
c            return
c         else
c            go to 10
c         endif
c      endif


c      istat = ivwrite(4+line,1,array,79)
c      errmess = 'everything OK'
c      return


      end












