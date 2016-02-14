	subroutine grbeamfep (nc,np,nc1,ncn,ip,np1,
     *                    m,mr,
     *                    x1,x2,s1,s,s2,
     *                    x,bm,bs,
     *                    cr,bar,disp,ier,
     *                    xaw,sw1,sw,
     *                    iswit,
     *                    ya,af,n1,n2,iopt)

c****************************************************************
c	multichannel data normalizing, 
c	input, output in multiplex form
c	beam computations:
c	1	choosing of channels numbers (ncn)
c	2	calculating mean and variance of each channel
c	3	computation of beam, calculation of it's mean and variance
c	4	whitenning of beam with preeviously estimated parameters
c	5	filtration of test signal by the whitennig filter
c	6	detection of wave phase at whitening beam trace
c******************************************************************
c	           programmers: v.i.pinsky,a.f.kushnir
c****************************************************************

       real x(*),x1(*),x2(*),s1(*),bar(*),xaw(*)
       real af(*),ya(*),cr(*),s(*)
       real wtd(500),cwt(50)
       integer ncn(*)

c**************** parameters *************************************
c  nc - input number of channels
c  np - input number of data points
c  nc1 - output number of channels
c  ncn(nc1) - array of choosed numbers of channels
c  ip -  initial point of processed data in input array
c  np1 - output number of data points, np1<np-ip-1
c  x1(nc*np) - input data (from 'rzm' programm)
c  x2(nc1*np1) - output data (for 'wzknor','fpwnor' programms)
c  s1(nc1) - mean values of channels
c  s(nc1) - variances of channels
c  s2 - average variance of channels
c  x(np1) - beam time series
c  bm - mean value of beam
c  bs - variance of beam
c  m - order of ar model for witening filtration
c  cr(m) - correlatin function of beam
c  af(m*(m+1)/2) - correlation matrix of beam
c  bar(m+1) - vector of autoregressive parameters of beam
c  disp - dispersion of ar residials
c  xaw(np1) - whitened beam time series
c  sw(np1) - whitened beam power
c  sw1 - mean value of whitened beam
c  iswit - new adaptation of whitening filter? yes=1,no=0
c  ya(m+1) - auxiliary array
c  mr - length of moving window for phase detection
c  n1 - beginning sample for noise definition
c  n2 - ending sample for noise definition
c  iopt - termination flag (1) after stack (2) after whitening (3) complete
c*********************************************************************
        s2=0.
       do 5 i=1,np1
5       x(i)=0.
c  calculation of mean and variance of each channel: s1(i),s(i)
c  computation of averaged channel (beam): x(j)
c  calculation of averaged variance of channels: s2
c
c  Begin Loop over channels
       do 10 i=1,nc1
        s1(i)=0.
        do 20 j=1,np1
          ij=i+(j-1)*nc1
c  x2 is partially demuxed x1 (data for channel i)
          x2(ij)=x1(ncn(i)+(j+ip-2)*nc)
          s1(i)=s1(i)+x2(ij)
c----check on indices
c	  write(*,*) ij, ncn(i)+(j+ip-2)*nc, x2(ij)
c	  read(*,*) pause
20     continue
c----s1 is now the average value for channel i
       s1(i)=s1(i)/np1
       do 30 j=1,np1
          ij=i+(j-1)*nc1
c-----demean x2
          x2(ij)=x2(ij)-s1(i)
c-----x now is a stack for all demeaned channels
          x(j)=x(j)+x2(ij)/nc1
c	  write(*,*) ij, j, x(j), x2(ij)
c	  read(*,*) pause
30     continue
       s(i)=0.
       do 60 j=1,np1
          ij=i+(j-1)*nc1
60     s(i)=s(i)+x2(ij)**2
c----s is now the variance of channel i
       s(i)=s(i)/np1
       s2=s2+s(i)
c   end of loop for channels
10     continue
c----s2 is the average variance over all channels
       s2=s2/nc1
c----if iopt=1, then just reassign the stack and return
	if (iopt.eq.1) then
	  do 61 j = 1, np1
	    xaw(j) = x(j)
61	  continue
	  return
	end if
c	write(*,*) ' iswit = ', iswit
c	open(1,file='grbeam.in')
c	write(1,*)(i, x(i),i=1,np1)
c	close (1)
c	write(*,*)'output data is written to file grbeam.in'
c******************************************************************
       if(iswit.eq.0)then
c    whitening of beam with preeviously estimated parameters
c          print *,'call grfilt'
         call tgrfilt(x,np1,bar,m,xaw)
      else
c    adaptation and whitening filtration of beam
c          print *,'call grdetest1'
c        call tgrdetest1(x,np1,m,ier,af,bar,disp,xaw,ya,cr)
        call tgrdetest3(x,np1,n1,n2,m,ier,af,bar,disp,xaw,ya,cr)
       end if
	if (iopt.eq.2) return
c	open(1,file='grbeam.white')
c	write(1,*)(i, xaw(i),i=1,np1)
c1000	format(i5, f15.5)
c	close (1)
c	write(*,*)'output data is written to file grbeam.white'
c*******************************************************************
c
c    calculation of beams mean values an variances
       bs=0.
       bm=0.
       sw=0.
       sw1=0.
       ssw=0.
       ssw1=0.
      do 70 j=1,np1
            xaw(j)=xaw(j)/sqrt(disp)
            sw=sw+xaw(j)**2
            sw1=sw1+xaw(j)
           bm=bm+x(j)
70         bs=bs+x(j)**2
       bm=bm/np1
       bs=bs/np1
       sw=sw/(np1-m)
       sw1=sw1/(np1-m)
	write(*,*) ' bm, bs, sw, sw1 = ', bm, bs, sw, sw1
c**********************************************************************
c
c detection of wave phases at the whitened beam trace
	mw=m
c	mr=50
	do 111 j=1,mr+mw
	wtd(j)=xaw(j)
111	continue
c
	do 112 j=mr+mw+1,np1
c
	do 113 jr=2,mr+mw
113	wtd(jr-1)=wtd(jr)
	wtd(mr+mw)=xaw(j)
c
	dd=1.
	call tgrwhdet(mr,mw,dd,wxid,wtd,cwt)
	xaw(j-int(mr/2))=wxid
c
112	continue
c
	ke=np1-int(mr/2)+1
	do 114 j=ke,np1
	xaw(j)=xaw(ke-1)
c	xaw(j)=0.
114	continue
	kb=mr+mw-int(mr/2)
	do 115 j=1,kb
	xaw(j)=xaw(kb+1)
c	xaw(j)=0.
115	continue
c*******************************************************************
c
c	open(1,file='grbeam.out')
c	do 1010 i = 1, np1
c	write(1,*) i, xaw(i)
c1010	continue
c	close (1)
c	write(*,*)'output data is written to file grbeam.out'

       return
       end
c
      subroutine tgrccov(x,n,y,c,m)
c************************************************************************
c  estimation of two time series crosscovariance function
c            programmer: i.v.savin
c*************************************************************************
c  x(n),y(n) - input time series
c  n - length data in points
c  c(m+1) - output crosscovarianse function
c  m - number of nonzerro lags
c
c  subroutines called:
c	none
c
c*************************************************************************
      real x(*),c(*),y(*)
      fn=n
      m1=m+1
      do 20 l=1,m1
         s=0.0
         do 10 i=l,n
            j=i-l+1
            s=s+x(i)*y(j)
   10    continue
         c(l)=s/fn
   20 continue
      return
      end
c
      subroutine tgrdetest1(x,n1,m,ier,a,b,d,e,y,cr)
c*********************************************************************
c   ar - adoptation using of the noise realisation  for the
c  "real time" detector of seismic events (for 1-component
c   record)
c              programmer: a.f.kushnir
c*********************************************************************
c   x -input record of noise (size > n1)
c   n1- length of the noise estimated
c   m - order of ar-model
c   b - output - vector of ar parameters of the noise (size m)
c   d - output - dispersion of the ar residials
c   a - output - fisher matrix of ar parameters (m+1)*(m+1) size
c   e - output - whitened noise - array of n1 size
c   y,cr - auxiliary arrays of m+1 size
c
c	subroutines called:
c		tgrsinv, tgrmprd, tgrfilt, tgrccov
c
c***********************************************************************
      real x(*),e(*),a(*),b(*),y(*),cr(*)
      data eps/1.0e-6/
c
      ier=0
      m1=m+1
c      write(6,*)'enter grccov'
	write(6,*) 'x(1), n1, y(1), m' ,x(1), n1, y(1), m
      call tgrccov(x,n1,x,y,m)
      do 5 i=1,m1
5          cr(i)=y(i)/y(1)
c      write(6,*)'fisher matrix',(a(i),i=1,m*(m+1)/2)
      do 10 i=1,m
      do 10 j=i,m
         iy=j-i+1
         ia=i+(j*j-j)/2
         a(ia)=y(iy)
   10 continue
	write(6,*)'inv. fisher.matrix',(a(i),i=1,m*(m+1)/2)
c      write(6,*)'enter grsinv'
      call tgrsinv(a,m,eps,ier0)
      if (ier0.lt.0) goto 100
      do 20 iy=1,m
20    y(iy)=y(iy+1)
c     write(6,*)'solution of Yule-Walker equations'
c      write(6,*)'enter grmprd'
      call tgrmprd(a,y,b,m,m,1,0,1)
	write(6,*)'output grmprd',(b(i),i=1,m)
      do 30 i=1,m
      b(m+2-i)=(-1)*b(m+1-i)
30    continue
      b(1)=1.
c     array x whitening filtration using ar-parameters estimated
c      write(6,*)'enter grfilt'
      call tgrfilt(x,n1,b,m,e)
	do 22 i = 1, 20
	  write(*,*) i, e(i)
22	continue
c     calculation of the innovation variance
c      write(6,*)'enter grccov'
      call tgrccov(e,n1,e,d,0)
      d1=d*n1
      d=d1/(n1-m)
      return
100   ier=-1
      return
      end
      subroutine tgrdetest2(x,n,n1,m,ier,a,b,d,e,y,cr)
c*********************************************************************
c   ar - adoptation using of the noise realisation  for the
c  "real time" detector of seismic events (for 1-component
c   record)
c              programmer: a.f.kushnir
c*********************************************************************
c   x -input record of noise (size > n1)
c   n - length of the time record (x)
c   n1- length of the record used for noise estimation
c   m - order of ar-model
c   b - output - vector of ar parameters of the noise (size m)
c   d - output - dispersion of the ar residials
c   a - output - fisher matrix of ar parameters (m+1)*(m+1) size
c   e - output - whitened trace - array of n size
c   y,cr - auxiliary arrays of m+1 size
c
c	subroutines called:
c		tgrsinv, tgrmprd, tgrfilt, tgrccov
c
c***********************************************************************
      real x(*),e(*),a(*),b(*),y(*),cr(*)
      data eps/1.0e-6/
c
      ier=0
      m1=m+1
c      write(6,*)'enter grccov'
c      write(6,*) 'x(1), n1, y(1), m' ,x(1), n1, y(1), m
      call tgrccov(x,n1,x,y,m)
      do 5 i=1,m1
5          cr(i)=y(i)/y(1)
c      write(6,*)'fisher matrix',(a(i),i=1,m*(m+1)/2)
      do 10 i=1,m
      do 10 j=i,m
         iy=j-i+1
         ia=i+(j*j-j)/2
         a(ia)=y(iy)
   10 continue
c      write(6,*)'inv. fisher.matrix',(a(i),i=1,m*(m+1)/2)
c      write(6,*)'enter grsinv'
      call tgrsinv(a,m,eps,ier0)
      if (ier0.lt.0) goto 100
      do 20 iy=1,m
20    y(iy)=y(iy+1)
c     write(6,*)'Solution of Yule-Walker equations'
c      write(6,*)'enter grmprd'
      call tgrmprd(a,y,b,m,m,1,0,1)
c     write(6,*)'output grmprd',(b(i),i=1,m)
      do 30 i=1,m
      b(m+2-i)=(-1)*b(m+1-i)
30    continue
      b(1)=1.
c     array x whitening filtration using ar-parameters estimated
c      write(6,*)'enter grfilt'
c      call tgrfilt(x,n1,b,m,e)
      call tgrfilt(x,n,b,m,e)
c     calculation of the innovation variance
c      write(6,*)'enter grccov'
      call tgrccov(e,n1,e,d,0)
      d1=d*n1
      d=d1/(n1-m)
      return
100   ier=-1
      return
      end
      subroutine tgrdetest3(x,n,n1,n2,m,ier,a,b,d,e,y,cr)
c*********************************************************************
c   ar - adoptation using of the noise realisation  for the
c  "real time" detector of seismic events (for 1-component
c   record)
c              programmer: a.f.kushnir
c*********************************************************************
c   x -input record of noise (size > n1)
c   xaux - auxilliary array for holding noise
c   n - length of the time record (x)
c   n1- beginning of part of record used for noise estimation
c   n2- end of part of record used for noise estimation
c   m - order of ar-model
c   b - output - vector of ar parameters of the noise (size m)
c   d - output - dispersion of the ar residials
c   a - output - fisher matrix of ar parameters (m+1)*(m+1) size
c   e - output - whitened trace - array of n size
c   y,cr - auxiliary arrays of m+1 size
c
c	subroutines called:
c		tgrsinv, tgrmprd, tgrfilt, tgrccov
c
c***********************************************************************
      real x(*),e(*),a(*),b(*),y(*),cr(*)
      data eps/1.0e-6/
c
      ier=0
      m1=m+1
c      write(6,*)'enter grccov'
	noi = n2 - n1 + 1
c	write(6,*) 'x(1), n1, n2, noi, n, m' ,x(1), n1, n2, noi, n, m
c----y is output as the cross-covariance function
      call tgrccov(x(n1),noi,x(n1),y,m)
c      call tgrccov(x,n,x,y,m)
c----cr is the cross-covariance normalized by the square of x (=y(1))
      do 5 i=1,m1
5          cr(i)=y(i)/y(1)
c      write(6,*)'fisher matrix',(a(i),i=1,m*(m+1)/2)
c---a() is a lower triangular cross-covariance matrix
c
c		y(1)
c	a = 	y(2) y(1)
c		y(3) y(2) y(1)
c

      do 10 i=1,m
      do 10 j=i,m
         iy=j-i+1
         ia=i+(j*j-j)/2
         a(ia)=y(iy)
   10 continue
c	write(6,*)'inv. fisher.matrix',(a(i),i=1,m*(m+1)/2)
c      write(6,*)'enter grsinv'
c---most likely this inverts a
      call tgrsinv(a,m,eps,ier0)
      if (ier0.lt.0) goto 100
      do 20 iy=1,m
20    y(iy)=y(iy+1)
c     write(6,*)'Solution of Yule-Walker equations'
c      write(6,*)'enter grmprd'
c----b is a-1*y
      call tgrmprd(a,y,b,m,m,1,0,1)
c	write(6,*)'output grmprd',(b(i),i=1,m)
c----advance the b index by 1
      do 30 i=1,m
      b(m+2-i)=(-1)*b(m+1-i)
30    continue
      b(1)=1.
c     array x whitening filtration using ar-parameters estimated
c      write(6,*)'enter grfilt'
c      call tgrfilt(x,n1,b,m,e)
c----apply the b (whitening) filter to x; output is e
      call tgrfilt(x,n,b,m,e)
c	do 22 i = 1, 20
c	  write(*,*) i, e(i)
c22	continue
c     calculation of the innovation variance
c      write(6,*)'enter grccov'
      call tgrccov(e(n1),noi,e(n1),d,0)
c      call tgrccov(e,n,e,d,0)
      d1=d*noi
      d=d1/(noi-m)
      return
100   ier=-1
      return
      end
c
      subroutine tgrfilt(x,n,b,m,e)
c***************************************************************
c   filtration of data by causal filter
c       programmer: a.f.kushnir
c******************************************************************
c  x(n) - input time series
c  n - data length
c  b(m+1) - filter impulse responce
c  m+1 - length of filter imp. responce
c  e(n) - output filtered data
c
c  subroutines called:
c	none
c
c*******************************************************************
      real x(*),b(*),e(*)
      m1=m+1
      ie=0
      do 5 i=1,m1
5     e(i)=0.0
      do 20 jx=m1,n
         w=0.0
         do 10 ib=1,m1
            ix=jx-ib+1
            w=w+b(ib)*x(ix)
10          continue
         ie=ie+1
         e(ie+m)=w
20    continue
      return
      end
c
      subroutine tgrloc(i,j,ir,n,m,ms)
c-------------------------------------
c  subroutines called:
c	none
c-------------------------------------
      ix=i
      jx=j
      if (ms-1) 10,20,30
   10 irx=n*(jx-1)+ix
      go to 36
   20 if (ix-jx) 22,24,24
   22 irx=ix+(jx*jx-jx)/2
      go to 36
   24 irx=jx+(ix*ix-ix)/2
      go to 36
   30 irx=0
      if (ix-jx) 36,32,36
   32 irx=ix
   36 ir=irx
      return
      end
c
      subroutine tgrmfsd(a,n,eps,ier)
      dimension a(*)
      double precision dpiv,dsum
      if(n-1) 12,1,1
    1 ier=0
      kpiv=0
      do 11 k=1,n
      kpiv=kpiv+k
      ind=kpiv
      lend=k-1
      tol=abs(eps*a(kpiv))
      do 11 i=k,n
      dsum=0.d0
      if(lend) 2,4,2
    2 do 3 l=1,lend
      lanf=kpiv-l
      lind=ind-l
    3 dsum=dsum+dble(a(lanf)*a(lind))
    4 dsum=dble(a(ind))-dsum
      if(i-k) 10,5,10
    5 if(sngl(dsum)-tol) 6,6,9
    6 if(dsum) 12,12,7
    7 if(ier) 8,8,9
    8 ier=k-1
    9 dpiv=dsqrt(dsum)
      a(kpiv)=dpiv
      dpiv=1.d0/dpiv
      go to 11
   10 a(ind)=dsum*dpiv
   11 ind=ind+i
      return
   12 ier=-1
      return
	end
c
      subroutine tgrmprd(a,b,r,n,m,msa,msb,l)
c--------------------------------------------
c  subroutines called:
c	tgrloc
c--------------------------------------------
      dimension a(*),b(*),r(*)
c         1
      ms=msa*10+msb
      if(ms-22) 30,10,30
   10 do 20 i=1,n
   20 r(i)=a(i)*b(i)
      return
c         2
   30 ir=1
      do 90 k=1,l
      do 90 j=1,n
      r(ir)=0.0
      do 80 i=1,m
      if(ms) 40,60,40
   40 call tgrloc(j,i,ia,n,m,msa)
      call tgrloc(i,k,ib,m,l,msb)
      if(ia) 50,80,50
   50 if(ib) 70,80,70
   60 ia=n*(i-1)+j
      ib=m*(k-1)+i
   70 r(ir)=r(ir)+a(ia)*b(ib)
   80 continue
   90 ir=ir+1
      return
      end
c
      subroutine tgrsinv(a,n,eps,ier)
c------------------------------------------
c	subroutines called:
c		tgrmfsd
c-------------------------------------------
      dimension a(*)
      double precision din,work
c         1
      call tgrmfsd(a,n,eps,ier)
      if(ier) 9,1,1
c         2
    1 ipiv=n*(n+1)/2
      ind=ipiv
c         3
      do 6 i=1,n
      din=1.d0/dble(a(ipiv))
      a(ipiv)=din
      min=n
      kend=i-1
      lanf=n-kend
      if(kend) 5,5,2
    2 j=ind
c         4
      do 4 k=1,kend
      work=0.d0
      min=min-1
      lhor=ipiv
      lver=j
c         5
      do 3 l=lanf,min
      lver=lver+1
      lhor=lhor+l
    3 work=work+dble(a(lver)*a(lhor))
c         6
      a(j)=-work*din
    4 j=j-min
c         7
    5 ipiv=ipiv-min
    6 ind=ind-1
c         8
c         9
c         10
c         11
      do 8 i=1,n
      ipiv=ipiv+i
      j=ipiv
c         12
      do 8 k=i,n
      work=0.d0
      lhor=j
c         13
      do 7 l=k,n
      lver=lhor+k-i
      work=work+dble(a(lhor)*a(lver))
    7 lhor=lhor+l
c         14
      a(j)=work
    8 j=j+k
c         15
    9 return
      end
c
      subroutine tgrsovd(x,n,y,c,m)
c--------------------------------------
c  estimation of cross-covariances
c  with norming at sqrt(n)
c  x(n)- first trace
c  y(n)- second trace
c  c(m+1)- cross-covariances
c
c  subroutines called:
c	none
c****************************************
      real x(*),c(*),y(*)
      fn=n
      m1=m+1
      do 20 l=1,m1
         s=0.0
         do 10 i=l,n
            j=i-l+1
            s=s+x(i)*y(j)
   10    continue
         c(l)=s/sqrt(fn)
   20 continue
      return
      end
c
c
      subroutine tgrwhdet(l,m,d,xisqt,e,cwt)
c----------------------------------------------------------------
c   chi-square-detector for a whitened trace
c----------------------------------------------------------------
c   l - length of moving window
c   m - order of ar model
c   d - dispersion of the ar model residials
c   xisqt - output, current value of the meansq after whitening
c   e - input- array of l+m size, moving window of the whitened trace
c   cwt- array of size m+1 curent vlue of the auto-correlations e
c
c  subroutines called:
c	tgrsovd
c
c****************************************************************
       real e(*),cwt(*)
       rl=l
c estimation of autocovariences of whitened data
       call tgrsovd(e,l,e,cwt,m)
       fstat=(cwt(1)-sqrt(rl)*d)**2/d**2
       w1=0.
       do 10 i=2,m+1
10     w1=w1+cwt(i)**2
       xisqt=w1/d**2+fstat
       return
       end
