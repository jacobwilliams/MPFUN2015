program tpphix3

!   This computes the Poisson phi function, as described in paper by Bailey,
!   Borwein, Crandall and Zucker, and finds the associated minimal polynomial.
!   Three-level PSLQ version.

!   David H Bailey   19 Apr 2017

!   Parameters set in parameter statements below:
!     ipal  0: The palindromic technique (for even kd) is NOT implemented.
!           1: The palindromic technique (for even kd) IS implemented.
!           NOTE: If ipal = 1, the parameter n below = 1 + degree of half
!           polynomial. In other words, n = size of relation in call to PSLQM3.
!           The full output polynomial has degree 2*n-2.
!     lcx   Size of chr1 array; default = 64.
!           ***Must be a multiple of 64.  Must be <= lmx.
!           ***Must be >= 1 + size in digits of largest coefficient.
!     lmx   Length of line1; default = 131072.
!     mpi   Multiplier of pi on RHS of conjectured identity; default = 8.
!     kd    Denominator of arguments of phi_2; default = 30.
!     kp    Numerator of first argument; default = 1.
!     kq    Numerator of second argument; default = 1.

!   PSLQM3 parameters:
!     idb   Debug level (0 - 3); default = 2.
!     n     Integer relation vector length = 1 + polynomial degree; default = 33.
!     ndp   Full precision level in digits; default = 1100.
!           ***Must be <= mpipl in module MPFUNF.
!           ***Must be >= ndpm + precision level required to find relation.
!     ndpm  Intermediate precision level in digits; default = 160.
!            Increase ndpm if run aborts in MP initialization.
!     nep   Log10 of full precision epsilon for detections; default = 30 - ndp.
!           ***Must not be smaller than the accuracy of input data.  In other
!           words, if data is accurate to 10^(-200), then nep >= -200.
!     nepm  Log10 of intermediate precision epsilon; default = 20 - ndpm.
!     rb    Log10 of max size (Euclidean norm) of acceptable relation;
!           default = 100. Run will abort if this is exceeded.
!     nwds  Precision level in words; default = int (ndp / mpdpw + 2).
!     nwdsm Intermediate level in words; default = int (ndpm / mpdpw + 2).

use mpmodule
implicit none
real (8) rb
integer idb, ipal, lcx, lmx, mpi, nsq
parameter (idb = 2, ipal = 1, lcx = 64, lmx = 131072, mpi = 8, rb = 100.d0)
integer i, iq, i1, i2, j, j1, k, kd, kp, kq, l1, m, n, n1, &
  ndp, ndpm, nep, nepm, nwds, nwdsm
parameter (kd = 30, kp = 1, kq = 1, n = 33, ndpm = 160, ndp = 1100, &
  nep = 30 - ndp, nwds = int (ndp / mpdpw + 2), nepm = 20 - ndpm, &
  nwdsm = int (ndpm / mpdpw + 2))
real (8) d1, d2, second, tm, tm0, tm1, tm2, tm3
external second
type (mp_real) al, beta, eps, pi, qq, xx, yy
type (mp_real) epsm
type (mp_complex) zz, c1, c2, c3, c4, theta1, theta2, theta3, theta4
external theta1, theta2, theta3, theta4
integer lnm(n), lnm2(2*n)
character(1) chr1(lcx)
character(64) form4, nam(n), nam2(2*n), namx
character(lmx) line1
real (8) r1(n)
type (mp_real) r(n), r2(2*n), x(n)
save

!   Check to see if precision level is high enough.

if (ndp > mpipl) then
  write (6, 1) ndp
1 format ('Increase precision level in module MPFUNF to at least'/ &
    i8,' digits.')
  stop
endif

!   Check whether the palindromic option is appropriate.

if (ipal /= 0 .and. kd /= 2 * int (kd/2)) then
  write (6, 2)
2 format ('The palindromic technique may not be used when kd is odd.')
  stop
endif

write (6, 3) ipal, kd, kp, kq, mpi, n, ndp, ndpm, nep, nepm
3 format ('Poisson Phi_2 computation and analysis:'/ &
  'ipal =', i4,'; kd =',i4,'; kp =',i4,'; kq =',i4,'; mpi =',i4,'; n =',i4/ &
  'ndp =',i6,'; ndpm =',i6,'; nep =',i6,'; nepm =',i6)

pi = mppi (nwds)
eps = mpreald (10.d0, nwds) ** nep
epsm = mpreald (10.d0, nwdsm) ** nepm
xx = mpreal (dble (kp), nwds) / mpreal (dble (kd), nwds)
yy = mpreal (dble (kq), nwds) / mpreal (dble (kd), nwds)

!   Evaluate formula (43) of Section 5.

tm0 = second ()
qq = exp (- pi)
zz = 0.5d0 * pi * mpcmplx (yy, xx, nwds)
c1 = theta1 (zz, qq, eps, nwds)
c2 = theta2 (zz, qq, eps, nwds)
c3 = theta3 (zz, qq, eps, nwds)
c4 = theta4 (zz, qq, eps, nwds)
al = abs (c2 * c4 / (c1 * c3)) ** (mpi / 2)
tm1 = second ()
write (6, 4) tm1 - tm0
4 format ('Alpha CPU time =',f12.2/'Alpha =')
call mpwrite (6, ndp+20, ndp, al)

!   Construct input x vector for PSLQM3.

if (ipal == 0) then
  x(1) = mpreal (1.d0, nwds)

  do i = 2, n
    x(i) = al * x(i-1)
  enddo

  nam(1) = '1'
  lnm(1) = 1

  do i = 2, n
    if (i <= 10) then
      write (namx, '("al^",i1)') i - 1
      nam(i) = namx(1:4)
      lnm(i) = 4
    elseif (i <= 100) then
      write (namx, '("al^",i2)') i - 1
      nam(i) = namx(1:5)
      lnm(i) = 5
    elseif (i <= 1000) then
      write (namx, '("al^",i3)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    elseif (i <= 10000) then
      write (namx, '("al^",i4)') i - 1
      nam(i) = namx(1:7)
      lnm(i) = 7
    else
      stop
    endif   
  enddo

  do i = 1, n
    lnm2(i) = lnm(i)
    nam2(i) = nam(i)
  enddo
elseif (ipal == 1) then
  beta = al + 1.d0 / al
  x(1) = mpreal (1.d0, nwds)

  do i = 2, n
    x(i) = beta * x(i-1)
  enddo

  nam(1) = '1'
  lnm(1) = 1

  do i = 2, n
    if (i <= 10) then
      write (namx, '("beta^",i1)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    elseif (i <= 100) then
      write (namx, '("beta^",i2)') i - 1
      nam(i) = namx(1:7)
      lnm(i) = 7
    elseif (i <= 1000) then
      write (namx, '("beta^",i3)') i - 1
      nam(i) = namx(1:8)
      lnm(i) = 8
    elseif (i <= 10000) then
      write (namx, '("beta^",i4)') i - 1
      nam(i) = namx(1:9)
      lnm(i) = 9
    else
      stop
    endif   
  enddo

  nam2(1) = '1'
  lnm2(1) = 1

  do i = 2, 2 * n
    if (i <= 10) then
      write (namx, '("al^",i1)') i - 1
      nam2(i) = namx(1:4)
      lnm2(i) = 4
    elseif (i <= 100) then
      write (namx, '("al^",i2)') i - 1
      nam2(i) = namx(1:5)
      lnm2(i) = 5
    elseif (i <= 1000) then
      write (namx, '("al^",i3)') i - 1
      nam2(i) = namx(1:6)
      lnm2(i) = 6
    elseif (i <= 10000) then
      write (namx, '("al^",i4)') i - 1
      nam2(i) = namx(1:7)
      lnm2(i) = 7
    else
      stop
    endif   
  enddo
endif

!   Perform relation search.

tm2 = second ()
call pslqm3 (idb, n, nwds, nwdsm, rb, eps, epsm, x, iq, r)
tm3 = second ()
write (6, 5) tm3 - tm2, tm3 - tm0
5 format ('PSLQM3 CPU time =',f12.2/'Total CPU time =',f12.2)

!   Output relation in two formats.

if (iq == 1) then

!   Produce format used below for output.

  form4 = '(64a1'
  i1 = 5

  do i = 65, lcx, 64
    form4(i1+1:i1+9) =  '/64a1'
    i1 = i1 + 5
  enddo

  form4(i1+1:i1+9) = '," * ",a)'
  i1 = i1 + 9
  
  if (ipal == 1) then
    write (6, 6)
6   format (/'Recovered half relation: 0 =')
  
    do i = 1, n
      if (r(i) .ne. 0.d0) then
        call mpfform (r(i), lcx, 0, chr1)

        do j = 1, lcx
          if (chr1(j) /= ' ') goto 110
        enddo

110     continue

        j1 = j      
        if (chr1(j1) /= '-') then
          if (j1 > 1) then
            j1 = j1 - 1
            chr1(j1) = '+'
          else
            do j = 1, lcx
              chr1(j) = '*'
            enddo
          endif  
        endif    
        write (6, form4) (chr1(j), j = 1, lcx), nam(i)(1:lnm(i))
      endif
    enddo
  endif

!    If palindromic property is used, expand to full polynomial.

  if (ipal == 0) then
    n1 = n

    do i = 1, n
      r2(i) = r(i)
    enddo
  elseif (ipal == 1) then
    n1 = 2 * n - 1
    call doublep (nwds, 2*n - 1, r, r2)
  endif

  write (6, 7)
7 format (/'Recovered full relation: 0 =')
  l1 = 0

!   Output full polynomial in Mathematica form.

  do i = 1, n1
    if (r2(i) .ne. 0.d0) then
      call mpfform (r2(i), lcx, 0, chr1)

      do j = 1, lcx
        if (chr1(j) /= ' ') goto 120
      enddo

120   continue

      j1 = j      
      if (chr1(j1) /= '-') then
        if (j1 > 1) then
          j1 = j1 - 1
          chr1(j1) = '+'
        else
          do j = 1, lcx
            chr1(j) = '*'
          enddo
        endif  
      endif    
      write (6, form4) (chr1(j), j = 1, lcx), nam2(i)(1:lnm2(i))

      if (l1 + 100 > lmx) then
        write (6, '("Error: line1 too long",2i8)') l1 + k, lmx
        stop
      endif
      line1(l1+1:l1+1) = ' '
      l1 = l1 + 1
      k = lcx - j1
  
      do j = 1, k
        line1(l1+j:l1+j) = chr1(j+j1-1)
      enddo
  
      l1 = l1 + k
      line1(l1+1:l1+lnm2(i)+1) = '*' // nam2(i)(1:lnm2(i))
      l1 = l1 + lnm2(i) + 1
    endif
  enddo

  i1 = 1
  write (6, 8)
8 format ('Output polynomial in Mathematica notation:')

130 continue

  if (i1 + 64 > l1) goto 140
  i2 = index (line1(i1+65:l1), ' ')
  if (i2 == 0) goto 140
  i2 = i1 + 64 + i2
  write (6, '(a)') line1(i1:i2) // '\ '
  i1 = i2 + 1
  goto 130

140 continue

  write (6, '(a)') line1(i1:l1)
endif

stop
end

subroutine doublep (nwds, n, polyh, poly)
use mpmodule
implicit none
integer i, j, k, n, nwds
type (mp_real) polyh(0:n/2), poly(-n/2:n/2), x(-n/2-1:n/2+1), y(-n/2-1:n/2+1)

do i = -n/2-1, n/2+1
  x(i) = mpreal (0.d0, nwds)
enddo

x(0) = polyh(n/2)

do j = 0, n / 2 - 1
  y(-j-1) = x(-j)
  y(j+1) = x(j)
  
  do k = -j, j
    y(k) = x(k-1) + x(k+1)
  enddo
  
  y(0) = y(0) + polyh(n/2-j-1)

  do k = -j-1, j+1
    x(k) = y(k)
  enddo
enddo

do k = -n/2, n/2
  poly(k) = x(k)
enddo

return
end

function theta1 (zz, qq, eps, nwds)
use mpmodule
implicit none
integer k1, k2, n, nwds
real (8) ds
type (mp_real) eps, qq, r1
type (mp_complex) theta1, t0, t1, t2, t3, t4, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = sqrt (sqrt (qq))
ds = -1.d0

do k1 = 1, 10001, 2
  ds = - ds
  k2 = k1 ** 2
  t2 = ds * r1 ** k2 * sin (dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta1: loop end error'
stop

100 continue

theta1 = 2.d0 * t0
return
end

function theta2 (zz, qq, eps, nwds)
use mpmodule
implicit none
integer k1, k2, n, nwds
type (mp_real) eps, qq, r1
type (mp_complex) theta2, t0, t1, t2, t3, t4, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = sqrt (sqrt (qq))

do k1 = 1, 10001, 2
  k2 = k1 ** 2
  t2 = r1 ** k2 * cos (dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta2: loop end error'
stop

100 continue

theta2 = 2.d0 * t0
return
end

function theta3 (zz, qq, eps, nwds)
use mpmodule
implicit none
integer k1, k2, n, nwds
type (mp_real) eps, qq, r1
type (mp_complex) theta3, t0, t1, t2, t3, t4, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = qq

do k1 = 1, 10000
  k2 = k1 ** 2
  t2 = r1 ** k2 * cos (2.d0 * dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta3: loop end error'
stop

100 continue

theta3 = cmplx (1.d0, 0.d0, mprknd) + 2.d0 * t0
return
end

function theta4 (zz, qq, eps, nwds)
use mpmodule
implicit none
integer k1, k2, n, ndp, nwds
real (8) ds
type (mp_real) eps, qq, r1
type (mp_complex) theta4, t0, t1, t2, t3, t4, zz

t0 = mpcmplx (cmplx (0.d0, 0.d0, mprknd), nwds)
r1 = qq
ds = 1.d0

do k1 = 1, 10000
  ds = - ds
  k2 = k1 ** 2
  t2 = ds * r1 ** k2 * cos (2.d0 * dble (k1) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta4: loop end error'
stop

100 continue

theta4 = cmplx (1.d0, 0.d0, mprknd) + 2.d0 * t0
return
end

!------------------------------

!   The following code performs the three-level, multi-pair PSLQ algorithm.
!   David H. Bailey    19 Apr 2017

subroutine pslqm3 (idb, n, nwds, nwdsm, rb, eps, epsm, x, iq, r)

!   Arguments are as follows:
!     Name  Type    Description
!     idb   int*4   Debug flag (0-4); increasing idb produces more output.
!     n     int*4   Length of input vector x and output relation vector r.
!     nwds  int*4   Full precision level in words. This must be sufficient
!                     to recover the relation, plus nwdsm words.
!     nwdsm int*4   Medium pecision level in words. This must be sufficient
!                     to handle the typically large dynamic range at start.
!     rb    real*8  Log10 of max size (Euclidean norm) of acceptable relation.
!                   NOTE: Run is aborted when this is exceeded.
!     eps   mp_real Tolerance for full precision relation detection. This is
!                     typically set to 2^(100-mpnbt*nwds) or so.
!     epsm  mp_real Tolerance for intermediate precision detection. This is
!                     typically set to 2^(70-mpnbt*nwdsm) or so.
!     x     mp_real Input mp_real vector.
!     iq    int*4   Output flag: 0 (unsuccessful) or 1 (successful).
!     r     mp_real Output integer relation vector, if successful.

!   The following parameters are set in this routine:
!     ipi   int*4   Iteration print interval when idb >= 2; default = 500.
!     ipm   int*4   Iteration check interval for MPM iterations; default = 10.
!     itm   int*4   Maximum iteration count; default = 10^7.
!                   NOTE: Run is aborted when this is exceeded. If itm >= 10^8,
!                     change all "i8" in formats to "i9" or as high as needed.
!     nrs   int*4   0: Restart file is not read or written; default = 0.
!                   1: Start new run, writing new restart file.
!                   2: Read restart file and then run, writing restart file.
!                   NOTE: When nrs = 2, input restart file is overwritten,
!                   so make a backup copy before running this program.
!     nsq   int*4   Size of tables used in iterdp and itermpw; default = 8.
!     deps  real*8  Tolerance for dynamic range check; default = 1d-10.

use mpmodule
implicit none
integer i, idb, imq, ipi, ipm, iq, it, itm, its, izd, izm, izmm, j, j1, n, &
  n1, n2, n3, n4, nrs, nsq, nwds, nwdsm
parameter (ipi = 500, ipm = 10, itm = 10000000, nrs = 0, nsq = 8)
real (8) d1, d2, d3, d4, deps, dplog10, dplog10m, rb, second, &
  tm0, tm1, times(6)
parameter (deps = 1d-10)
real (8) da(n,n), db(n,n), dh(n,n), dq(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsyq(n,nsq), dy(n), dsy(n)
type (mp_real) b(n,n), h(n,n), r(n), x(n), y(n), eps, t1, t2
type (mp_real) bound, dynrange, dynrangem, epsm, wa(n,n), wb(n,n), wh(n,n), &
  wsyq(n,nsq), wy(n), wn, w1, w2
external bound, dplog10, dplog10m, dynrange, dynrangem, second

!   Initialize.

if (idb >= 2) write (6, 1) n
1 format ('PSLQM3 integer relation detection: n =',i5)
iq = 0
it = 0
its = 0
imq = 0
wn = mpreald (0.d0, nwdsm)

if (idb >= 2) write (6, 2) it
2 format ('Iteration',i8,3x,'MP initialization')

!   If nrs >= 1, open restart file; if nrs = 2, read restart file. 

if (nrs >= 1) open (12, file = 'pslqm3.rst', form = 'unformatted')
if (nrs <= 1) then
  do i = 1, 6
    times(i) = 0.d0
  enddo

  tm0 = second ()
  call initmp (idb, n, nwds, b, h, x, y)
  tm1 = second ()
  times(1) = tm1 - tm0
else
  rewind 12
  read (12) it
  read (12) times
  read (12) y
  read (12) b
  read (12) h
  rewind 12
endif

100 continue

!   Check if dynamic range of y vector is too great for MPM iterations;
!   if so, then abort. In the main calling program, the medium precision
!   level nwdsm must be set to at least high enough to handle the initial
!   dynamic range for a given problem.

w1 = dynrange (n, nwdsm, y)
if (idb >= 2) then
  call decmdm (w1, d1, n1)
  write (6, 3) it, d1, n1
3 format ('Iteration',i8,3x,'Min/max ratio in y =',f11.6,'e',i6)
endif
if (w1 < epsm) then
  if (idb >= 1) write (6, 4)
4 format ('Run aborted.')
  goto 160
endif

!   Initialize MPM arrays from MP arrays.

if (idb >= 3) write (6, 5) it
5 format ('Iteration',i8,3x,'MPM initialization')
tm0 = second ()
call initmpm (idb, n, nsq, nwdsm, wa, wb, wh, wy, wsyq, h, y)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

110 continue

!   Check if dynamic range of wy vector is too great for DP iterations
!   (which is often the case at the start of the run). If so, then
!   perform MPM iterations instead of DP iterations.

w1 = dynrangem (n, nwdsm, wy)
if (idb >= 3) then
  call decmdm (w1, d1, n1)
  write (6, 6) it, d1, n1
6 format ('Iteration',i8,3x,'Min/max ratio in wy =',f11.6,'e',i6)
endif
if (w1 < mpreald (deps, nwdsm)) then
  goto 130
endif

!   DP iterations.
!   Initialize DP arrays from MPM arrays.

if (idb >= 3) write (6, 7) it
7 format ('Iteration',i8,3x,'Start DP iterations')
call initdp (idb, n, nsq, nwdsm, da, db, dh, dy, dsyq, wh, wy)

!   Save DP arrays and iteration count.

call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
its = it

!   Perform an LQ decomposition on DH.

call lqdp (n, n - 1, dh)

120 continue

!   Perform one DP iteration.

it = it + 1
if (idb >= 4 .or. idb >= 2 .and. mod (it, ipi) == 0) write (6, 8) it
8 format ('Iteration',i8)
tm0 = second ()
call iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)
tm1 = second ()
times(3) = times(3) + (tm1 - tm0)

!   Test conditions on iterdp output flag izd:
!   0: Iteration was uneventful; periodically save arrays and continue.
!   1: Relation found or DP precision exhausted; perform MPM update.
!   2: Very large value appeared in DA or DB; abort DP iter and do MPM update.

if (izd == 0) then
  if (mod (it - its, ipm) == 0) then
    call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
    its = it
  endif
  goto 120
else

!   Check if DP iteration was aborted above; if so, then revert to previous data.

  if (izd == 2) then
    it = its
    call savedp (n, dsa, dsb, dsh, dsy, da, db, dh, dy)
  endif

!   Update the MPM arrays from the DP arrays.

  if (idb >= 3) write (6, 9) it
9 format ('Iteration',i8,3x,'MPM update')
  tm0 = second ()
  call updtmpm (idb, it, n, nwdsm, epsm, da, db, wa, wb, wh, wy, izmm)
  tm1 = second ()
  times(4) = times(4) + (tm1 - tm0)

!   Test conditions on updtmpm output flag izmm:
!   0: MPM update was uneventful; if izd = 2, goto MPM iter; otherwise continue.
!   1: MPM found relation or MPM precision is exhausted; perform MP update.
!   2: Very large value arose in wa or wb arrays; abort run.

  if (izmm == 0) then

!   If DP iteration was aborted above, then go perform MPM iterations,
!   Otherwise restart DP iterations and continue.

    if (izd == 2) then
      goto 130
    else
      goto 110
    endif
  elseif (izmm == 1) then

!   MPM update above either found relation, or precision was exhausted;
!   either way, update the MP arrays from the MPM arrays.

    if (idb >= 2) write (6, 10) it
10  format ('Iteration',i8,3x,'MP update')
    tm0 = second ()
    call updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)
    tm1 = second ()
    times(5) = times(5) + (tm1 - tm0)

!   If checkpoint-restart is enabled, then save the latest MP arrays.

    if (nrs >= 1) then
      rewind 12
      write (12) it
      write (12) times
      write (12) y
      write (12) b
      write (12) h
      rewind 12
    endif

!   Compute norm bound, using MPM precision.

    call lqmpm (n, n - 1, nwdsm, wh)
    w1 = bound (n, nwdsm, wh)
    call decmdm (w1, d3, n3)
    wn = max (wn, w1)
    call decmdm (wn, d4, n4)
    if (idb >= 2) then
      write (6, 11) it, d3, n3, d4, n4
11    format ('Iteration',i8,3x,'Norm bound =',f11.6,'e',i5,3x, &
        'Max. bound =',f11.6,'e',i5)
    endif

!   Check if iteration limit or norm bound limit is exceeded; if so, quit.

    if (it > itm) then
      if (idb >= 1) write (6, 12) itm
12    format ('Iteration limit exceeded',i8)
      goto 160
    endif
    if (dplog10m (wn) > rb) then
      if (idb >= 1) write (6, 13) rb
13    format ('Norm bound limit exceeded.',1p,d15.6)
      goto 160
    endif

!   Test conditions on updtmp output flag izm:
!   0: MP update was uneventful; goto tag 100 above.
!   1: A small value was found in MP update; go output relation.
!   2: Precision is exhausted; quit.

    if (izm == 0) then
      goto 100
    elseif (izm == 1) then
      goto  150
    elseif (izm == 2) then
      goto 160
    endif
  elseif (izmm == 2) then

!   A very large value was found in wa or wb -- abort.

    goto 160
  endif
endif

130 continue

!   MPM iterations.

if (idb >= 2) write (6, 14) it, nint (nwdsm * mpdpw)
14 format ('Iteration',i8,3x,'Start MPM iterations: precision =',i6,' digits')

!   Perform LQ decomposition using MPM precision.

tm0 = second ()
call lqmpm (n, n - 1, nwdsm, wh)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Perform MPM iterations.

140 continue

!   Perform one MPM iteration.

it = it + 1
if (idb >= 2) write (6, 15) it
15 format ('Iteration',i8)
tm0 = second ()
call itermpm (idb, it, n, nsq, nwdsm, epsm, wa, wb, wh, wsyq, wy, imq, izmm)
tm1 = second ()
times(6) = times(6) + (tm1 - tm0)

!   Test conditions on itermpm output flag izmm:
!   0: Iteration was uneventful; periodically check for DP; continue.
!   1: A small value was found or MPM precision exhausted; do MP update.
!   2: A very large value appeared in WA or WB; abort run.

if (izmm == 0) then

!   Periodically check to see if DP iterations can be resumed (but perform
!   at least IPM iterations in MPM).

  if (mod (it - its, ipm) == 0) then
    w1 = dynrangem (n, nwdsm, wy)
    if (idb >= 2) then
      call decmdm (w1, d1, n1)
      write (6, 16) it, d1, n1
16    format ('Iteration',i8,3x,'Min/max ratio in wy =',f11.6,'e',i6)
    endif
    if (w1 > mpreald (deps, nwdsm)) then 
      if (idb >= 2) write (6, 17) it
17    format ('Iteration',i8,3x,'Return to DP iterations')
      goto 110
    endif
  endif
  goto 140
elseif (izmm == 1) then

!   Update the MP arrays from the MPM arrays.

  if (idb >= 2) write (6, 18) it
18 format ('Iteration',i8,3x,'MP update')
  tm0 = second ()
  call updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)
  tm1 = second ()
  times(5) = times(5) + (tm1 - tm0)

!   Test conditions on updtmp output flag izm:
!   0: MP update was uneventful; goto tag 100 above.
!   1: A small value was found in MP update; go output relation.
!   2: Precision is exhausted; quit.

  if (izm == 0) then
    goto 100
  elseif (izm == 1) then
    goto 150
  elseif (izm == 2) then
    goto 160
  endif
elseif (izmm == 2) then

!   During itermpm, a very large entry was produced in wa or wb -- abort.

  goto 160
endif

150 continue

!   A relation has been detected.  Output the final norm bound and other info.

t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
enddo

!   The norm calculation here is performed in medium precision.

w1 = mpreald (0.d0, nwdsm)

do i = 1, n
  w1 = w1 + mpreal (r(i), nwdsm) ** 2
enddo

w1 = sqrt (w1)
call decmdm (w1, d3, n3)

!   Output the final norm bound and other information.

if (idb >= 1) then
  call lqmpm (n, n - 1, nwdsm, wh)
  w2 = bound (n, nwdsm, wh)
  wn = max (wn, w2)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmdm (wn, d4, n4)
  write (6, 19) it, d1, n1, d2, n2, d4, n4
19 format ('Iteration',i8,3x,'Relation detected'/ &
  'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6/'Max. bound =',f11.6,'e',i6)
  write (6, 20) j1, d3, n3, d1, n1
20 format ('Index of relation =',i4,3x,'Norm =',f11.5,'e',i5,3x, &
  'Residual =',f11.6,'e',i6)
endif

!   If run was successful, set iq = 1.

if (dplog10m (w1) <= rb) then
  iq = 1
else
  if (idb >= 2) write (6, 21)
21 format ('Relation is too large.')
endif

160 continue

!   Output CPU run times and return.

if (idb >= 2) write (6, 22) times
22 format ('CPU run times:'/(6f12.2))

return
end

!------------------------------

!   First-level subroutines.

function dynrange (n, nwdsm, y)

!   This returns the dynamic range of y, i.e., ratio of min|y_k| to max|y_k|,
!   using MPM precision.

use mpmodule
implicit none
integer i, n, n1, nwdsm
real (8) d1
type (mp_real) dynrange, t1, t2, t3
type (mp_real) y(n)

t1 = mpreald (1.d300, nwdsm)
t2 = mpreald (0.d0, nwdsm)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (mpreal (y(i), nwdsm))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

dynrange = t1 / t2
return
end

function dynrangem (n, nwdsm, wy)

!   This returns the dynamic range of wy, i.e., ratio of min|wy_k| to max|wy_k|,
!   using MPM precision.

use mpmodule
implicit none
integer i, n, n1, nwdsm
real (8) d1
type (mp_real) dynrangem, wy(n), t1, t2, t3

t1 = mpreald (1.d300, nwdsm)
t2 = mpreald (0.d0, nwdsm)

!   Find the min and max absolute value in the wy vector.

do i = 1, n
  t1 = min (t1, abs (wy(i)))
  t2 = max (t2, abs (wy(i)))
enddo

dynrangem = t1 / t2
return
end

subroutine initdp (idb, n, nsq, nwdsm, da, db, dh, dy, dsyq, wh, wy)

!   This initializes the DP arrays from the MPM arrays.
!   This is performed in medium precision.
!   Input:  idb, n, nsq, nwdsm, wh, wy.
!   Output: da, db, dh, dy, dsyq.

use mpmodule
implicit none
integer i, idb, j, n, nsq, nwdsm
real (8) da(n,n), db(n,n), dh(n,n), dy(n), dsyq(n,nsq)
type (mp_real) wh(n,n), wy(n), t1, t2, t3, t4

t2 = mpreald (0.d0, nwdsm)

!   Find the max absolute value in the wy vector.

do i = 1, n
  t2 = max (t2, abs (wy(i)))
enddo

!   Set dy to be the scaled wy vector.

t1 = 1.d0 / t2

do i = 1, n
  dy(i) = t1 * wy(i)
enddo

!   Find the maximum absolute value of the wh matrix diagonals.

t2 = mpreald (0.d0, nwdsm)

do j = 1, n - 1
  t2 = max (t2, abs (wh(j,j)))
enddo

!   Set dh to be the scaled wh matrix.

t1 = 1.d0 / t2

do j = 1, n - 1
  do i = 1, n
    dh(i,j) = t1 * wh(i,j)
  enddo
enddo

!   Set da and db to the identity.

do j = 1, n
  do i = 1, n
    da(i,j) = 0.d0
    db(i,j) = 0.d0
  enddo

  da(j,j) = 1.d0
  db(j,j) = 1.d0
enddo

!   Zero the dsyq array.

do j = 1, nsq
  do i = 1, n
    dsyq(i,j) = 0.d0
  enddo
enddo

if (idb >= 4) then
  write (6, 2)
2 format ('initdp: Scaled dy vector:')
  call matoutdp (1, n, dy)
  write (6, 3)
3 format ('initdp: Scaled dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine initmp (idb, n, nwds, b, h, x, y)

!   Initializes MP arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nwds.
!   Output: b, h, x, y.

use mpmodule
implicit none
integer i, idb, j, n, nwds
integer ix(n)
real (8) dx(n)
type (mp_real) b(n,n), h(n,n), s(n), x(n), y(n), t1

if (idb >= 4) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = mpreald (0.d0, nwds)
  enddo

  b(j,j) = mpreald (1.d0, nwds)
enddo

t1 = mpreald (0.d0, nwds)

!   Compute the x vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.d0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = mpreald (0.d0, nwds)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo

if (idb >= 4) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine initmpm (idb, n, nsq, nwdsm, wa, wb, wh, wy, wsyq, h, y)

!   This initializes the MPM arrays from the MP arrays.
!   This is performed in medium precision.
!   Input: idb, n, nsq, nwdsm, h, y.
!   Output: wa, wb, wh, wy, wsyq.

use mpmodule
implicit none
integer i, idb, j, n, n1, nsq, nwdsm
real (8) d1
integer ix(n)
real (8) dx(n)
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wy(n), wsyq(n,nsq), t1, t2, t3
type (mp_real) h(n,n), y(n)

t2 = mpreald (0.d0, nwdsm)

!   Find the max absolute value in the y vector.

do i = 1, n
  t3 = abs (mpreal (y(i), nwdsm))
  t2 = max (t2, t3)
enddo

!   Set wy to the scaled y vector.

t1 = 1.d0 / t2

do i = 1, n
  wy(i) = t1 * mpreal (y(i), nwdsm)
enddo

!   Set wh to h.

do j = 1, n - 1
  do i = 1, n
    wh(i,j) = mpreal (h(i,j), nwdsm)
  enddo
enddo

!   Set wa and wb to the identity.

do j = 1, n
  do i = 1, n
    wa(i,j) = mpreald (0.d0, nwdsm)
    wb(i,j) = mpreald (0.d0, nwdsm)
  enddo

  wa(j,j) = mpreald (1.d0, nwdsm)
  wb(j,j) = mpreald (1.d0, nwdsm)
enddo

!   Zero the wsyq array.

do j = 1, nsq
  do i = 1, n
    wsyq(i,j) = mpreald (0.d0, nwdsm)
  enddo
enddo

if (idb >= 4) then
  write (6, 3)
3 format ('initmpm: wy vector:')
  call matoutmdm (1, n, ix, dx, wy)
  write (6, 4)
4 format ('initmpm: Factored wh matrix:')
  call matoutmdm (n, n - 1, ix, dx, wh)
endif

100 continue

return
end

subroutine iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)


!   This performs one iteration of the PSLQ algorithm using DP arithmetic.
!   Input: idb, it, n, nsq, da, db, dh, dsyq, dy.
!   Output: da, db, dh, dsyq, dy, imq, izd.

!   NOTE: Parameter tmx2 = 2^52, not 2^53, so as to ensure that values > 2^53
!   never arise, even as intermediate values, in the update loop below.

implicit none
integer i, idb, ii, ij, im, im1, imq, it, izd, j, j1, j2, k, mpr, mq, n, nsq
real (8) deps, gam, t1, t2, t3, t4, tmx1, tmx2
integer ip(n), ir(n), is(n)
real (8) da(n,n), db(n,n), dh(n,n), dq(n), dsyq(n,nsq), &
  dt(n,n), dy(n)
parameter (tmx1 = 1.d13, tmx2 = 2.d0**52, deps = 1.d-14)

izd = 0
mpr = nint (0.4d0 * n)
gam = sqrt (4.d0 / 3.d0)

!   Compute dq vector = {gam^i * |dh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  dq(i) = gam ** i * abs (dh(i,i))
enddo

call qsortdp (n - 1, dq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest dq(i).

do i = 1, n
  is(i) = 0
enddo

!  If the imq flag is set, perform this iteration with mq = 1.

if (imq == 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii == 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of dy, and rows of da, db and dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = dy(im)
  dy(im) = dy(im1)
  dy(im1) = t1

  do i = 1, n
    t1 = da(im,i)
    da(im,i) = da(im1,i)
    da(im1,i) = t1
    t1 = db(im,i)
    db(im,i) = db(im1,i)
    db(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = dh(im,i)
    dh(im,i) = dh(im1,i)
    dh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
    t1 = dh(im,im)
    t2 = dh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = dh(i,im)
      t4 = dh(i,im1)
      dh(i,im) = t1 * t3 + t2 * t4
      dh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on dh, using the diagonal scheme.  Multipliers are
!   saved in the dt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      dh(ij,j) = dh(ij,j) - dt(ij,k) * dh(k,j)
    enddo

    dt(ij,j) = anint (dh(ij,j) / dh(j,j))
    dh(ij,j) = dh(ij,j) - dt(ij,j) * dh(j,j)
  enddo
enddo

!   Update dy, using the dt array.  Find min absolute value of dy.

t1 = abs (dy(n))

do j = 1, n - 1
  do i = j + 1, n
    dy(j) = dy(j) + dt(i,j) * dy(i)
  enddo

  t1 = min (t1, abs (dy(j)))
enddo

!   Update da and db, using the dt array.  Find the max absolute value of
!   da and db entries as they are calculated (not merely at the end).

t2 = 0.d0

do k = 1, n
  dq(k) = 0.d0

  do j = 1, n - 1
    do i = j + 1, n
      da(i,k) = da(i,k) - dt(i,j) * da(j,k)
      db(j,k) = db(j,k) + dt(i,j) * db(i,k)
      dq(k) = max (dq(k), abs (da(i,k)), abs (db(j,k)))
    enddo
  enddo
enddo

do k = 1, n
  t2 = max (t2, dq(k))
enddo

if (t1 <= deps) then
  if (idb >= 3) write (6, 1) it, t1
1 format ('Iteration',i8,3x,'iterdp: Small value in dy =',1pd15.6)
  izd = 1
endif

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 3) write (6, 2) it, t2
2 format ('Iteration',i8,3x,'iterdp: Large value in da or db =',1pd15.6)
  izd = 1
elseif (t2 > tmx2) then
  if (idb >= 2) write (6, 3) it, t2
3 format ('Iteration',i8,3x,'iterdp: Very large value in da or db =',1pd15.6)
  izd = 2
  return
endif

!   Compare the dy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = 0.d0

  do i = 1, n
    t1 = max (t1, abs (dy(i) - dsyq(i,j)))
  enddo

  if (t1 <= deps) then
    if (idb >= 2) write (6, 4) it, j
4   format ('Iteration',i8,3x,'iterdp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector dy in the table dsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  dsyq(i,k) = dy(i)
enddo

if (idb >= 4) then
  write (6, 5)
5 format ('iterdp: Updated dy:')
  call matoutdp (1, n, dy)
  write (6, 6)
6 format ('iterdp: Updated da matrix:')
  call matoutdp (n, n, da)
  write (6, 7)
7 format ('iterdp: Updated db matrix:')
  call matoutdp (n, n, db)
  write (6, 8)
8 format ('iterdp: Updated dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine itermpm (idb, it, n, nsq, nwdsm, epsm, wa, wb, wh, wsyq, wy, imq, izmm)

!   This performs one iteration of the PSLQM algorithm using MPM arithmetic.
!   This is performed in medium precision.
!   Input: idb, it, n, nsq, nwdsm, epsm, imq.
!   Output: wa, wb, wh, wsyq, wy, imq, izmm.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izmm, j, j1, j2, k, mpr, mq, n, n1, &
  nsq, ntl, nwdsm
parameter (ntl = 72)
real (8) d1
type (mp_real) gam, t1, t2, t3, t4, epsm, tmx1, tmx2
integer ip(n), ir(n), is(n)
real (8) dx(n)
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wq(n), wsyq(n,nsq), wt(n,n), wy(n)

tmx1 = 1.d0 / epsm
tmx2 = mpreald (2.d0, nwdsm) ** (nwdsm * mpnbt)
izmm = 0
mpr = nint (0.4d0 * n)
gam = sqrt (mpreald (4.d0, nwdsm) / mpreald (3.d0, nwdsm))

!   Compute wq vector = {gam^i * |wh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  wq(i) = gam ** i * abs (wh(i,i))
enddo

call qsortmpm (n - 1, wq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest wq(i).

do i = 1, n
  is(i) = 0
enddo

if (imq == 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii == 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of wy, and rows of wa, wb and wh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = wy(im)
  wy(im) = wy(im1)
  wy(im1) = t1

  do i = 1, n
    t1 = wa(im,i)
    wa(im,i) = wa(im1,i)
    wa(im1,i) = t1
    t1 = wb(im,i)
    wb(im,i) = wb(im1,i)
    wb(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = wh(im,i)
    wh(im,i) = wh(im1,i)
    wh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in wh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
    t1 = wh(im,im)
    t2 = wh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = wh(i,im)
      t4 = wh(i,im1)
      wh(i,im) = t1 * t3 + t2 * t4
      wh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on wh, using the diagonal scheme.  Multipliers are
!   saved in the wt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      wh(ij,j) = wh(ij,j) - wt(ij,k) * wh(k,j)
    enddo

    wt(ij,j) = anint (wh(ij,j) / wh(j,j))
    wh(ij,j) = wh(ij,j) - wt(ij,j) * wh(j,j)
  enddo
enddo

!   Update wy, using the wt array.  Find min absolute value of wy.

t1 = abs (wy(n))

do j = 1, n - 1
  do i = j + 1, n
    wy(j) = wy(j) + wt(i,j) * wy(i)
  enddo

  t1 = min (t1, abs (wy(j)))
enddo

!   Update wa and wb, using the wt array.  Find the max absolute value of
!   wa and wb entries as they are calculated (not merely at the end).

t2 = mpreald (0.d0, nwdsm)

do k = 1, n
  wq(k) = mpreald (0.d0, nwdsm)

  do j = 1, n - 1
    do i = j + 1, n
      wa(i,k) = wa(i,k) - wt(i,j) * wa(j,k)
      wb(j,k) = wb(j,k) + wt(i,j) * wb(i,k)
      wq(k) = max (wq(k), abs (wa(i,k)), abs (wb(j,k)))
    enddo
  enddo
enddo

do k = 1, n
  t2 = max (t2, wq(k))
enddo

if (t1 <= epsm) then
  if (idb >= 2) then
    call decmdm (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i8,3x,'itermpm: Small value in wy =',f11.6,'e',i6)
  endif
  izmm = 1
endif

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 2) then
    call decmdm (t2, d1, n1)
    write (6, 2) it, d1, n1
2   format ('Iteration',i8,3x,'itermpm: Large value in wa or wb =', &
      f11.6,'e',i6)
  endif
  izmm = 1
elseif (t2 > tmx2) then
  if (idb >= 1) then
    call decmdm (t2, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i8,3x,'itermpm: Very large value in wa or wb =', &
      f11.6,'e',i6/'Run aborted.')
  endif
  izmm = 2
  goto 200
endif

!   Compare the wy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = mpreald (0.d0, nwdsm)

  do i = 1, n
    t1 = max (t1, abs (wy(i) - wsyq(i,j)))
  enddo

  if (t1 <= epsm) then
    if (idb >= 2) write (6, 4) it, j
4   format ('Iteration',i8,3x,'itermpm: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector wy in the table wsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  wsyq(i,k) = wy(i)
enddo

if (idb >= 4) then
  write (6, 5)
5 format ('itermpm: Updated wy:')
  call matoutmdm (1, n, ip, dx, wy)
  write (6, 6)
6 format ('itermpm: Updated wa matrix:')
  call matoutmdm (n, n, ip, dx, wa)
  write (6, 7)
7 format ('itermpm: Updated wb matrix:')
  call matoutmdm (n, n, ip, dx, wb)
  write (6, 8)
8 format ('itermpm: Updated wh matrix:')
  call matoutmdm (n, n - 1, ip, dx, wh)
endif

200 continue

return
end

subroutine lqdp (n, m, dh)

!   This performs an LQ decomposition on the DP matrix dh.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   Input: n, m, dh.
!   Output: dh.

implicit none
integer i, j, l, lup, m, ml, n
real (8) dh(n,m), nrmxl, one, t, zero

zero = 0.d0
one = 1.d0
lup = min (m,n)

!   Perform the householder reduction of dh.

do l = 1, lup
  if (l == m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + dh(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl == zero) go to 270
  if (dh(l,l) .ne. zero) nrmxl = sign (nrmxl, dh(l,l))
  t = one / nrmxl

  do i = 0, ml
    dh(l,l+i) = t * dh(l,l+i)
  enddo

  dh(l,l) = one + dh(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + dh(l,l+i) * dh(j,l+i)
    enddo

    t = - t / dh(l,l)

    do i = 0, ml
      dh(j,l+i) = dh(j,l+i) + t * dh(l,l+i)
    enddo
  enddo

!   Save the transformation.

  dh(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero dh above the diagonal.

do j = 1, m
  do i = 1, j - 1
    dh(i,j) = 0.d0
  enddo
enddo

return
end

subroutine lqmpm (n, m, nwdsm, h)

!   This performs an LQ decomposition on the MPM matrix h.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   This is performed in medium precision.
!   Input: n, m, nwdsm, h.
!   Output: h.

use mpmodule
implicit none
integer i, j, l, lup, m, ml, n, nwdsm
type (mp_real) h(n,m), nrmxl, one, t, zero

zero = mpreald (0.d0, nwdsm)
one = mpreald (1.d0, nwdsm)
lup = min (m,n)

!   Perform the householder reduction of h.

do l = 1, lup
  if (l == m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + h(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl == zero) go to 270
  if (h(l,l) .ne. zero) nrmxl = sign (nrmxl, h(l,l))
  t = one / nrmxl

  do i = 0, ml
    h(l,l+i) = t * h(l,l+i)
  enddo

  h(l,l) = one + h(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + h(l,l+i) * h(j,l+i)
    enddo

    t = - t / h(l,l)

    do i = 0, ml
      h(j,l+i) = h(j,l+i) + t * h(l,l+i)
    enddo
  enddo

!   Save the transformation.

  h(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero h above the diagonal.

do j = 1, m
  do i = 1, j - 1
    h(i,j) = mpreald (0.d0, nwdsm)
  enddo
enddo

return
end

subroutine savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)

!   This saves the arrays dy, da, db, dh in case dp iterations must be aborted.
!   A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
!   exchanged, serves to restore these arrays.
!   Input: n, da, db, dh, dy.
!   Output: dsa, dsb, dsh, dsy.

implicit none
integer i, j, n
real (8) da(n,n), db(n,n), dh(n,n), dy(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsy(n)

do i = 1, n
  dsy(i) = dy(i)
enddo

do j = 1, n
  do i = 1, n
    dsa(i,j) = da(i,j)
    dsb(i,j) = db(i,j)
  enddo
enddo

do j = 1, n - 1
  do i = 1, n
    dsh(i,j) = dh(i,j)
  enddo
enddo

return
end

subroutine updtmp (idb, it, n, nwds, wa, wb, eps, b, h, y, izm)

!   This update the MP arrays from the MPM arrays.
!   This is performed in full precision.
!   Input: idb, it, n, nwds, wa, wb, eps, b, h, y.
!   Output: b, h, y, izm.

use mpmodule
implicit none 
integer i, i1, idb, it, izm, n, n1, n2, ntl, nwds
parameter (ntl = 72)
real (8) d1, d2
type (mp_real) eps, t1, t2, teps
integer ix(n)
real (8) dx(n)
type (mp_real) wa(n,n), wb(n,n)
type (mp_real) b(n,n), h(n,n), y(n)

teps = mpreald (2.d0, nwds) ** ntl * eps
izm = 0

!   Update y with wb.

call mxmm (n, 1, nwds, wb, y)
i1 = 0
t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

do i = 1, n
  if (abs (y(i)) < t1) then
    i1 = i
    t1 = abs (y(i))
  endif
  t2 = max (t2, abs (y(i)))
enddo

if (idb >= 2) then
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  write (6, 1) it, d1, n1, d2, n2
1  format ('Iteration',i8,3x,'updtmp: Min, max of y =',f11.6,'e',i6, &
    f11.6,'e',i6)
endif

!   Update b with wb.

call mxmm (n, n, nwds, wb, b)

!   Update h with wa.

call mxmm (n, n - 1, nwds, wa, h)

!   Find the largest entry of b in the same row as the smallest y.

t2 = mpreald (0.d0, nwds)

do i = 1, n
  t2 = max (t2, abs (b(i1,i)))
enddo

if (t1 <= t2 * teps) then
  if (idb >= 2) then
    call decmd (t1, d1, n1) 
    write (6, 2) it, d1, n1
2   format ('Iteration',i8,3x,'updtmp: Small value in y =',f11.6,'e',i6)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 3) it
3   format ('Iteration',i8,3x,'updtmp: Precision exhausted.')
    izm = 2
  endif
endif

if (idb >= 4) then
  write (6, 4)
4 format ('updtmp: Updated y:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 5)
5 format ('updtmp: Updated b matrix:')
  call matoutmd (n, n, ix, dx, b)
  write (6, 6)
6 format ('updtmp: Updated h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine updtmpm (idb, it, n, nwdsm, epsm, da, db, wa, wb, wh, wy, izmm)
 
!   This updates the MPM arrays from the DP arrays.
!   This is performed in medium precision.
!   Input: idb, it, n, nwdsm, epsm, da, db, wa, wb, wh, wy.
!   Output: wa, wb, wh, wy, izmm.

use mpmodule
implicit none 
integer i, idb, it, izmm, j, n, n1, n2, ntl, nwdsm
parameter (ntl = 72)
real (8) d1, d2
type (mp_real) t1, t2, epsm, tmx1, tmx2
integer ix(n)
real (8) da(n,n), db(n,n), dx(n)
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wy(n)

tmx1 = 1.d0 / epsm
tmx2 = mpreald (2.d0, nwdsm) ** (nwdsm * mpnbt)
izmm = 0
t1 = mpreald (1.d300, nwdsm)
t2 = mpreald (0.d0, nwdsm)

!   Update wy with db.

call mxmdm (n, 1, nwdsm, db, wy)

do i = 1, n
  t1 = min (t1, abs (wy(i)))
  t2 = max (t2, abs (wy(i)))
enddo

if (idb >= 3) then
  call decmdm (t1, d1, n1)
  call decmdm (t2, d2, n2)
  write (6, 1) it, d1, n1, d2, n2
1 format ('Iteration',i8,3x,'updtmpm: Min, max of wy =',f11.6,'e',i6, &
     f11.6,'e',i6)
endif
if (t1 <= epsm) then
  if (idb >= 2) then
    call decmdm (t1, d1, n1)
    write (6, 2) it, d1, n1
2   format ('Iteration',i8,3x,'updtmpm: Small value in wy =',f11.6,'e',i6)
  endif
  izmm = 1
endif

!   Update wa with da.

call mxmdm (n, n, nwdsm, da, wa)

!   Update wb with db.

call mxmdm (n, n, nwdsm, db, wb)
t2 = mpreald (0.d0, nwdsm)

do j = 1, n
  do i = 1, n
    t2 = max (t2, abs (wa(i,j)), abs (wb(i,j)))
  enddo
enddo

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 2) then
    call decmdm (t2, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i8,3x,'updtmpm: Large value in wa or wb =', &
      f11.6,'e',i6)
  endif
  izmm = 1
elseif (t2 > tmx2) then
  if (idb >= 1) then
    call decmdm (t2, d1, n1)
    write (6, 4) it, d1, n1
4   format ('updtmpm: Very large value in wa or wb =',f11.6,'e',i6/ &
      'Run aborted.')
  endif
  izmm = 2
  goto 100
endif

!   Update wh with da.

call mxmdm (n, n - 1, nwdsm, da, wh)

if (idb >= 4) then
  write (6, 5)
5 format ('updtmpm: Updated wy:')
  call matoutmdm (1, n, ix, dx, wy)
  write (6, 6)
6 format ('updtmpm: Updated wa matrix:')
  call matoutmdm (n, n, ix, dx, wa)
  write (6, 7)
7 format ('updtmpm: Updated wb matrix:')
  call matoutmdm (n, n, ix, dx, wb)
  write (6, 8)
8 format ('updtmpm: Updated wh matrix:')
  call matoutmdm (n, n - 1, ix, dx, wh)
endif

100 continue

return
end

!------------------------------

!   Second- and third-level subroutines.

function bound (n, nwdsm, wh)

!   This computes the norm bound using MPM arithmetic.

use mpmodule
implicit none
integer i, n, nwdsm
type (mp_real) wh(n,n), bound, t1

t1 = mpreald (0.d0, nwdsm)

do i = 1, n - 1
  t1 = max (t1, abs (wh(i,i)))
enddo

bound = 1.d0 / t1
return
end

function dplog10 (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
real (8) da, dplog10, t1
type (mp_real) a

call mpmdi (a, da, ia)
if (da == 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

function dplog10m (a)

!   For input MPM value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
real (8) da, dplog10m, t1
type (mp_real) a

call mpmdi (a, da, ia)
if (da == 0.d0) then
  dplog10m = -999999.d0
else
  dplog10m = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

subroutine decmd (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib
real (8) da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (mp_real) a

call mpmdi (a, da, ia)
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 < 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end

subroutine decmdm (a, b, ib)

!   For input MPM value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib
real (8) da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (mp_real) a

call mpmdi (a, da, ia)
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 < 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end

subroutine matoutdp (n1, n2, a)

!   This outputs the DP matrix a.

implicit none
integer i, j, n1, n2
real (8) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)
  write (6, 2) (a(i,j), j = 1, n2)
2 format (1p5d15.5)
enddo

return
end

subroutine matoutmd (n1, n2, ix, dx, a)

!   This outputs the MP matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
real (8) dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i6))
enddo

return
end

subroutine matoutmp (n1, n2, a)

!   This outputs the MP matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use mpmodule
implicit none
integer i, j, n1, n2
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpwrite (6, 80, 60, a(i,j))
  enddo
enddo

return
end

subroutine matoutmdm (n1, n2, ix, dx, a)

!   This outputs the MPM matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
real (8) dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmdm (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i6))
enddo

return
end

subroutine mxmdm (n1, n2, nwdsm, a, b)

!   This multiplies the DP square matrix a by the MPM matrix b, and the result
!   is placed in b.  n1, n2 are the matrix dimensions as indicated below.
!   This is performed in medium precision.
!   Input: n1, n2, nwdsm, a, b.
!   Output: b.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwdsm
real (8) a(n1,n1)
type (mp_real) b(n1,n2), c(n1,n2)

do j = 1, n2
  do i = 1, n1
    c(i,j) = mpreald (0.d0, nwdsm)
    
    do k = 1, n1
      c(i,j) = c(i,j) + mpprodd (b(k,j), a(i,k))
    enddo
  enddo
enddo

do j = 1, n2
  do i = 1, n1
    b(i,j) = c(i,j)
  enddo
enddo

return
end

subroutine mxmm (n1, n2, nwds, a, b)

!   This multiplies the MPM square matrix a by the MP matrix b, and the result
!   is placed in b.  n1, n2 are the matrix dimensions as indicated below.
!   This is performed in full precision.
!   Input: n1, n2, nwds, a, b.
!   Output: b.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwds
type (mp_real) a(n1,n1)
type (mp_real) b(n1,n2), c(n1,n2)

do j = 1, n2
  do i = 1, n1
    c(i,j) = mpreal (0.d0, nwds)

    do k = 1, n1
      c(i,j) = c(i,j) + mpreal (a(i,k), nwds) * b(k,j)
    enddo
  enddo
enddo

do j = 1, n2
  do i = 1, n1
    b(i,j) = c(i,j)
  enddo
enddo

return
end

subroutine qsortdp (n, a, ip)

!   This routine sorts the entries of the N-long DP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   Input: n, a.
!   Output: ip.

implicit none
integer i, iq, it, j, jq, jz, k, l, n
real (8) a(n), s0
integer ip(n), ik(50), jk(50)

do i = 1, n
  ip(i) = i
enddo

if (n == 1) return

k = 1
ik(1) = 1
jk(1) = n

130 continue

i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 < a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 > a(ip(l))) goto 180
enddo

j = i
goto 190

180 continue

j = l
if (i >= j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue

if (s0 >= a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 continue

k = k - 1
jz = 0
if (j == iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 continue

i = i + 1
if (i == jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz == 0) goto 220
if (j - iq >= jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 continue

if (k > 0) goto 130

return
end

subroutine qsortmpm (n, a, ip)

!   This routine sorts the entries of the N-long MPM vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   This is performed in medium precision.
!   Input: n, a, ip.
!   Output: ip.

use mpmodule
implicit none
integer i, iq, it, j, jq, jz, k, l, n
type (mp_real) a(n), s0
integer ip(n), ik(50), jk(50)

do i = 1, n
  ip(i) = i
enddo

if (n == 1) return

k = 1
ik(1) = 1
jk(1) = n

130 continue

i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 < a(ip(l))) goto 160
enddo

i = j
goto 190

160 continue

i = l

do l = j, i, -1
  if (s0 > a(ip(l))) goto 180
enddo

j = i
goto 190

180 continue

j = l
if (i >= j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue

if (s0 >= a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 continue

k = k - 1
jz = 0
if (j == iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 continue

i = i + 1
if (i == jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz == 0) goto 220
if (j - iq >= jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 continue

if (k > 0) goto 130

return
end
