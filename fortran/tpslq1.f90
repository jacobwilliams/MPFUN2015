program tpslq1

!   David H Bailey   19 Apr 2017

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2016 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program demonstrates pslq1, which performs the one-level standard
!   PSLQ algorithm on an input vector.  A variety of sample input vectors can
!   be generated as inputs to pslq1, as given in the parameters below.
!   For additional details, see:

!   David H. Bailey and David J. Broadhurst, "Parallel integer relation
!   detection: Techniques and applications," Mathematics of Computation,
!   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
!   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.

!   The pslq1 routine is 100% THREAD SAFE -- all requisite parameters and
!   arrays are passed through subroutine arguments.

!   These parameters are set in the parameter statement below.

!     idb   Debug level (0 - 3); default = 2.
!     itr   Number of trials when kq = 1, 2 or 3; default = 1.
!     lcx   Size of chr1 array; default = 64.
!           ***Must be a multiple of 64.
!           ***Must be >= 1 + size in digits of largest coefficient.
!     n     Integer relation vector length, = 1 + polynomial degree.
!           When kq = 0, set n = kr * ks + 1; this is the default.
!     kq    0: for the algebraic case [1, al, al^2, ... al^(n-1)], where
!             al = 3^(1/kr) - 2^(1/ks); this is the default.
!           1: for testing algebraic relations of a number read from a file.
!           2: for testing additive relations of numbers read from a file.
!           3: for testing multiplicative relations of numbers read from a file.
!           4: for custom input.
!     kr    Degree of root of 3 when kq = 0; default = 5.
!     ks    Degree of root of 2 when kq = 0; default = 6.
!     ndp   Full precision level in digits; default = 240.
!           ***Must be <= mpipl in module MPFUNF.
!           ***Must be >= ndpm + precision level required to find relation.
!     nep   Log10 of full precision epsilon for detections, typically 20 - ndp.
!           ***Must not be smaller than the accuracy of input data.  In other
!           words, if data is accurate to 10^(-200), then nep > -200.
!     rb    Log10 of max size (Euclidean norm) of acceptable relation;
!           default = 100. Run will abort if this is exceeded.
!     nwds  Precision level in words; default = int (ndp / mpdpw + 2).

use mpmodule
implicit none
real (8) d1, d2, rb, rm, second, tm, tm0, tm1
integer i, i1, idb, iq, itr, j, j1, k, kq, kr, ks, lcx, n, ndp, nep, nwds
parameter (idb = 2, itr = 1, kq = 0, kr = 5, ks = 6, lcx = 64, &
  n = kr * ks + 1, ndp = 250, nep = 20 - ndp, rb = 100.d0, &
  nwds = int (ndp / mpdpw + 2))
integer lnm(n)
character(1) chr1(lcx)
character(64) form4, nam(n), namx
real (8) dr(n), r1(n)
type (mp_real) al, eps, t1, t2, t3, r(n), x(n)
external second

!   Check to see if default precision is high enough.

if (ndp > mpipl) then
  write (6, '("Increase default precision in module MPFUNF.")')
  stop
endif

eps = mpreal (10.d0, nwds) ** nep
write (6, 1) itr, n, kq, kr, ks, rb, ndp, nep
1 format ('PSLQ1 Test Program'/ &
  'No. trials = ',i3,3x,'n =',i4,3x,'kq =',i2/ &
  'kr =',i2,3x,'ks =',i2,3x,'rb =',1p,d12.4/ &
  'Full precision level ndp =',i6,' digits'/ &
  'Full precision epsilon level nep = ',i6)

if (kq .eq. 1 .or. kq .eq. 2 .or. kq .eq. 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

do k = 1, itr
  write (6, 2) k
2 format (/'Start trial', i3)

  if (kq .eq. 0) then

!   This code generates al = 3^(1/kr) - 2^(1/ks).  al is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by al.

    al = mpnrt (mpreald (3.d0, nwds), kr) - mpnrt (mpreald (2.d0, nwds), ks)
  elseif (kq .eq. 1) then

!   Read an algebraic constant from a file.

    call mpread (11, al, nwds)
  elseif (kq .eq. 2) then

!   Read constants from a file for additive test.

    do i = 1, n
      call mpread (11, al, nwds)
      x(i) = al
      write (namx, '(''con'',i3.3)') i
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  elseif (kq .eq. 3) then

!   Read constants from a file for multiplicative test.

    do i = 1, n
      call mpread (11, al, nwds)
      x(i) = log (al)
      write (namx, '(''log(con'',i3.3,'')'')') i
      nam(i) = namx(1:11)
      lnm(i) = 11
    enddo
  elseif (kq .eq. 4) then

!   Produce X vector by a custom scheme.

  endif

!   If kq is 0 or 1, generate x = [1, al, al^2, ..., al^(n-1)].

  if (kq .eq. 0 .or. kq .eq. 1) then
    x(1) = mpreald (1.d0, nwds)
    nam(1) = '1'
    lnm(1) = 1

    do i = 2, n
      x(i) = al * x(i-1)
      write (namx, '(''al^'',i3)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  endif

!   Perform relation search.

  tm0 = second ()
  call pslq1 (idb, n, nwds, rb, eps, x, iq, r)
  tm1 = second ()

!   Output relation, if one was found.

  if (iq == 1) then
    write (6, 3)
3   format (/'Recovered relation: 0 =')
    form4 = '(64a1'
    i1 = 5

    do i = 65, lcx, 64
      form4(i1+1:i1+9) =  '/64a1'
      i1 = i1 + 5
    enddo

    form4(i1+1:i1+9) = '," * ",a)'
    i1 = i1 + 9
  
    do i = 1, n
      if (r(i) .ne. 0.d0) then
        call mpfform (r(i), lcx, 0, chr1)

        do j = 1, lcx
          if (chr1(j) /= ' ') goto 120
        enddo

120     continue

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

  write (6, 6) tm1 - tm0
6 format ('CPU Time =',f12.4)
enddo

stop
end

!------------------------------

!   The following code performs the one-level, standard PSLQ algorithm.
!   David H. Bailey     19 Apr 2017

subroutine pslq1 (idb, n, nwds, rb, eps, x, iq, r)

!   Arguments are as follows:
!     Name  Type    Description
!     idb   int*4   Debug flag (0-3); increasing idb produces more output.
!     n     int*4   Length of input vector x and output relation r.
!     nwds  int*4   Full precision level in words. This must be sufficient
!                     to recover the relation.
!     rb    real*8  Log10 of max size (Euclidean norm) of acceptable relation.
!                   NOTE: Run is aborted when this is exceeded.
!     eps   mp_real Tolerance for full precision relation detection. This is
!                     typically set to 2^(70-mpnbt*nwds) or so.
!     x     mp_real Input mp_real vector.
!     iq    int*4   Output flag: 0 (unsuccessful) or 1 (successful).
!     r     mp_real Output integer relation vector, if successful.

!   The following parameters are set in this routine:
!     ipi   int*4   Iteration print interval when idb >= 2; default = 100.
!     ipm   int*4   Iteration  check interval for MP iterations; default = 10.
!     itm   int*4   Maximum iteration count; default = 10^5.
!                   NOTE: Run is aborted when this is exceeded. If itm >= 10^7,
!                     change all "i7" in formats to "i8" or as high as needed.
use mpmodule
implicit none
integer i, idb, ipi, ipm, iq, it, itm, izm, j, j1, n, nwds, n1, n2, n3, n4
parameter (ipi = 100, ipm = 500, itm = 100000)
real (8) d1, d2, d3, d4, dplog10, rb, second, tm0, tm1, times(2)
type (mp_real) b(n,n), h(n,n), r(n), x(n), y(n), bound, dynrange, rn, eps, &
  t1, t2, t3, t4
external bound, dplog10, dynrange, second

!   Initialize.

if (idb .ge. 2) write (6, 1) n
1 format ('PSLQ1 integer relation detection: n =',i5)
iq = 0
it = 0
rn = mpreal (0.d0, nwds)

do i = 1, 2
  times(i) = 0.d0
enddo

if (idb .ge. 2) write (6, 2) it
2 format ('Iteration',i7,3x,'MP initialization')
tm0 = second ()
call initmp (idb, n, nwds, eps, b, h, x, y, izm)
tm1 = second ()
times(1) = tm1 - tm0
if (izm .ne. 0) goto 110

!   MP iterations.

if (idb .ge. 2) write (6, 3) it
3 format ('Iteration',i7,3x,'Start MP iterations')

100 continue

!   Perform one MP iteration.

it = it + 1
if (idb .eq. 3 .or. idb .ge. 2 .and. mod (it, ipi) .eq. 0) write (6, 4) it
4 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, it, n, nwds, eps, b, h, y, izm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Test conditions on itermp output flag izm:
!   0: MP update was uneventful; periodically output min, max, norm; continue.
!   1: A small value was found in MP update; go output relation.
!   2: Precision is exhausted; quit.

if (izm .eq. 0) then

!   Periodically output min, max and norm.

  if (mod (it, ipm) .eq. 0) then

!   Find min and max absolute value of y vector.

    call minmax (n, nwds, y, t1, t2)
    if (idb >= 2) then
      call decmd (t1, d1, n1)
      call decmd (t2, d2, n2)
      write (6, 5) it, d1, n1, d2, n2
5     format ('Iteration',i7,3x,'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6)
    endif

!   Compute norm bound.

    t1 = bound (n, nwds, h)
    rn = max (rn, t1)
    if (idb .ge. 2) then
      call decmd (t1, d1, n1)
      call decmd (rn, d2, n2)
      write (6, 6) it, d1, n1, d2, n2
6     format ('Iteration',i7,3x,'Norm bound =',f11.6,'e',i5,4x,'Max. bound =', &
        f11.6,'e',i5)
    endif

!   Check if iteration limit or norm bound limit is exceeded; if so, quit.

    if (it .gt. itm) then
      if (idb .ge. 1) write (6, 7) itm
7     format ('Iteration limit exceeded',i7)
      goto 120
    endif
    if (dplog10 (rn) .gt. rb) then
      if (idb .ge. 1) write (6, 8) rb
8     format ('Norm bound limit exceeded.',1p,d15.6)
      goto 120
    endif
  endif
  goto 100
elseif (izm .eq. 1) then
  goto 110
elseif (izm .eq. 2) then
  goto 120
endif

110 continue

!   A relation has been detected.

tm0 = second ()
t1 = mpreald (1.d300, nwds)
t2 = mpreal (0.d0, nwds)
t3 = mpreal (0.d0, nwds)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(i,j1)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.d0 ** n3

!   Output the final norm bound and other info.

if (idb .ge. 1) then
  t4 = bound (n, nwds, h)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmd (t3, d3, n3)
  call decmd (t4, d4, n4)
  write (6, 9) it, d1, n1, d2, n2, d4, n4
9 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'e',i5,f11.6,'e',i5/'Max. bound =',f11.6,'e',i5)
  write (6, 10) j1, d3, n3, d1, n1
10 format ('Index of relation =',i4,3x,'Norm =',f11.6,'e',i5,3x, &
  'Residual =',f11.6,'e',i5)
endif

!   If run was successful, set iq = 1.

if (dplog10 (t3) .le. rb) then
  iq = 1
else
  if (idb .ge. 2) write (6, 11)
11 format ('Relation is too large.')
endif

120 continue

!   Output CPU run times and return.

if (idb .ge. 2) write (6, 12) times
12 format ('CPU times:'/(5f12.2))

return
end

!------------------------------

!   First-level subroutines.

!   First-level subroutines.

subroutine minmax (n, nwds, y, y1, y2)

!   This returns min|y_k| and max|y_k| using MP precision.
!   Input: n, nwds, y.
!   Output: y1, y2.

use mpmodule
implicit none
integer i, n, n1, nwds
real (8) d1
type (mp_real) dynrange, t1, t2, t3
type (mp_real) y(n), y1, y2

t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (mpreal (y(i), nwds))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

y1 = t1
y2 = t2
return
end

subroutine initmp (idb, n, nwds, eps, b, h, x, y, izm)

!   This initializes MP arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nwds, eps.
!   Output: b, h, x, y, izm.

use mpmodule
implicit none
integer i, idb, izm, j, k, n, nwds, n1
real (8) d1
integer ix(n)
real (8) dx(n)
type (mp_real) b(n,n), h(n,n), s(n), x(n), y(n), eps, t1

if (idb .ge. 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif
izm = 0

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = mpreal (0.d0, nwds)
  enddo

  b(j,j) = mpreal (1.d0, nwds)
enddo

t1 = mpreal (0.d0, nwds)

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
    h(i,j) = mpreal (0.d0, nwds)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
!      t1 = y(i) * t1
      h(i,j) = - y(i) * t1
  enddo
enddo

!  Perform reduction on h, updating y, b and h.

do i = 2, n
  do j = i - 1, 1, -1
    t1 = anint (h(i,j) / h(j,j))
    if (t1 .ne. 0.d0) then
      y(j) = y(j) + t1 * y(i)

      do k = i, n
        b(k,j) = b(k,j) + t1 * b(k,i)
      enddo

      do k = 1, j
        h(i,k) = h(i,k) - t1 * h(j,k)
      enddo
    endif
  enddo
enddo

t1 = mpreald (1.d300, nwds)

do i = 1, n
  t1 = min (t1, abs (y(i)))
enddo

if (t1 .le. eps) then
  if (idb .ge. 2) then
    call decmd (t1, d1, n1) 
    write (6, 2) 0, d1, n1
2   format ('Iteration',i8,3x,'initmp: Small value in y =',f11.6,'e',i5)
  endif
  izm = 1
endif

if (idb .ge. 3) then
  write (6, 3)
3 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 4)
4 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine itermp (idb, it, n, nwds, eps, b, h, y, izm)

!   This performs one iteration of the PSLQ algorithm using MP arithmetic.
!   This is performed in medium precision.
!   Input: idb, it, n, nwds, eps.
!   Output: b, h, y, izm.

use mpmodule
implicit none
integer i, idb, im, im1, it, izm, j, j1, k, n, nwds, n1, ntl
parameter (ntl = 72)
real (8) d1, d2
integer ix(n)
real (8) dx(n)
type (mp_real) b(n,n), h(n,n), y(n)
type (mp_real) eps, gam, teps, t1, t2, t3, t4

teps = 2.d0 ** ntl * eps
izm = 0
t1 = mpreal (0.d0, nwds)
gam = sqrt (mpreald (4.d0, nwds) / mpreald (3.d0, nwds))

!   Find  im = i  such that  gam^i * |h(i,i)|  is maximized.

do i = 1, n - 1
  t2 = gam ** i * abs (h(i,i))
  if (t2 .gt. t1) then
    im = i
    t1 = t2
  endif
enddo

im1 = im + 1

! Exchange the im and im+1 entries of y, rows of h, and columns of b.

t1 = y(im)
y(im) = y(im1)
y(im1) = t1

do i = 1, n
  t1 = b(i,im)
  b(i,im) = b(i,im1)
  b(i,im1) = t1
enddo

do i = 1, n - 1
  t1 = h(im,i)
  h(im,i) = h(im1,i)
  h(im1,i) = t1
enddo

!  Update h with permutation produced above.

if (im .le. n - 2) then
  t1 = h(im,im)
  t2 = h(im,im1)
  t3 = sqrt (t1 ** 2 + t2 ** 2)
  t1 = t1 / t3
  t2 = t2 / t3

  do i = im, n
    t3 = h(i,im)
    t4 = h(i,im1)
    h(i,im) = t1 * t3 + t2 * t4
    h(i,im1) = - t2 * t3 + t1 * t4
  enddo
endif

t2 = mpreal (0.d0, nwds)

!  Perform reduction on h, updating y, a, b and h.

do i = im1, n
  j1 = min (i - 1, im1)

  do j = j1, 1, -1
    t1 = anint (h(i,j) / h(j,j))
    if (t1 .ne. 0.d0) then
      y(j) = y(j) + t1 * y(i)

      do k = 1, n
        b(k,j) = b(k,j) + t1 * b(k,i)
      enddo

      do k = 1, j
        h(i,k) = h(i,k) - t1 * h(j,k)
      enddo
    endif
  enddo
enddo

!   Find the min of y.

t1 = mpreald (1d300, nwds)

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!  Find the largest entry of b in the same column as the smallest y.

t2 = mpreal (0.d0, nwds)

do i = 1, n
  t2 = max (t2, abs (b(i,j1)))
enddo

if (t1 .le. t2 * teps) then
  if (idb .ge. 2) then
    call decmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'e',i5)
  endif
  if (t1 .le. t2 * eps) then
    izm = 1
  else
    if (idb .ge. 1) write (6, 2) it
2   format ('Iteration',i7,3x,'itermp: Precision exhausted.')
    izm = 2
  endif
endif

if (idb .ge. 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, ix, dx, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

!------------------------------

!   Second- and third-level subroutines.

function bound (n, nwds, h)

!   This computes the norm bound using DP arithmetic.

use mpmodule
implicit none
integer i, n, nwds
type (mp_real) bound, h(n,n), t1, t2

t1 = mpreal (0.d0, nwds)

do i = 1, n - 1
  t1 = max (t1, abs (h(i,i)))
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
if (da .eq. 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da)) + ia * log10 (2.d0)
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
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

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
2 format (4(f13.8,'e',i5))
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

