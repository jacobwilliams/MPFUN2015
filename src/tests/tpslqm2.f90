program tpslqm2

!   David H Bailey   19 Apr 2017

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2016 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program demonstrates pslqm2, which performs the one-level multipair
!   PSLQ algorithm on an input vector.  A variety of sample input vectors can
!   be generated as inputs to pslqm2, as given in the parameters below.  The
!   pslqm2 routine is suitable for relations up to degree 100 or so; above this
!   level pslqm3 should be used for significantly better performance.
!   For additional details, see:

!   David H. Bailey and David J. Broadhurst, "Parallel integer relation
!   detection: Techniques and applications," Mathematics of Computation,
!   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
!   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.

!   The pslqm2 routine is 100% THREAD SAFE -- all requisite parameters and
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
!     kr    Degree of root of 3 when kq = 0; default = 7.
!     ks    Degree of root of 2 when kq = 0; default = 8.
!     ndp   Full precision level in digits; default = 650.
!           ***Must be <= mpipl in module MPFUNF.
!           ***Must be >= precision level required to find relation.
!     nep   Log10 of full precision epsilon for detections; default = 20 - ndp.
!           ***Must not be smaller than the accuracy of input data.  In other
!           words, if data is accurate to 10^(-200), then nep > -200.
!     rb    Log10 of max size (Euclidean norm) of acceptable relation; 
!           default = 100. Run will abort if this is exceeded.
!     nwds  Precision level in words; default = int (ndp / mpdpw + 2).

use mpmodule
implicit none
real (8) d1, d2, rb, second, tm, tm0, tm1
integer i, i1, idb, iq, itr, j, j1, k, kq, kr, ks, lcx, n, ndp, nep, nwds
parameter (idb = 2, itr = 1, kq = 0, kr = 7, ks = 8, lcx = 64, &
  n = kr * ks + 1, ndp = 650, nep = 20 - ndp, rb = 100.d0, &
  nwds = int (ndp / mpdpw + 2))
integer lnm(n)
character(1) chr1(lcx)
character(64) form4, nam(n), namx
real (8) r1(n)
type (mp_real) al, eps, r(n), x(n)
external second

!   Check to see if default precision is high enough.

if (ndp > mpipl) then
  write (6, '("Increase default precision in module MPFUNF.")')
  stop
endif

eps = mpreal (10.d0, nwds) ** nep
write (6, 1) itr, n, kq, kr, ks, rb, ndp, nep
1 format ('PSLQM2 Test Program'/ &
  'No. trials = ',i3,3x,'n =',i4,3x,'kq =',i2/ &
  'kr =',i2,3x,'ks =',i2,3x,'rb =',1p,d12.4/ &
  'Full precision level ndp =',i6,' digits'/ &
  'Full precision epsilon level nep = ',i6)

if (kq == 1 .or. kq == 2 .or. kq == 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

do k = 1, itr
  write (6, 2) k
2 format (/'Start trial', i3)

  if (kq == 0) then

!   This code generates al = 3^(1/kr) - 2^(1/ks).  al is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by al.

    al = mpnrt (mpreald (3.d0, nwds), kr) - mpnrt (mpreald (2.d0, nwds), ks)
  elseif (kq == 1) then

!   Read an algebraic constant from a file.

    call mpread (11, al, nwds)
  elseif (kq == 2) then

!   Read constants from a file for additive test.

    do i = 1, n
      call mpread (11, al, nwds)
      x(i) = al
      write (namx, '(''con'',i3.3)') i
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  elseif (kq == 3) then

!   Read constants from a file for multiplicative test.

    do i = 1, n
      call mpread (11, al, nwds)
      x(i) = log (al)
      write (namx, '(''log(con'',i3.3,'')'')') i
      nam(i) = namx(1:11)
      lnm(i) = 11
    enddo
  elseif (kq == 4) then

!   Produce X vector by a custom scheme.

  endif

!   If kq is 0 or 1, generate x = [1, al, al^2, ..., al^(n-1)].

  if (kq == 0 .or. kq == 1) then
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
  call pslqm2 (idb, n, nwds, rb, eps, x, iq, r)
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

!   The following code performs the two-level, multi-pair PSLQ algorithm.
!   David H. Bailey     19 Apr 2017

subroutine pslqm2 (idb, n, nwds, rb, eps, x, iq, r)

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
!     itm   int*4   Maximum iteration count; default = 10^6.
!                   NOTE: Run is aborted when this is exceeded. If itm >= 10^7,
!                     change all "i7" in formats to "i8" or as high as needed.
!     nsq   int*4   Size of tables used in iterdp and itermpw; default = 8.
!     deps  real*8  Tolerance for dynamic range check; default = 1d-10.

use mpmodule
implicit none
integer i, idb, imq, ipi, ipm, iq, it, itm, its, izd, izm, j, j1, n, n1, n2, &
  n3, nsq, nwds
parameter (ipi = 100, ipm = 10, itm = 1000000, nsq = 8)
real (8) bounddp, d1, d2, d3, d4, rb, rn, second, tm0, tm1, &
  times(4)
real (8) da(n,n), db(n,n), dh(n,n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsyq(n,nsq), dy(n), dsy(n), deps
parameter (deps = 1d-10)
type (mp_real) b(n,n), h(n,n), syq(n,nsq), r(n), x(n), &
  y(n), dynrange, eps, t1, t2, t3
external bounddp, dynrange, second

!   Initialize.

if (idb >= 2) write (6, 1) n
1 format ('PSLQM2 integer relation detection: n =',i5)
iq = 0
it = 0
its = 0
imq = 0
rn = 0.d0

do i = 1, 4
  times(i) = 0.d0
enddo

if (idb >= 2) write (6, 2) it
2 format ('Iteration',i7,3x,'MP initialization')
tm0 = second ()
call initmp (idb, n, nsq, nwds, b, h, syq, x, y)
tm1 = second ()
times(1) = tm1 - tm0

100 continue

!   Check if dynamic range of y vector is too great for DP iterations
!   (which is often the case at the start of the run). If so, then
!   perform MP iterations instead of DP iterations.

t1 = dynrange (n, nwds, y)
if (idb >= 2) then
  call decmd (t1, d1, n1)
  write (6, 3) it, d1, n1
3 format ('Iteration',i7,3x,'Min/max ratio in y =',f11.6,'e',i6)
endif
if (t1 < mpreald (deps, nwds)) then
  goto 120
endif

!   DP iterations.
!   Initialize DP arrays from MP arrays.

if (idb >= 3) write (6, 4) it
4 format ('Iteration',i7,3x,'Start DP iterations')

call initdp (idb, n, nsq, da, db, dh, dy, dsyq, h, y)

!   Save DP arrays and iteration count.

call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
its = it

!   Perform an LQ decomposition on dh.

call lqdp (n, n - 1, dh)

110 continue

!   Perform one DP iteration.

it = it + 1
if (idb >= 3 .or. idb >= 2 .and. mod (it, ipi) == 0) write (6, 5) it
5 format ('Iteration',i7)
tm0 = second ()
call iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Test conditions on iterdp output flag izd:
!   0: Iteration was uneventful; periodically save arrays and continue.
!   1: Relation found or DP precision exhausted; perform MP update.
!   2: Very large value appeared in DA or DB; abort DP iter and do MP update.

if (izd == 0) then
  if (mod (it - its, ipm) == 0) then
    call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
    its = it
  endif
  goto 110
else

!   Check if DP iteration was aborted above; if so, then revert to previous data.

  if (izd == 2) then

!   DP iteration was aborted -- revert to previous data.

    it = its
    call savedp (n, dsa, dsb, dsh, dsy, da, db, dh, dy)
  endif

!   Update the MP arrays from the DP arrays.

  if (idb >= 2) write (6, 6) it
6 format ('Iteration',i7,3x,'MP update')
  tm0 = second ()
  call updtmp (idb, it, n, nwds, da, db, eps, b, h, y, izm)
  tm1 = second ()
  times(3) = times(3) + (tm1 - tm0)

!   Compute norm bound using DP.

  do j = 1, n - 1
    do i = 1, n
      dh(i,j) = h(i,j)
    enddo
  enddo
    
  call lqdp (n, n - 1, dh)
  d3 = bounddp (n, dh)
  if (d3 == -1.d0) goto 150
  rn = max (rn, d3)
  if (idb >= 2) then
    write (6, 7) it, d3, rn
7   format ('Iteration',i7,3x,'Norm bound =',1p,d15.6,4x,'Max. bound =', &
      1p,d15.6)
  endif

!   Check if iteration limit or norm bound limit is exceeded; if so, quit.

  if (it > itm) then
    if (idb >= 1) write (6, 8) itm
8   format ('Iteration limit exceeded',i7)
    goto 150
  endif
  if (log10 (rn) > rb) then
    if (idb >= 1) write (6, 9) rb
9   format ('Norm bound limit exceeded.',1p,d15.6)
    goto 150
  endif

!   Test conditions on updtmp output flag izm:
!   0: MP update was uneventful; goto tag 100 above, provided izd = 0;
!        if izd = 2, then go perform MP iterations.
!   1: A small value was found in MP update; go output relation.
!   2: Precision is exhausted; quit.

  if (izm == 0) then
    if (izd == 2) then
      goto 120
    else
      goto 100
    endif
  elseif (izm == 1) then
    goto 140
  elseif (izm == 2) then
    goto 150
  endif
endif

120 continue

!   MP iterations.

if (idb >= 2) write (6, 10) it
10 format ('Iteration',i7,3x,'Start MP iterations')

!   Perform an LQ decomposition on h.

tm0 = second ()
call lqmp (n, n - 1, nwds, h)
tm1 = second ()
times(1) = times(1) + (tm1 - tm0)

130 continue

!   Perform one MP iteration.

it = it + 1
if (idb >= 2) write (6, 11) it
11 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, it, n, nsq, nwds, eps, b, h, syq, y, imq, izm)
tm1 = second ()
times(4) = times(4) + (tm1 - tm0)

!   Test conditions on itermp output flag izm:
!   0: Iteration was uneventful; periodically check for DP; continue.
!   1: A small value was found; go output relation.
!   2: Precision is exhausted; quit.

if (izm == 0) then

!   Periodically check to see if DP iterations can be resumed (but perform
!   at least IPM iterations in MP).

  if (mod (it - its, ipm) == 0) then
    t1 = dynrange (n, nwds, y)
    if (idb >= 2) then
      call decmd (t1, d1, n1)
      write (6, 12) it, d1, n1
12    format ('Iteration',i7,3x,'Min/max ratio in y =',f11.6,'e',i6)
    endif
    if (t1 > mpreald (deps, nwds)) then
      if (idb >= 2) write (6, 13) it
13    format ('Iteration',i7,3x,'Return to DP iterations')
      goto 100
    endif
  endif
  goto 130
elseif (izm == 1) then
  goto 140
elseif (izm == 2) then
  goto 150
endif

!   A relation has been detected.  Output the final norm bound and other info.

140 continue

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
  r(i) = b(j1,i)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.d0 ** n3

!   Output the final norm bound and other info.

if (idb >= 1) then
  do j = 1, n - 1
    do i = 1, n
      dh(i,j) = h(i,j)
    enddo
  enddo

  call lqdp (n, n - 1, dh)
  d4 = bounddp (n, dh)
  rn = max (rn, d4)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  write (6, 14) it, d1, n1, d2, n2, rn
14 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'e',i5,f11.6,'e',i5/'Max. bound =',1p,d15.6)
  write (6, 15) j1, d3, d1, n1
15 format ('Index of relation =',i4,3x,'Norm =',1p,d15.6,3x, &
  'Residual =',0p,f11.6,'e',i5)
endif

!   If run was successful, set iq = 1.

if (log10 (d3) <= rb) then
  iq = 1
else
  if (idb >= 2) write (6, 16)
16 format ('Relation is too large.')
endif

150 continue

!   Output CPU run times and return.

if (idb >= 2) write (6, 17) times
17 format ('CPU times:'/(4f12.2))

return
end

!------------------------------

!   First-level subroutines.

function dynrange (n, nwds, y)

!   This returns the dynamic range of y, i.e., ratio of min|y_k| to max|y_k|,
!   using MP precision.

use mpmodule
implicit none
integer i, n, n1, nwds
real (8) d1
type (mp_real) dynrange, t1, t2, t3
type (mp_real) y(n)

t1 = mpreald (1.d300, nwds)
t2 = mpreald (0.d0, nwds)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (mpreal (y(i), nwds))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

dynrange = t1 / t2
return
end

subroutine initdp (idb, n, nsq, da, db, dh, dy, dsyq, h, y)

!   This initializes the DP arrays from the MP arrays.
!   This is performed in four-word precision.
!   Input:  idb, n, nsq, h, y.
!   Output: da, db, dh, dy, dsyq.

use mpmodule
implicit none
integer i, idb, j, n, n1, nsq, nwx
real (8) da(n,n), db(n,n), dh(n,n), dy(n), dsyq(n,nsq)
type (mp_real) h(n,n), y(n), t1, t2, t3, t4
parameter (nwx = 4)

t2 = mpreal (0.d0, nwx)

!   Find the max absolute value in the y vector.

do i = 1, n
  t2 = max (t2, mpreal (abs (y(i)), nwx))
enddo

!   Set dy to be the scaled y vector.

t1 = 1.d0 / t2

do i = 1, n
  dy(i) = t1 * mpreal (y(i), nwx)
enddo

!   Find the maximum absolute value of the h matrix diagonals.

t2 = mpreald (0.d0, nwx)

do j = 1, n - 1
  t2 = max (t2, mpreal (abs (h(j,j)), nwx))
enddo

!   Set dh to be the scaled h matrix.

t1 = 1.d0 / t2

do j = 1, n - 1
  do i = 1, n
    dh(i,j) = t1 * mpreal (h(i,j), nwx)
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

if (idb >= 3) then
  write (6, 2)
2 format ('initdp: Scaled dy vector:')
  call matoutdp (1, n, dy)
  write (6, 3)
3 format ('initdp: Scaled dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine initmp (idb, n, nsq, nwds, b, h, syq, x, y)

!   This initializes MP arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nsq, nwds.
!   Output: b, h, syq, x, y.

use mpmodule
implicit none
integer i, idb, j, n, nsq, nwds
integer ix(n)
real (8) dx(n)
type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), x(n), y(n), t1

if (idb >= 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif

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
    h(i,j) = - y(i) * t1
  enddo
enddo

!   Zero the syq array.

do j = 1, nsq
  do i = 1, n
    syq(i,j) = mpreal (0.d0, nwds)
  enddo
enddo

if (idb >= 3) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine iterdp (idb, it, n, nsq, da, db, dh, dsyq, dy, imq, izd)

!   This performs one iteration of the PSLQ algorithm using DP arithmetic.
!   Input: idb, it, n, nsq, da, db, dh, dsyq, dy.
!   Output: da, db, dh, dsyq, dy, imq, izd.

!   NOTE: Parameter tmx2 = 2^52, not 2^53, so as to ensure that values > 2^53
!   never arise, even as intermediate values, in the update loop below.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izd, j, j1, j2, k, mpr, mq, n, &
  nsq
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
1 format ('Iteration',i7,3x,'iterdp: Small value in dy =',1pd15.6)
  izd = 1
endif

if (t2 > tmx1 .and. t2 <= tmx2) then
  if (idb >= 3) write (6, 2) it, t2
2 format ('Iteration',i7,3x,'iterdp: Large value in da or db =',1pd15.6)
  izd = 1
elseif (t2 > tmx2) then
  if (idb >= 2) write (6, 3) it, t2
3 format ('Iteration',i7,3x,'iterdp: Very large value in da or db =',1pd15.6)
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
4   format ('Iteration',i7,3x,'iterdp: Duplicate found, j =',i6)
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

if (idb >= 3) then
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

subroutine itermp (idb, it, n, nsq, nwds, eps, b, h, syq, y, imq, izm)

!   This performs one iteration of the PSLQM algorithm using MP arithmetic.
!   This is performed in medium precision.
!   Input: idb, it, n, nsq, nwds, eps, imq.
!   Output: b, h, syq, y, imq, izm.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izm, j, j1, j2, k, mpr, mq, n, n1, &
  nsq, ntl, nwds
parameter (ntl = 72)
real (8) d1
type (mp_real) eps, gam, t1, t2, t3, t4, teps
integer ip(n), ir(n), is(n)
real (8) dx(n)
type (mp_real) b(n,n), h(n,n), q(n), syq(n,nsq), t(n,n), y(n)

teps = 2.d0 ** ntl * eps
izm = 0
mpr = nint (0.4d0 * n)
gam = sqrt (mpreald (4.d0, nwds) / mpreald (3.d0, nwds))

!   Compute q vector = {gam^i * |h(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  q(i) = gam ** i * abs (h(i,i))
enddo

call qsortmp (n - 1, q, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest q(i).

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

!   Exchange the pairs of entries of y, and rows of b and h.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = y(im)
  y(im) = y(im1)
  y(im1) = t1

  do i = 1, n
    t1 = b(im,i)
    b(im,i) = b(im1,i)
    b(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = h(im,i)
    h(im,i) = h(im1,i)
    h(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in h.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
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
enddo

!   Perform reduction on h, using the diagonal scheme.  Multipliers are
!   saved in the t array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      h(ij,j) = h(ij,j) - t(ij,k) * h(k,j)
    enddo

    t(ij,j) = anint (h(ij,j) / h(j,j))
    h(ij,j) = h(ij,j) - t(ij,j) * h(j,j)
  enddo
enddo

!   Update y, using the t array.  Find min absolute value of y.

t1 = abs (y(n))
j1 = n

do j = 1, n - 1
  do i = j + 1, n
    y(j) = y(j) + t(i,j) * y(i)
  enddo

  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!   Update b, using the t array.

do k = 1, n
  do j = 1, n - 1
    do i = j + 1, n
      b(j,k) = b(j,k) + t(i,j) * b(i,k)
    enddo
  enddo
enddo

!  Find the largest entry of b in the same row as the smallest y.

t2 = mpreal (0.d0, nwds)

do i = 1, n
  t2 = max (t2, abs (b(j1,i)))
enddo

if (t1 <= t2 * teps) then
  if (idb >= 2) then
    call decmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'e',i5)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 2) it
 2  format ('Iteration',i7,3x,'itermp: Precision exhausted.')
    izm = 2
  endif
endif

!   Compare the y vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = mpreal (0.d0, nwds)

  do i = 1, n
    t1 = max (t1, abs (y(i) - syq(i,j)))
  enddo

  if (t1 <= t2 * teps) then
    if (idb >= 2) write (6, 3) it, j
3   format ('Iteration',i7,3x,'itermp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector y in the table syq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  syq(i,k) = y(i)
enddo

if (idb >= 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
  call matoutmd (1, n, ip, dx, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, ip, dx, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, ip, dx, h)
endif

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

subroutine lqmp (n, m, nwds, h)

!   This performs an LQ decomposition on the MP matrix h.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   Input: n, m, nwds, h.
!   Output: h.

use mpmodule
implicit none
integer i, j, l, lup, m, ml, n, nwds
type (mp_real) h(n,m), nrmxl, one, t, zero

zero = mpreal (0.d0, nwds)
one = mpreal (1.d0, nwds)
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
    h(i,j) = mpreal (0.d0, nwds)
  enddo
enddo

return
end

subroutine savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)

!   This saves the arrays dy, da, db, dh in case dp iterations must be aborted.
!   A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
!   exchanged, serves to restore these arrays.

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

subroutine updtmp (idb, it, n, nwds, da, db, eps, b, h, y, izm)

!   This updates the MP arrays from the DP arrays.
!   Input: idb, it, n, nwds, wa, wb, eps, b, h, y.
!   Output: b, h, y, izm.

use mpmodule
implicit none 
integer i, i1, idb, it, izm, n, n1, n2, ntl, nwds
parameter (ntl = 72)
integer ix(n)
real (8) da(n,n), db(n,n), dx(n), d1, d2
type (mp_real) b(n,n), h(n,n), y(n), eps, t1, t2, teps

if (idb >= 2) write (6, 1) it
1 format ('Iteration',i7,3x,'updtmp: MP update from DP arrays.')
teps = 2.d0 ** ntl * eps
izm = 0

!   Update y with db.

call mxmdm (n, 1, nwds, db, y)
i1 = 0
t1 = mpreald (1d300, nwds)
t2 = mpreal (0.d0, nwds)

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
  write (6, 2) it, d1, n1, d2, n2
2  format ('Iteration',i7,3x,'updtmp: Min, max of y =',f11.6,'e',i5, &
     f11.6,'e',i5)
endif

!   Update b with db.

call mxmdm (n, n, nwds, db, b)

!   Update h with da.  There is no need to perform a LQ decomposition on h.

call mxmdm (n, n - 1, nwds, da, h)

!   Find the largest entry of b in the same row as the smallest y.

t2 = mpreal (0.d0, nwds)

do i = 1, n
  t2 = max (t2, abs (b(i1,i)))
enddo

if (t1 <= t2 * teps) then
  if (idb >= 2) then
    call decmd (t1, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i7,3x,'updtmp: Small value in y =',f11.6,'e',i5)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 4) it
4   format ('Iteration',i7,3x,'updtmp: Precision exhausted.')
    izm = 2
  endif
endif

if (idb >= 3) then
  write (6, 5)
5 format ('updtmp: Updated y:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 6)
6 format ('updtmp: Updated b matrix:')
  call matoutmd (n, n, ix, dx, b)
  write (6, 7)
7 format ('updtmp: Updated h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

!------------------------------

!   Second- and third-level subroutines.

function bounddp (n, dh)

!   This computes the norm bound using DP arithmetic.

implicit none
integer i, n
real (8) bounddp, t1
real (8) dh(n,n)

t1 = 0.d0

do i = 1, n - 1
  t1 = max (t1, abs (dh(i,i)))
enddo

if (t1 < 1.d-300) then
  write (6, 1)
1 format ('bound: h matrix too small--use 1-level or 3-level pslqm program.')
  bounddp = -1.d0
else
  bounddp = 1.d0 / t1
endif

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
2 format (4(f13.8,'D',i5))
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

subroutine mxmdm (n1, n2, nwds, a, b)

!   This multiplies the DP square matrix a by the MP matrix b, and the result
!   is placed in b.  n1, n2 are the matrix dimensions as indicated below.
!   Input: n1, n2, nwds, a, b.
!   Output: b.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwds
real (8) a(n1,n1)
type (mp_real) b(n1,n2), c(n1)

do j = 1, n2
  do i = 1, n1
    c(i) = mpreal (0.d0, nwds)
    
    do k = 1, n1
       c(i) = c(i) + mpprodd (b(k,j), a(i,k))
    enddo
  enddo

  do i = 1, n1
    b(i,j) = c(i)
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
integer ip(n), ik(50), jk(50)
real (8) a(n), s0

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

subroutine qsortmp (n, a, ip)

!   This routine sorts the entries of the N-long MP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
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
