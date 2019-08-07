program tquadgs

!   David H. Bailey     15 Sep 2018

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2018 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   NOTE: This program runs very long (800 seconds with MPFUN-MPFR or 2800
!   seconds with MPFUN-Fort, on the author's computer), compared with tquad.f90
!   and the other sample application programs. For this reason it is NOT
!   included in the script mpfun-tests.scr.

!   This program demonstrates the quadrature routine 'quadgs', which employs
!   Gaussian quadrature. The function quadgs is suitable to integrate a
!   function that is continuous, infinitely differentiable and integrable on an
!   open interval containing the given finite interval. It can also be used for
!   integrals on the entire real line, by making a variable transformation,
!   which can be achieved by merely computing a modified set of weights and
!   abscissas -- see subroutine initqgs and the examples below.

!   The quadrature computation proceeds level-by-level, with the computational
!   cost (and, often, the accuracy) approximately doubling with each iteration,
!   until a target accuracy (500-digit precision in this case), based on an
!   accuracy estimate computed by the program, is achieved. In most cases
!   quadgs is significantly faster than quadts (tanh-sinh quadrature) for
!   regular functions on finite intervals. However, quadgs is not effective
!   in cases where the function or one of its higher-order derivatives has a
!   singularity at one or both endpoints. Also, the computational cost of
!   generating weight-abscissa data in quadgs is MUCH greater than for quadts,
!   quades or quadss; on the other hand, these data can be computed once and
!   stored in a file -- see below. For additional details on Gaussian
!   quadrature versus tanh-sinh schemes, see:

!   David H. Bailey, Xiaoye S. Li and Karthik Jeyabalan, "A comparison of
!   three high-precision quadrature schemes," Experimental Mathematics, 
!   vol. 14 (2005), no. 3, pg. 317-329, preprint available at
!   http://www.davidhbailey.com/dhbpapers/quadrature-em.pdf.

!   The function to be integrated must be defined in an external function
!   subprogram (see samples below), and the name of the function must be
!   included in a "type (mp_real)" and an "external" statement. Prior to
!   calling the quadrature routine, the initialization routine initgs must
!   be called to initialize the weight-abscissa arrays wkgs, xkgs, wkgu and
!   xkgu, or else these arrays must be read from a file -- see idata  below.

!   All of these routines are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Specific instructions:

!   First call initqgs with the arrays wkgs, xkgs, xkgu and xkgu as the last
!   four arguments. The initialization routine initqgs may run for a long
!   time, depending on precision: for nq1 = 11 and ndp1 = 500, initqgs runs
!   12 minutes on the author's computer. In the code below, one has the option
!   of calling initqgs and then writing wkgs, xkgs, wkgu and xkgu to a file,
!   or else reading the previously computed arrays from a file -- see the
!   parameter idata. 

!   For finite interval problems, call quadgs with the arrays wkgs and xkgs 
!   as the last two arguments. Set x1 and x2, the limits of integration,
!   in executable statements, using high precision (nwds2 words) if possible.
!   Note in the examples below, x1 and x2 and set using high-precision values
!   of zero, one and pi (the variables zero, one and mppic). Using high
!   precision values for x1 and x2 facilitates more accurate argument scaling
!   in quadgs and, if necessary, in the function routines. See examples below.

!   For problems over the entire real line, call quadgs with wkgu and xkgu as
!   the last two arguments, and set x1 = -one and x2 = one. See examples below.

!   The problems performed below are a subset of the problems performed in
!   tquad.f90; the problem numbers here correspond to those in tquad.f90.

!   The following inputs are set in the parameter statement below:
!   idata  0: compute abscissas and weights from scratch and write to file
!             before performing quadrature problems.
!          1: read abscissas and weights from a file.
!          By default, idata = 0; the file name is "gauss-mm-nnnn.dat", where
!          "mm" is the value of nq1 and "nnnn" is the value of ndp1.
!   ndp1   Primary ("low") precision in digits; this is the working precision
!          and also the target accuracy of quadrature results.
!   ndp2   Secondary ("high") precision in digits. By default, ndp2 = 2*ndp1.
!   neps1  Log10 of the primary tolerance. By default, neps1 = - ndp1.
!   neps2  Log10 of the secondary tolerance. By default, neps2 = -ndp2.
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!          (possibly doubles) the number of accurate digits in the result,
!          but also roughly doubles the run time. nq1 must be at least 3.
!   nq2    Space parameter for wk and xk arrays in the calling program.  By
!          default it is set to 6 * 2^nq1 + 100.
!   nwds1  Low precision in words. By default nwds1 = int (ndp1 / mpdpw + 2).
!   nwds2  High precision in words. By default nwds2 = int (ndp2 / mpdpw + 2).

use mpmodule
implicit none
integer i, i1, i2, idata, ndp1, ndp2, neps1, neps2, nq1, nq2, nwds1, nwds2, n1
parameter (idata = 0, ndp1 = 500, ndp2 = 1000, neps1 = -ndp1, neps2 = -ndp2, &
  nq1 = 11, nq2 = 6 * 2 ** nq1 + 100, nwds1 = int (ndp1 / mpdpw + 2), &
  nwds2 = int (ndp2 / mpdpw + 2))
character*32 chr32
real (8) dplog10q, d1, d2, second, tm0, tm1, tm2
type (mp_real) err, quadgs, fun01, fun02, fun03, fun04, &
  fun15, fun16, fun17, fun18, one, t1, t2, t3, t4, wkgs(-1:nq2), xkgs(-1:nq2), &
  wkgu(-1:nq2), xkgu(-1:nq2), zero
type (mp_real) mppic, mpl02, x1, x2
external fun01, fun02, fun03, fun04, fun15, fun16, fun17, fun18, &
  quadgs, second

!   Check to see if default precision is high enough.

if (ndp2 > mpipl) then
  write (6, '("Increase default precision in module MPFUNF.")')
  stop
endif

!   Compute pi and log(2) to high precision (nwds2 words).

one = mpreal (1.d0, nwds2)
zero = mpreal (0.d0, nwds2)
mppic = mppi (nwds2)
mpl02 = mplog2 (nwds2)

write (6, 1) nq1, ndp1, ndp2, neps1, neps2
1 format ('Quadgs test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,'  Epsilon2 =',i6)

!   Initialize weights and abscissas: wkgs, xkgs, wkgu, xkgu.

write (chr32, '("gauss-"i2.2,"-",i4.4,".dat")') nq1, ndp1
open (11, file = chr32)
rewind 11
tm0 = second ()
if (idata == 0) then

!   Generate xkgs, xkgs, wkgu and xkgu from scratch and write to file.

  call initqgs (nq1, nq2, nwds1, wkgs, xkgs, wkgu, xkgu)
  n1 = dble (xkgs(-1))

  do i = -1, n1
    write (11, '(2i6)') i, n1
    call mpwrite (11, ndp1 + 30, ndp1 + 10, wkgs(i), xkgs(i))
    call mpwrite (11, ndp1 + 30, ndp1 + 10, wkgu(i), xkgu(i))
  enddo
else

!   Read xkgs, xkgs, wkgu and xkgu from file.

  read (11, '(2i6)') i1, i2
  call mpread (11, wkgs(-1), xkgs(-1), nwds1)
  call mpread (11, wkgu(-1), xkgu(-1), nwds1)
  n1 = dble (xkgs(-1))
  if (i1 /= -1 .or. i2 /= n1) stop

  do i = 0, n1
    read (11, '(2i6)') i1, i2
    if (i1 /= i .or. i2 /= n1) stop
    call mpread (11, wkgs(i), xkgs(i), nwds1)
    call mpread (11, wkgu(i), xkgu(i), nwds1)
  enddo
endif
close (11)
tm1 = second ()
tm2 = tm1 - tm0
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.4)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite intervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadgs (fun01, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
3 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (0.25d0, nwds1)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1
4 format ('Actual error =',f10.6,'x10^',i6)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadgs (fun02, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = (mppic - 2.d0 + 2.d0 * mpl02) / 12.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 13)
13 format (/'Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadgs (fun03, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * (exp (0.5d0 * mppic) - 1.d0)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 14)
14 format (/ &
  'Problem 4: Int_0^1 arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2)) dt = 5*Pi^2/96'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadgs (fun04, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 5.d0 * mppic**2 / 96.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 25)
25 format (/ &
   'Functions on the entire real line:'// &
   'Problem 15: Int_-inf^inf 1/(1+t^2) dt = Pi'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun15, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mppic
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 26)
26 format (/'Problem 16: Int_-inf^inf 1/(1+t^4) dt = Pi/Sqrt(2)'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun16, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mppic / sqrt (mpreal (2.d0, nwds1))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 27)
27 format (/'Problem 17: Int_-inf^inf e^(-t^2/2) dt = sqrt (2*Pi)'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun17, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (2.d0 * mppic)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 28)
28 format (/'Problem 18: Int_-inf^inf e^(-t^2/2) cos(t) dt = sqrt (2*Pi/e)'/)
x1 = -one
x2 = one
tm0 = second ()
t1 = quadgs (fun18, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgu, xkgu)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (2.d0 * mppic / exp (mpreal (1.d0, nwds1)))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 99) tm2
99 format ('Total CPU time =',f12.6)

stop
end

function fun01 (t, nwds1, nwds2)

!   fun01(t) = t * log(1+t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun01, t1, t2
type (mp_real) t

t1 = mpreal (t, nwds1)
fun01 = t1 * log (1.d0 + t1)
return
end

function fun02 (t, nwds1, nwds2)

!   fun02(t) = t^2 * arctan(t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun02, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun02 = t1 ** 2 * atan (t1)
return
end

function fun03 (t, nwds1, nwds2)

!   fun03(t) = e^t * cos(t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun03, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun03 = exp(t1) * cos(t1)
return
end

function fun04 (t, nwds1, nwds2)

!   fun04(t) = arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2))

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun04, t1, t2
type (mp_real) t

t1 = mpreal (t, nwds1)
t2 = sqrt (2.d0 + t1**2)
fun04 = atan(t2) / ((1.d0 + t1**2) * t2)
return
end

function fun15 (t, nwds1, nwds2)
use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun15

fun15 = 1.d0 / (1.d0 + t**2)
return
end

function fun16 (t, nwds1, nwds2)
use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun16

fun16 = 1.d0 / (1.d0 + t**4)
return
end

function fun17 (t, nwds1, nwds2)
use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun17

fun17 = exp (-0.5d0 * t**2)
return
end

function fun18 (t, nwds1, nwds2)
use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) t, fun18

fun18 = exp (-0.5d0 * t**2) * cos (t)
return
end

subroutine initqgs (nq1, nq2, nwds1, wkgs, xkgs, wkgu, xkgu)

!   This subroutine initializes the quadrature arrays wkgs and xkgs for standard
!   Gaussian quadrature, and also wkgu and xkgu for quadrature over the real line.
!   It employs a Newton iteration scheme with a dynamic precision level that
!   approximately doubles with each iteration.  The argument nq2, which is the
!   space allocated for wkgs, xkgs, wkgu and xkgu in the calling program, should
!   be at least 6 * 2^nq1 + 100.

!   The wkgu and xkgu arrays are computed from wkgs and xkgs based on the
!   transformation t = tan (pi/2 * x), which transforms an integral on
!   (-infinity, infinity) to an integral on (-1, 1). In particular,
!     xkgu(i) = tan (pi/2 * xkgs(i))
!     wkgu(i) = pi/2 * wkgs(i) / cos (pi/2 * xkgs(i))^2

!   David H Bailey    15 Sep 2018

use mpmodule
implicit none
integer i, ierror, ik0, iprint, j, j1, k, n, ndebug, nq1, nq2, nwds1, nwp
real (8) dpi
parameter (ik0 = 100, iprint = 1, ndebug = 2, dpi = 3.141592653589793238d0)
type (mp_real) eps, pi2, r, t1, t2, t3, t4, t5, wkgs(-1:nq2), xkgs(-1:nq2), &
  wkgu(-1:nq2), xkgu(-1:nq2), zero

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqgs: Gaussian quadrature initialization')
endif

pi2 = 0.5d0 * mppi (nwds1)
zero = mpreal (0.d0, nwds1)
wkgs(-1) = mpreal (dble (nq1), nwds1)
xkgs(-1) = zero
wkgs(0) = zero
xkgs(0) = zero
wkgs(1) = mpreal (dble (nq1), nwds1)
xkgs(1) = mpreal (dble (ik0), nwds1)
wkgu(-1) = mpreal (dble (nq1), nwds1)
xkgu(-1) = zero
wkgu(0) = zero
xkgu(0) = zero
wkgu(1) = mpreal (dble (nq1), nwds1)
xkgu(1) = mpreal (dble (ik0), nwds1)
i = ik0

do j = 2, ik0
  wkgs(j) = zero
  xkgs(j) = zero
  wkgu(j) = zero
  xkgu(j) = zero
enddo

do k = 1, nq1
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, i, nq2
  n = 3 * 2 ** (k + 1)

  do j = 1, n / 2

!   Set working precision = 4 words, and compute a DP estimate of the root.

    nwp = 4
    eps = mpreal (2.d0 ** (48 - nwp * mpnbt), nwds1)
    r = mpreald (cos ((dpi * (j - 0.25d0)) / (n + 0.5d0)), nwp)

!   Compute the j-th root of the n-degree Legendre polynomial using Newton's
!   iteration.

100 continue

!   Perform the next 11 lines with working precision nwp words.

    t1 = mpreal (1.d0, nwp)
    t2 = mpreal (0.d0, nwp)
    r = mpreal (r, nwp)

    do j1 = 1, n
      t3 = t2
      t2 = t1
      t1 = ((2 * j1 - 1) * r * t2 - (j1 - 1) * t3) / j1
    enddo

    t4 = n * (r * t1 - t2) / (r ** 2  - 1.d0)
    t5 = r
    r = r - t1 / t4
    
!   Once convergence is achieved at nwp = 4, then start doubling (almost) the
!   working precision level at each iteration until full precision is reached.

    if (nwp == 4) then
      if (abs (r - t5) > eps) goto 100
      nwp = min (2 * nwp - 1, nwds1)
      goto 100
    elseif (nwp < nwds1) then
      nwp = min (2 * nwp - 1, nwds1)
      goto 100
    endif

    i = i + 1
    if (i > nq2) goto 110
    xkgs(i) = r
    t4 = n * (r * t1 - t2) / (r ** 2  - 1.d0)
    wkgs(i) = 2.d0 / ((1.d0 - r ** 2) * t4 ** 2)
    xkgu(i) = tan (pi2 * r)
    wkgu(i) = pi2 * wkgs(i) / cos (pi2 * r)**2
  enddo

!   Save i (starting index for the next block) in the first 100 elements of 
!   xkgs and xkgu.

  xkgs(k+1) = mpreal (dble (i), nwds1)
  xkgu(k+1) = mpreal (dble (i), nwds1)
enddo

!   Save the size of the arrays in index -1 of xkgs and xkgu.

xkgs(-1) = mpreal (dble (i), nwds1)
xkgu(-1) = mpreal (dble (i), nwds1)
if (ndebug >= 2) then
  write (6, 2) i
2 format ('initqgs: Table spaced used =',i8)
endif
goto 130

110 continue

write (6, 3) nq2
3 format ('initqgs: Table space parameter is too small; value =',i8)
stop

130 continue

return
end

function quadgs (fun, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkgs, xkgs)

!   This routine computes the integral of the function fun on the interval
!   (x1, x2) with a target tolerance of 10^neps1.  The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wkgs and xkgs arrays, which must first be initialized
!   by calling initqgs. The function fun is not evaluated at the endpoints
!   x1 and x2.  If quadgs outputs the message "Increase Quadlevel" or "Terms
!   too large", adjust nq1 and neps2 as necessary in the call to initqgs.

!   Both initqgs and quadgs are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   For some functions, it is important that the endpoints x1 and x2 be
!   computed to high precision (nwds2 words) in the calling program, that
!   these high-precision values be passed to quadgs (where scaled abscissas
!   are calculated to high precision), and that the function definition
!   itself uses these high-precision scaled abscissas in any initial
!   subtractions or other sensitive operations involving the input argument.
!   Otherwise the accuracy of the quadrature result might only be half as
!   high as it otherwise could be.  See the function definitions of fun06,
!   fun07, fun09 and fun10 for examples on how this is done.  Otherwise the
!   function evaluation can and should be performed with low precision
!   (nwds1 words) for faster run times.  The two precision levels (nwds1
!   and nwds2) are specified by the user in the calling program.

!   David H Bailey   15 Sep 2018

use mpmodule
implicit none
integer i, ierror, ik0, k, j, n, ndebug, nds, nq1, nq2, &
  neps1, nerr, nwds1, nwds2, nwp
logical log1, log2
double precision d1, d2, d3, dfrac, dplog10q
parameter (dfrac = 100.d0, ik0 = 100, ndebug = 2)
type (mp_real) ax, bx, c10, eps1, eps2, epsilon1, err, fun, &
  quadgs, tsum, s1, s2, s3, t1, t2, tw1, tw2, twmx, wkgs(-1:nq2), xkgs(-1:nq2), &
  x1, x2, xx1, xx2
external fun, dplog10q

!  These two lines are performed in high precision (nwds2 words).

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)

!  The remaining initialization is performed in low precision (nwds1 words).

epsilon1 = dfrac * mpreal (10.d0, nwds1) ** neps1
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkgs(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadgs: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif

do k = 1, nq1
  n = 3 * 2 ** (k + 1)
  s3 = s2
  s2 = s1
  twmx = mpreal (0.d0, nwds1)
  tsum = mpreal (0.d0, nwds1)
  i = dble (xkgs(k))
  
  do j = 1, n / 2
    i = i + 1

!   These two lines are performed in high precision.

    xx1 = - ax * xkgs(i) + bx
    xx2 = ax * xkgs(i) + bx

    t1 = fun (xx1, nwds1, nwds2)
    tw1 = t1 * wkgs(i)

    if (j + k > 2) then
      t2 = fun (xx2, nwds1, nwds2)
      tw2 = t2 * wkgs(i)
    else
      t2 = mpreal (0.d0, nwds1)
      tw2 = mpreal (0.d0, nwds1)
    endif

    tsum = tsum + tw1 + tw2
    twmx = max (twmx, abs (tw1), abs (tw2))
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.

  s1 =  mpreal (ax, nwds1) * tsum
  eps1 = twmx * epsilon1
  d1 = dplog10q (abs ((s1 - s2) / s1))
  d2 = dplog10q (abs ((s1 - s3) / s1))
  d3 = dplog10q (eps1) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 .eq. -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3)))
  endif
  
!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadgs: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quadgs: Estimated error = 10^',i7)
    goto 140
  endif
enddo

140 continue

quadgs = s1
return
end

function dplog10q (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
real (8) da, dplog10q, t1
type (mp_real) a

call mpmdi (a, da, ia)
if (da .eq. 0.d0) then
  dplog10q = -999999.d0
else
  dplog10q = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

subroutine decmdq (a, b, ib)

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

