program tquad

!   David H Bailey   16 Sep 2018

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2018 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program demonstrates three quadrature routines:

!   quadts:  Implements the tanh-sinh quadrature scheme of Takahashi and Mori,
!     for functions on a finite interval such as (0,1).
!   quades:  Implements the exp-sinh scheme, for functions on a
!     semi-infinite interval such as (0, infinity).
!   quadss:  Implements the sinh-sinh scheme, for functions on the entire
!     real line.

!   The separate program tquadgs.f90 implements Gaussian quadrature.

!   These schemes have some desirable properties, such as the cost of computing
!   weight-abscissa pairs only increases linearly with the number of evaluation
!   points (instead of quadratically, as with Gaussian quadrature). Also,
!   quadts in particular often works well even when the function has a blow-up
!   singularity or infinite derivative at one or both endpoints. The quadrature
!   computations proceed level-by-level, with the computational cost (and,
!   often, the accuracy) approximately doubling with each iteration, until a
!   target accuracy (500-digit precision in this case), based on an accuracy
!   estimate computed by the program, is achieved. The exp-sinh and sinh-sinh
!   schemes are variants of the tanh-sinh scheme, which is described here:

!   David H. Bailey, Xiaoye S. Li and Karthik Jeyabalan, "A comparison of
!   three high-precision quadrature schemes," Experimental Mathematics, 
!   vol. 14 (2005), no. 3, pg. 317-329, preprint available at
!   http://www.davidhbailey.com/dhbpapers/quadrature-em.pdf.

!   The function to be integrated must be defined in an external function
!   subprogram (see samples below), and the name of the function must be
!   included in a "type (mp_real)" and in an "external" statement.  Prior to
!   calling the quadrature routine, the corresponding initialization routine
!   must be called to initialize the abscissa and weight arrays:  initqts,
!   initqes and initqss, corresponding to the quadrature routines quadts,
!   quades and quadss, respectively.

!   All of these routines are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Here are some specific instructions for the individual quadrature schemes:

!   quadts:  

!   First call initqts with the arrays wkts and xkts as the last two arguments,
!   and then call quadts with these same two arrays as the last two arguments. 
!   Set x1 and x2, the limits of integration, in executable statements,
!   using high precision (nwds2 words) if possible. For some integrand
!   functions, it is important that these endpoints be computed to higher
!   precision (nwds2 words) in the calling program, that these higher
!   precision values be passed to quadts (where scaled abscissas are
!   calculated), and that the function definition itself uses these higher
!   precision scaled abscissas in any initial subtractions or other
!   sensitive operations involving the input argument. Otherwise the
!   accuracy of the quadrature result might only be half as high as it
!   otherwise could be. Note in the examples below, x1 and x2 and set
!   using high-precision values of zero, one and pi (the variables zero,
!   one and mppic). See the function definitions of fun06, fun07, fun09
!   and fun10 for examples of how to scale within the function routine. Once
!   the input argument has been accurately scaled, the function evaluation
!   itself can and should be performed with standard precision (nwds1 words)
!   for faster run times.

!   In the initialization routine initqts, weight-abscissa pairs are
!   calculated until the weights are smaller than 10^(neps2), where neps2
!   is the high precision epsilon (typically twice the low precision
!   epsilon -- see below). In some problems it may be necessary to adjust
!   neps2 to a more negative value to obtain the best accuracy.

!   quades:

!   First call initqes with the arrays wkes and xkes as the last two arguments,
!   and then call quades with these same two arrays as the last two arguments. 
!   It is assumed that the left endpoint is x1, and the right endpoint is
!   +infinity. The comment about computing the endpoints to high precision
!   in the note for quadts above also applies here, as does the comment
!   about computing weight-abscissa pairs based on neps2.

!   quadss:

!   First call initqss with the arrays wkss and xkss as the last two arguments,
!   and then call quadss with these same two arrays as the last two arguments. 
!   No endpoints are specified here -- the integral is performed over the
!   entire real line. However, the comment about computing weight-abscissa
!   pairs, based on neps2, also applies here.

!   These inputs are set in the parameter statement below:

!   ndp1   Primary ("low") precision in digits; this is the target accuracy
!          of quadrature results.
!   ndp2   Secondary ("high") precision in digits. By default, ndp2 = 2*ndp1.
!   neps1  Log10 of the primary tolerance. By default, neps1 = - ndp1.
!   neps2  Log10 of the secondary tolerance. By default, neps2 = -ndp2.
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!          (possibly doubles) the number of accurate digits in the result,
!          but also roughly doubles the run time. nq1 must be at least 3.
!   nq2    Space parameter for wk and xk arrays in the calling program.  By
!          default it is set to 12 * 2^nq1. Increase nq2 if directed by a 
!          message produced in initqts.
!   nwds1  Low precision in words. By default nwds1 = int (ndp1 / mpdpw + 2).
!   nwds2  High precision in words. By default nwds2 = int (ndp2 / mpdpw + 2).

use mpmodule
implicit none
integer i, ndp1, ndp2, neps1, neps2, nq1, nq2, nwds1, nwds2, n1
parameter (ndp1 = 500, ndp2 = 1000, neps1 = -ndp1, neps2 = -ndp2, &
  nq1 = 11, nq2 = 12 * 2 ** nq1, nwds1 = int (ndp1 / mpdpw + 2), &
  nwds2 = int (ndp2 / mpdpw + 2))
real (8) dplog10q, d1, d2, second, tm0, tm1, tm2
type (mp_real) err, quades, quadss, quadts, fun01, fun02, fun03, fun04, &
  fun05, fun06, fun07, fun08, fun09, fun10, fun11, fun12, fun13, fun14, &
  fun15, fun16, fun17, fun18, one, t1, t2, t3, t4, wkes(-1:nq2), &
  xkes(-1:nq2), wkss(-1:nq2), xkss(-1:nq2), wkts(-1:nq2), xkts(-1:nq2), zero
type (mp_real) mppic, mpl02, x1, x2
external fun01, fun02, fun03, fun04, fun05, fun06, fun07, fun08, &
  fun09, fun10, fun11, fun12, fun13, fun14, fun15, fun16, fun17, fun18, &
  quades, quadss, quadts, second

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
1 format ('Quadts test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,'  Epsilon2 =',i6)

!   Initialize quadrature tables wk and xk (weights and abscissas).

tm0 = second ()
call initqes (nq1, nq2, nwds1, neps2, wkes, xkes)
call initqss (nq1, nq2, nwds1, neps2, wkss, xkss)
call initqts (nq1, nq2, nwds1, neps2, wkts, xkts)
tm1 = second ()
tm2 = tm1 - tm0
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite intervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun01, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
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
t1 = quadts (fun02, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
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
t1 = quadts (fun03, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
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
t1 = quadts (fun04, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 5.d0 * mppic**2 / 96.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 15)
15 format (/&
  'Continuous functions on finite intervals, but non-diff at an endpoint:'// &
  'Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun05, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (-4.d0, nwds1) / 9.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 16)
16 format (/'Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun06, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.25d0 * mppic
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 17)
17 format (/&
  'Functions on finite intervals with integrable singularity at an endpoint:'//&
  'Problem 7: Int_0^1 sqrt(t)/sqrt(1-t^2) dt = 2*sqrt(pi)*gamma(3/4)/gamma(1/4)'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun07, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 2.d0 * sqrt (mpreal (mppic, nwds1)) * gamma (mpreal (0.75d0, nwds1)) &
  / gamma (mpreal (0.25d0, nwds1))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 18)
18 format (/'Problem 8: Int_0^1 log(t)^2 dt = 2'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun08, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (2.d0, nwds1)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 19)
19 format (/'Problem 9: Int_0^(pi/2) log(cos(t)) dt = -pi*log(2)/2'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun09, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = -0.5d0 * mppic * mpl02
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 20)
20 format (/'Problem 10: Int_0^(pi/2) sqrt(tan(t)) dt = pi*sqrt(2)/2'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun10, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * mppic * sqrt (mpreal (2.d0, nwds1))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 21)
21 format (/&
  'Functions on a semi-infinite interval:'//&
  'Problem 11: Int_1^inf 1/(1+t^2) dt = pi/4'/)
x1 = one
tm0 = second ()
t1 = quades (fun11, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.25d0 * mppic
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 22)
22 format (/'Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)'/)
x1 = zero
tm0 = second ()
t1 = quades (fun12, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (mppic)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 23)
23 format (/'Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)'/)
x1 = zero
tm0 = second ()
t1 = quades (fun13, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (0.5d0 * mppic)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 24)
24 format (/'Problem 14: Int_pi^inf e^(-t)*cos(t) dt = -1/2 * exp(-pi)'/)
x1 = mppic
tm0 = second ()
t1 = quades (fun14, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = -0.5d0 * exp (- mpreal (mppic, nwds1))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 25)
25 format (/ &
   'Functions on the entire real line:'// &
   'Problem 15: Int_-inf^inf 1/(1+t^2) dt = Pi'/)
tm0 = second ()
t1 = quadss (fun15, nq1, nq2, nwds1, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mppic
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 26)
26 format (/'Problem 16: Int_-inf^inf 1/(1+t^4) dt = Pi/Sqrt(2)'/)
tm0 = second ()
t1 = quadss (fun16, nq1, nq2, nwds1, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mppic / sqrt (mpreal (2.d0, nwds1))
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 27)
27 format (/'Problem 17: Int_-inf^inf e^(-t^2/2) dt = sqrt (2*Pi)'/)
tm0 = second ()
t1 = quadss (fun17, nq1, nq2, nwds1, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (2.d0 * mppic)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 28)
28 format (/'Problem 18: Int_-inf^inf e^(-t^2/2) cos(t) dt = sqrt (2*Pi/e)'/)
tm0 = second ()
t1 = quadss (fun18, nq1, nq2, nwds1, neps1, wkss, xkss)
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

function fun05 (t, nwds1, nwds2)

!    fun05(t) = sqrt(t)*log(t)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun05, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun05 = sqrt (t1) * log (t1)
return
end

function fun06 (t, nwds1, nwds2)

!    fun06(t) = sqrt(1-t^2)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun06, t1, t2
type (mp_real) t

!   The subtraction to compute t2 must be performed using high precision
!   (nwds2), but after the subtraction its low precision value is fine.

t1 = mpreal (t, nwds1)
t2 = mpreal (1.d0 - t**2, nwds1)
fun06 = sqrt (t2)
return
end

function fun07 (t, nwds1, nwds2)

!   fun07(t) = sqrt (t) / sqrt(1-t^2)

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun07, t1, t2
type (mp_real) t

!   The subtraction to compute t2 must be performed using high precision
!   (nwds2), but after the subtraction its low precision value is fine.

t1 = mpreal (t, nwds1)
t2 = mpreal (1.d0 - t, nwds1)
fun07 = sqrt (t1) / sqrt (t2 * (1.d0 + t1))
return
end

function fun08 (t, nwds1, nwds2)

!   fun08(t) = log(t)^2

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun08, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun08 = log (t1) ** 2
return
end

function fun09 (t, nwds1, nwds2)

!   fun09(t) = log (cos (t))

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun09, pi, t1, t2, t3, t4
type (mp_real) t

!   The subtraction to compute t2 must be performed using high precision
!   (nwds2), but after the subtraction its low precision value is fine.

t1 = mpreal (t, nwds1)
pi = mppi (nwds2)
t3 = mpreal (0.25d0 * pi, nwds1)
t2 = mpreal (0.5d0 * pi - t, nwds1)

if (t1 < t3) then
  t4 = cos (t1)
else
  t4 = sin (t2)
endif
fun09 = log (t4)

return
end

function fun10 (t, nwds1, nwds2)

!   fun10(t) = sqrt(tan(t))

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun10, pi, t1, t2, t3, t4
type (mp_real) t

!   The subtraction to compute t3 must be performed using high precision
!   (nwds2), but after the subtraction its low precision value is fine.

t1 = mpreal (t, nwds1)
pi = mppi (nwds2)
t3 = mpreal (0.25d0 * pi, nwds1)
t2 = mpreal (0.5d0 * pi - t, nwds1)

if (t1 < t3) then
  fun10 = sqrt (tan (t1))
else
  fun10 = 1.d0 / sqrt (tan (t2))
endif
return
end

function fun11 (t, nwds1, nwds2)

!   1/(1 + t^2) on (0, Inf).

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun11, t1
type (mp_real) t

t1 = mpreal (t, nwds1)
fun11 = 1.d0 / (1.d0 + t1 ** 2)
return
end

function fun12 (t, nwds1, nwds2)

!   e^(-t)/sqrt(t) on (0, inf).

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun12, t1, t2
type (mp_real) t

t1 = mpreal (t, nwds1)
fun12 = exp (-t1) / sqrt (t1)
return
end

function fun13 (t, nwds1, nwds2)

!   e^(-t^2/2) on (0, inf).

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun13, t1, t2
type (mp_real) t

t1 = mpreal (t, nwds1)
fun13 = exp (-0.5d0 * t1 ** 2)
return
end

function fun14 (t, nwds1, nwds2)

!   e^(-t) cos(t) on (0, inf).

use mpmodule
implicit none
integer nwds1, nwds2
type (mp_real) fun14, t1, t2
type (mp_real) t

t1 = mpreal (t, nwds1)
fun14 = exp (-t1) * cos (t1)
return
end

function fun15 (t, nwds1)
use mpmodule
implicit none
integer nwds1
type (mp_real) t, fun15

fun15 = 1.d0 / (1.d0 + t**2)
return
end

function fun16 (t, nwds1)
use mpmodule
implicit none
integer nwds1
type (mp_real) t, fun16

fun16 = 1.d0 / (1.d0 + t**4)
return
end

function fun17 (t, nwds1)
use mpmodule
implicit none
integer nwds1
type (mp_real) t, fun17

fun17 = exp (-0.5d0 * t**2)
return
end

function fun18 (t, nwds1)
use mpmodule
implicit none
integer nwds1
type (mp_real) t, fun18

fun18 = exp (-0.5d0 * t**2) * cos (t)
return
end

subroutine initqes (nq1, nq2, nwds1, neps2, wkes, xkes)

!   This subroutine initializes the quadrature arrays xkes and wkes for quades.
!   The argument nq2 is the space allocated for wkes and xkes in the calling
!   program.  By default it is set to 12 * 2^nq1.  If initqes outputs the message
!   "Table space parameter is too small", adjust nq2.  Also, if quades outputs
!   the message "Terms too large", adjust nq1 and neps2 as necessary in the call
!   to initqes. The argument neps2 controls termination of the loop below, which
!   ends when wkes(k) * 10^(neps2) > 1.

!   The wkes and xkes arrays are computed based on the transformation
!   t = exp (sinh (pi/2 * x).  See comments below.

!   Both initqes and quades are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey   16 Sep 2018

use mpmodule
implicit none
integer i, ierror, iprint, j, k, k1, ndebug, neps2, nq1, nq2, nwds1
parameter (iprint = 1024, ndebug = 2)
type (mp_real) eps2, h, p2, t1, t2, t3, t4, u1, u2, &
  wkes(-1:nq2), xkes(-1:nq2)

write (6, 1)
1 format ('initqes: Exp-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwds1) ** neps2
p2 = 0.5d0 * mppi (nwds1)
h = mpreal (0.5d0 ** nq1, nwds1)
wkes(-1) = mpreal (dble (nq1), nwds1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
    t1 = mpreal (dble (k) * h, nwds1)

!   xkes(k) = exp (u1)
!   wkes(k) = exp (u1) * u2
!   where u1 = pi/2 * sinh (t1) and u2 = pi/2 * cosh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  xkes(k) = exp (u1)
  wkes(k) = xkes(k) * u2

  if (wkes(k) * eps2 > 1.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqes: Table space parameter is too small; value =',i8)
stop

100 continue

xkes(-1) = mpreal (dble (k), nwds1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqes: Table spaced used =',i8)
endif

return
end

function quades (fun, x1, nq1, nq2, nwds1, nwds2, neps1, wkes, xkes)

!   This routine computes the integral of the function fun on the interval
!   (x1, inf) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wkes and xkes arrays, which must first be 
!   initialized by calling initqes. If quades outputs the message "Terms too
!   large", adjust nq1 and neps2 as necessary in the call to initqes.

!   Both initqes and quades are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   For some functions, it is important that the endpoint x1 be
!   computed to high precision (nwds2 words) in the calling program, that
!   these high-precision values be passed to quadges (where scaled abscissas
!   are calculated to high precision), and that the function definition
!   itself uses these high-precision scaled abscissas in any initial
!   subtractions or other sensitive operations involving the input argument.
!   Otherwise the accuracy of the quadrature result might only be half as
!   high as it otherwise could be. See the function definitions of fun06,
!   fun07, fun09 and fun10 for examples on how this is done.  Otherwise the
!   function evaluation can and should be performed with low precision
!   (nwds1 words) for faster run times. The two precision levels (nwds1
!   and nwds2) are specified by the user in the calling program.

!   David H Bailey  16 Sep 2018

use mpmodule
implicit none
integer i, ierror, ip(0:100), iz1, iz2, izx, j, k, k1, k2, n, ndebug, &
  nds, neps1, nq1, nq2, nqq1, nwds1, nwds2
parameter (izx = 5, ndebug = 2)
logical log1
real (8) d1, d2, d3, d4, dplog10q
type (mp_real) c10, eps1, eps2, epsilon1, err, fun, h, &
  quades, tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx, &
  wkes(-1:nq2), xkes(-1:nq2)
type (mp_real) x1, xx1, xx2
external fun, dplog10q

epsilon1 = mpreal (10.d0, nwds1) ** neps1
tsum = mpreal (0.d0, nwds1)
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
h = mpreal (1.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkes(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quades: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wkes(-1))
n = dble (xkes(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = mpreal (0.d0, nwds1)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwds2 words) to ensure full accuracy in the quadrature
!   results, even though the abscissas xkes(i) were computed in low precision.

      xx1 = x1 + mpreal (xkes(i), nwds2)
      xx2 = x1 + 1.d0 / mpreal (xkes(i), nwds2)
      log1 = xx1 > x1
  
!   The remaining computations are performed in low precision (nwds1 words).

      if (iz1 < izx) then
        t1 = fun (xx1, nwds1)
        tw1 = t1 * wkes(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = mpreal (0.d0, nwds1)
        tw1 = mpreal (0.d0, nwds1)
      endif

      if (i > 0 .and. log1 .and. iz2 < izx) then
        t2 = fun (xx2, nwds1)
        tw2 = t2 * wkes(i) / xkes(i)**2
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = mpreal (0.d0, nwds1)
        tw2 = mpreal (0.d0, nwds1)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  h * tsum
  eps1 = twmx * epsilon1
  eps2 = abs (max (twi1, twi2) / s1)
  d1 = dplog10q (abs ((s1 - s2) / s1))
  d2 = dplog10q (abs ((s1 - s3) / s1))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 .eq. -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quades: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quades: Terms too large -- adjust neps2 in call to initqes.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quades: Estimated error = 10^',i7)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10q (abs (err)))
5   format ('quades: Estimated error = 10^',i7/&
    'Adjust nq1 and neps2 in initqes for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quades = s1
return
end

subroutine initqss (nq1, nq2, nwds1, neps2, wkss, xkss)

!   This subroutine initializes the quadrature arrays xkss and wkss for quadss.
!   The argument nq2 is the space allocated for wkss and xkss in the calling
!   program. By default it is set to 12 * 2^nq1. If initqss outputs the message
!   "Table space parameter is too small", adjust nq2. Also, if quadss outputs
!   the message "Terms too large", adjust nq1 and neps2 as necessary in the call
!   to initqss. The argument neps2 controls termination of the loop below, which
!   ends when wkss(k) * 10^(neps2) > 1.

!   The wkss and xkss arrays are computed based on the transformation
!   t = sinh (sinh (pi/2 * x).  See comments below.

!   Both initqss and quadss are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey   16 Sep 2018

use mpmodule
implicit none
integer i, ierror, iprint, j, k, k1, ndebug, neps2, nq1, nq2, nwds1
parameter (iprint = 1024, ndebug = 2)
type (mp_real) eps2, h, p2, t1, t2, t3, t4, u1, u2, &
  wkss(-1:nq2), xkss(-1:nq2)

write (6, 1)
1 format ('initqss: Sinh-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwds1) ** neps2
p2 = 0.5d0 * mppi (nwds1)
h = mpreal (0.5d0 ** nq1, nwds1)
wkss(-1) = mpreal (dble (nq1), nwds1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
    t1 = mpreal (dble (k) * h, nwds1)

!   xkss(k) = sinh (u1)
!   wkss(k) = cosh (u1) * u2
!   where u1 = pi/2 * sinh (t1) and u2 = pi/2 * cosh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  t3 = exp (u1)
  xkss(k) = 0.5d0 * (t3 - 1.d0 / t3)
  wkss(k) = 0.5d0 * (t3 + 1.d0 / t3) * u2

  if (wkss(k) * eps2 > 1.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqss: Table space parameter is too small; value =',i8)
stop

100 continue

xkss(-1) = mpreal (dble (k), nwds1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqss: Table spaced used =',i8)
endif

return
end

function quadss (fun, nq1, nq2, nwds1, neps1, wkss, xkss)

!   This routine computes the integral of the function fun on the interval
!   (-inf, inf) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wkss and xkss arrays, which must first be 
!   initialized by calling initqss. If quadss outputs the message "Terms too
!   large", adjust nq1 and neps2 as necessary in the call to initqss.

!   Both initqss and quadss are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey  16 Sep 2018

use mpmodule
implicit none
integer i, ierror, ip(0:100), iz1, iz2, izx, j, k, k1, k2, n, ndebug, &
  nds, neps1, nq1, nq2, nqq1, nwds1
parameter (izx = 5, ndebug = 2)
real (8) d1, d2, d3, d4, dplog10q
type (mp_real) c10, eps1, eps2, epsilon1, err, fun, h, &
  quadss, tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx, &
  wkss(-1:nq2), xkss(-1:nq2)
external fun, dplog10q

epsilon1 = mpreal (10.d0, nwds1) ** neps1
tsum = mpreal (0.d0, nwds1)
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
h = mpreal (1.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkss(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadss: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wkss(-1))
n = dble (xkss(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = mpreal (0.d0, nwds1)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then
      if (iz1 < izx) then
        t1 = fun (xkss(i), nwds1)
        tw1 = t1 * wkss(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = mpreal (0.d0, nwds1)
        tw1 = mpreal (0.d0, nwds1)
      endif

      if (i > 0 .and. iz2 < izx) then
        t2 = fun (-xkss(i), nwds1)
        tw2 = t2 * wkss(i)
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = mpreal (0.d0, nwds1)
        tw2 = mpreal (0.d0, nwds1)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  h * tsum
  eps1 = twmx * epsilon1
  eps2 = abs (max (twi1, twi2) / s1)
  d1 = dplog10q (abs ((s1 - s2) / s1))
  d2 = dplog10q (abs ((s1 - s3) / s1))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 .eq. -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadss: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quadss: Terms too large -- adjust neps2 in call to initqss.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quadss: Estimated error = 10^',i7)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10q (abs (err)))
5   format ('quadss: Estimated error = 10^',i7/&
    'Adjust nq1 and neps2 in initqss for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quadss = s1
return
end

subroutine initqts (nq1, nq2, nwds1, neps2, wkts, xkts)

!   This subroutine initializes the quadrature arrays xkts and wkts for quadts.
!   The argument nq2 is the space allocated for wkts and xkts in the calling
!   program.  By default it is set to 12 * 2^nq1.  If initqts outputs the message
!   "Table space parameter is too small", adjust nq2.  Also, if quadts outputs
!   the message "Terms too large", adjust nq1 and neps2 as necessary in the call
!   to initqts. The argument neps2 controls termination of the loop below, which
!   ends when wkts(k) < 10^(neps2).

!   The wkts and xkts arrays are computed based on the transformation
!   t = tanh (sinh (pi/2 * x).  Note however that xkts contains one minus the
!   conventional abscissas, in order to conserve precision. See comments below.

!   Both initqts and quadts are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey   16 Sep 2018

use mpmodule
implicit none
integer i, ierror, iprint, j, k, k1, ndebug, neps2, nq1, nq2, nwds1
parameter (iprint = 1024, ndebug = 2)
type (mp_real) eps2, h, p2, t1, t2, t3, t4, u1, u2, wkts(-1:nq2), xkts(-1:nq2)

write (6, 1)
1 format ('initqts: Tanh-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwds1) ** neps2
p2 = 0.5d0 * mppi (nwds1)
h = mpreal (0.5d0 ** nq1, nwds1)
wkts(-1) = mpreal (dble (nq1), nwds1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = mpreal (dble (k) * h, nwds1)

!   xkts(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wkts(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  t3 = exp (u2)
  t4 = 0.5d0 * (t3 + 1.d0 / t3)
  xkts(k) = 1.d0 / (t3 * t4)
  wkts(k) = u1 / t4 ** 2

  if (wkts(k) < eps2) goto 100
enddo

write (6, 2) nq2
2 format ('initqts: Table space parameter is too small; value =',i8)
stop

100 continue

xkts(-1) = mpreal (dble (k), nwds1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
endif

return
end

function quadts (fun, x1, x2, nq1, nq2, nwds1, nwds2, neps1, wkts, xkts)

!   This routine computes the integral of the function fun on the interval
!   (x1, x2) with a target tolerance of 10^neps1. The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wkts and xkts arrays, which must first be 
!   initialized by calling initqts. If quadts outputs the message "Terms too
!   large", adjust nq1 and neps2 as necessary in the call to initqts.

!   Both initqts and quadts are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   For some functions, it is important that the endpoint x1 be
!   computed to high precision (nwds2 words) in the calling program, that
!   these high-precision values be passed to quadts (where scaled abscissas
!   are calculated to high precision), and that the function definition
!   itself uses these high-precision scaled abscissas in any initial
!   subtractions or other sensitive operations involving the input argument.
!   Otherwise the accuracy of the quadrature result might only be half as
!   high as it otherwise could be. See the function definitions of fun06,
!   fun07, fun09 and fun10 for examples on how this is done.  Otherwise the
!   function evaluation can and should be performed with low precision
!   (nwds1 words) for faster run times. The two precision levels (nwds1
!   and nwds2) are specified by the user in the calling program.

!   David H Bailey  16 Sep 2018

use mpmodule
implicit none
integer i, ierror, ip(0:100), iz1, iz2, izx, j, k, k1, k2, n, ndebug, &
  nds, neps1, nq1, nq2, nqq1, nwds1, nwds2
parameter (izx = 5, ndebug = 2)
logical log1, log2
real (8) d1, d2, d3, d4, dplog10q
type (mp_real) c10, eps1, eps2, epsilon1, err, fun, h, &
  quadts, tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx, &
  wkts(-1:nq2), xkts(-1:nq2)
type (mp_real) ax, bx, x1, x2, xki, xt1, xx1, xx2
external fun, dplog10q

!  These two lines are performed in high precision (nwds2 words).

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)

!  The remaining initialization is performed in low precision (nwds1 words).

epsilon1 = mpreal (10.d0, nwds1) ** neps1
tsum = mpreal (0.d0, nwds1)
s1 = mpreal (0.d0, nwds1)
s2 = mpreal (0.d0, nwds1)
h = mpreal (1.d0, nwds1)
c10 = mpreal (10.d0, nwds1)

if (wkts(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadts: Quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wkts(-1))
n = dble (xkts(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = mpreal (0.d0, nwds1)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1  
    if (mod (i, k2) /= 0 .or. k == 1) then

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwds2 words) to ensure full accuracy in the quadrature
!   results, even though the abscissas xkts(i) were computed in low precision.

      xki = xkts(i)
      xt1 = 1.d0 - mpreal (xki, nwds2)
      xx1 = - ax * xt1 + bx
      xx2 = ax * xt1 + bx
      log1 = xx1 > x1
      log2 = xx2 < x2      

!   The remaining computations are performed in low precision (nwds1 words).

      if (log1 .and. iz1 < izx) then
        t1 = fun (xx1, nwds1, nwds2)
        tw1 = t1 * wkts(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = mpreal (0.d0, nwds1)
        tw1 = mpreal (0.d0, nwds1)
      endif

      if (i > 0 .and. log2 .and. iz2 < izx) then
        t2 = fun (xx2, nwds1, nwds2)
        tw2 = t2 * wkts(i)
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = mpreal (0.d0, nwds1)
        tw2 = mpreal (0.d0, nwds1)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  mpreal (ax, nwds1) * h * tsum
  eps1 = twmx * epsilon1
  eps2 = abs (max (twi1, twi2) / s1)
  d1 = dplog10q (abs ((s1 - s2) / s1))
  d2 = dplog10q (abs ((s1 - s3) / s1))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = mpreal (1.d0, nwds1)
  elseif (d1 .eq. -999999.d0) then
    err = mpreal (0.d0, nwds1)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadts: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call mpwrite (6, 80, 60, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quadts: Terms too large -- adjust neps2 in call to initqts.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quadts: Estimated error = 10^',i7)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10q (abs (err)))
5   format ('quadts: Estimated error = 10^',i7/&
    'Adjust nq1 and neps2 in initqts for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quadts = s1
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
