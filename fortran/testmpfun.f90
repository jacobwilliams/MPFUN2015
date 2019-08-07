program testmpfun

!   This briefly tests most individual MPFUN operations and functions
!   (including all mixed mode arithmetic and comparison operations), by
!   comparing each result with the equivalent operation performed in double
!   precision or double complex as appropriate.  This is clearly not an
!   exhaustive or fully high-precision test, but it often detects software bugs
!   and compiler issues.  All multiprecision results are output to 500 digits,
!   so that the output file can be compared with a reference file.

!   David H Bailey  19 Apr 2017

use mpmodule
implicit none
integer i, ndp, nwds
parameter (ndp = 500, nwds = ndp / mpdpw + 2)
character(500) chr500
character(32) chr32
character(1) chr1(ndp+30), chrpi(ndp+30)
integer i1, i2, i3, i4
logical l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
real (8) d1, d2, d3, d4, e1, e2
complex (8) dc1, dc2, dc3, dc4, ec1, ec2
type (mp_real) t1, t2, t3, t4
type (mp_complex) z1, z2, z3, z4

write (6, '(a)') 'MPUN2015 quick check of operations and functions'

!   Define a few sample data values.

chr500 = &
'3.14159265358979323846264338327950288419716939937510582097494459230781&
&6406286208998628034825342117067982148086513282306647093844609550582231&
&7253594081284811174502841027019385211055596446229489549303819644288109&
&7566593344612847564823378678316527120190914564856692346034861045432664&
&8213393607260249141273724587006606315588174881520920962829254091715364&
&3678925903600113305305488204665213841469519415116094330572703657595919&
&5309218611738193261179310511854807446237996274956735188575272489122793&
&81830119491'

do i = 1, ndp
  chrpi(i) = chr500(i:i)
enddo

t1 = mppi (nwds)
t2 = - mplog2 (nwds)
d1 = t1
d2 = t2
e1 = 3141.d0 / 8192.d0
e2 = 6931.d0 / 8192.d0
z1 = mpcmplx (0.5d0 * mppi (nwds), exp (mpreal (0.5d0, nwds)), nwds)
z2 = mpcmplx (- gamma (mpreal (0.5d0, nwds)), cos (mpreal (1.d0, nwds)), nwds)
dc1 = z1
dc2 = z2
ec1 = cmplx (e1, e2, mprknd)
ec2 = cmplx (-e2, e1, mprknd)
i1 = 5
i2 = -3

write (6, '(/a/)') 'Data arrays:'
write (6, '(a)') 't1 = pi:'
call mpwrite (6, ndp + 20, ndp, t1)
write (6, '(a)') 't2 = -log(2):'
call mpwrite (6, ndp + 20, ndp, t2)
write (6, '(a)') 'z1 = (0.5*pi, exp(0.5)):'
call mpwrite (6, ndp + 20, ndp, z1)
write (6, '(a)') 'z2 = (-Gamma(0.5), Cos(1)):'
call mpwrite (6, ndp + 20, ndp, z2)

write (6, '(a)') 'e1 = 3141/8192:'
write (6, '(1p,d25.15)') e1
write (6, '(a)') 'e2 = 6931 / 8192:'
write (6, '(1p,d25.15)') e2
write (6, '(a)')  'ec1 = (3141/8192, 6931/8192)'
write (6, '(1p,2d25.15)') ec1
write (6, '(a)') 'ec2 = (6931/8192, 3141/8192):'
write (6, '(1p,2d25.15)') ec2

write (6, '(/a/)') 'Real data operations:'

write (6, '(a)') 'addition: t1+t2 ='
call mpwrite (6, ndp + 20, ndp, t1 + t2)
write (6, 1) d1 + d2
1 format (1p,2d25.15)
call checkdp (t1 + t2, d1 + d2, ndp)

write (6, '(a)') 'addition: t1+e2 ='
call mpwrite (6, ndp + 20, ndp, t1 + e2)
write (6, 1) d1 + e2
call checkdp (t1 + e2, d1 + e2, ndp)

write (6, '(a)') 'addition: e1+t2 ='
call mpwrite (6, ndp + 20, ndp, e1 + t2)
write (6, 1) e1 + d2
call checkdp (e1 + t2, e1 + d2, ndp)

write (6, '(a)') 'subtraction: t1-t2 ='
call mpwrite (6, ndp + 20, ndp, t1 - t2)
write (6, 1) d1 - d2
call checkdp (t1 - t2, d1 - d2, ndp)

write (6, '(a)') 'subtraction: t1-e2 ='
call mpwrite (6, ndp + 20, ndp, t1 - e2)
write (6, 1) d1 - e2
call checkdp (t1 - e2, d1 - e2, ndp)

write (6, '(a)') 'subtraction: e1-t2 ='
call mpwrite (6, ndp + 20, ndp, e1 - t2)
write (6, 1) e1 - d2
call checkdp (e1 - t2, e1 - d2, ndp)

write (6, '(a)') 'multiplication: t1*t2 ='
call mpwrite (6, ndp + 20, ndp, t1 * t2)
write (6, 1) d1 * d2
call checkdp (t1 * t2, d1 * d2, ndp)

write (6, '(a)') 'multiplication: t1*e2 ='
call mpwrite (6, ndp + 20, ndp, t1 * e2)
write (6, 1) d1 * e2
call checkdp (t1 * e2, d1 * e2, ndp)

write (6, '(a)') 'multiplication: e1*t2 ='
call mpwrite (6, ndp + 20, ndp, e1 * t2)
write (6, 1) e1 * d2
call checkdp (e1 * t2, e1 * d2, ndp)

write (6, '(a)') 'division: t1/t2 ='
call mpwrite (6, ndp + 20, ndp, t1 / t2)
write (6, 1) d1 / d2
call checkdp (t1 / t2, d1 / d2, ndp)

write (6, '(a)') 'division: t1/e2 ='
call mpwrite (6, ndp + 20, ndp, t1 / e2)
write (6, 1) d1 / e2
call checkdp (t1 / e2, d1 / e2, ndp)

write (6, '(a)') 'division: e1/t2 ='
call mpwrite (6, ndp + 20, ndp, e1 / t2)
write (6, 1) e1 / d2
call checkdp (e1 / t2, e1 / d2, ndp)

write (6, '(a)') 'exponentiation: t1**i1 ='
call mpwrite (6, ndp + 20, ndp, t1 ** i1)
write (6, 1) d1 ** i1
call checkdp (t1 ** i1, d1 ** i1, ndp)

write (6, '(a)') 'exponentiation: t1**t2 ='
call mpwrite (6, ndp + 20, ndp, t1 ** t2)
write (6, 1) d1 ** d2
call checkdp (t1 ** t2, d1 ** d2, ndp)

write (6, '(a)') 'equal test: t1 == t2, e1 == t2, t1 == e2'
l1 = t1 == t2; l2 = d1 == d2; l3 = e1 == t2; l4 = e1 == d2; l5 = t1 == e2; l6 = d1 == e2
write (6, '(6l4)') l1, l2, l3, l4, l5, l6
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)

write (6, '(a)') 'not-equal test: t1 /= t2, e1 /= t2, t1 =/ e2'
l1 = t1 /= t2; l2 = d1 /= d2; l3 = e1 /= t2; l4 = e1 /= d2; l5 = t1 /= e2; l6 = d1 /= e2
write (6, '(6l4)') l1, l2, l3, l4, l5, l6
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)

write (6, '(a)') 'less-than-or-equal test: t1 <= t2, e1 <= t2, t1 <= e2'
l1 = t1 <= t2; l2 = d1 <= d2; l3 = e1 <= t2; l4 = e1 <= d2; l5 = t1 <= e2; l6 = d1 <= e2
write (6, '(6l4)') l1, l2, l3, l4, l5, l6
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)

write (6, '(a)') 'greater-than-or-equal test: t1 >= t2, e1 >= t2, t1 >= e2'
l1 = t1 >= t2; l2 = d1 >= d2; l3 = e1 >= t2; l4 = e1 >= d2; l5 = t1 >= e2; l6 = d1 >= e2
write (6, '(6l4)') l1, l2, l3, l4, l5, l6
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)

write (6, '(a)') 'less-than test: t1 < t2, e1 < t2, t1 < e2'
l1 = t1 < t2; l2 = d1 < d2; l3 = e1 < t2; l4 = e1 < d2; l5 = t1 < e2; l6 = d1 < e2
write (6, '(6l4)') l1, l2, l3, l4, l5, l6
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)

write (6, '(a)') 'greater-than test: t1 > t2, e1 > t2, t1 > e2'
l1 = t1 > t2; l2 = d1 > d2; l3 = e1 > t2; l4 = e1 > d2; l5 = t1 > e2; l6 = d1 > e2
write (6, '(6l4)') l1, l2, l3, l4, l5, l6
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)

write (6, '(a)') 'abs(t2) ='
call mpwrite (6, ndp + 20, ndp, abs (t2))
write (6, '(1p,2d25.15)') abs (d2)
call checkdp (abs (t2), abs (d2), ndp)

write (6, '(a)') 'acos(t2) ='
call mpwrite (6, ndp + 20, ndp, acos (t2))
write (6, '(1p,2d25.15)') acos (d2)
call checkdp (acos (t2), acos (d2), ndp)

write (6, '(a)') 'aint(t1) ='
call mpwrite (6, ndp + 20, ndp, aint (t1))
write (6, '(1p,2d25.15)') aint (d1)
call checkdp (aint (t1), aint (d1), ndp)

write (6, '(a)') 'anint(t1) ='
call mpwrite (6, ndp + 20, ndp, anint (t1))
write (6, '(1p,2d25.15)') anint (d1)
call checkdp (anint (t1), anint (d1), ndp)

write (6, '(a)') 'asin(t2) ='
call mpwrite (6, ndp + 20, ndp, asin (t2))
write (6, '(1p,2d25.15)') asin (d2)
call checkdp (asin (t2), asin (d2), ndp)

write (6, '(a)') 'atan(t1) ='
call mpwrite (6, ndp + 20, ndp, atan (t1))
write (6, '(1p,2d25.15)') atan (d1)
call checkdp (atan (t1) , atan (d1), ndp)

write (6, '(a)') 'atan2(t1,t2) ='
call mpwrite (6, ndp + 20, ndp, atan2 (t1,t2))
write (6, '(1p,2d25.15)') atan2 (d1,d2)
call checkdp (atan2 (t1, t2), atan2 (d1, d2), ndp)

write (6, '(a)') 'cos(t2) ='
call mpwrite (6, ndp + 20, ndp, cos (t2))
write (6, '(1p,2d25.15)') cos (d2)
call checkdp (cos (t2), cos (d2), ndp)

write (6, '(a)') 'cosh(t1) ='
call mpwrite (6, ndp + 20, ndp, cosh (t1))
write (6, '(1p,2d25.15)') cosh (d1)
call checkdp (cosh (t1), cosh (d1), ndp)

write (6, '(a)') 'exp(t1) ='
call mpwrite (6, ndp + 20, ndp, exp (t1))
write (6, '(1p,2d25.15)') exp (d1)
call checkdp (exp (t1), exp (d1), ndp)

write (6, '(a)') 'gamma(t1) ='
call mpwrite (6, ndp + 20, ndp, gamma (t1))
write (6, '(1p,2d25.15)') gamma (d1)
call checkdp (gamma (t1), gamma (d1), ndp)

write (6, '(a)') 'log(t1) ='
call mpwrite (6, ndp + 20, ndp, log (t1))
write (6, '(1p,2d25.15)') log (d1)
call checkdp (log (t1), log (d1), ndp)

write (6, '(a)') 'max(t1,t2) ='
call mpwrite (6, ndp + 20, ndp, max (t1,t2))
write (6, '(1p,2d25.15)') max (d1,d2)
call checkdp (max (t1, t2), max (d1, d2), ndp)

write (6, '(a)') 'min(t1,t2) ='
call mpwrite (6, ndp + 20, ndp, min (t1,t2))
write (6, '(1p,2d25.15)') min (d1,d2)
call checkdp (min (t1, t2), min (d1, d2), ndp)

write (6, '(a)') 'mpcssh(t1) ='
call mpcssh (t1, t3, t4)
call mpwrite (6, ndp + 20, ndp, t3, t4)
d3 = cosh (d1)
d4 = sinh (d1)
write (6, '(1p,2d25.15)') d3, d4
call checkdp (t3, d3, ndp)
call checkdp (t4, d4, ndp)

write (6, '(a)') 'mpcssn(t2) ='
call mpcssn (t2, t3, t4)
call mpwrite (6, ndp + 20, ndp, t3, t4)
d3 = cos (d2)
d4 = sin (d2)
write (6, '(1p,2d25.15)') d3, d4
call checkdp (t3, d3, ndp)
call checkdp (t4, d4, ndp)

write (6, '(a)') 'mpeform(t1) ='
call mpeform (t1, ndp + 20, ndp, chr1)
write (6, '(80a1)') (chr1(i), i = 1, ndp + 20)
write (6, '(1p,d25.15)') d1
call checkdp (mpreal (chr1, ndp, nwds), d1, nwds)

write (6, '(a)') 'mpfform(t1) ='
call mpfform (t1, ndp + 20, ndp, chr1)
write (6, '(80a1)') (chr1(i), i = 1, ndp + 20)
write (6, '(f25.15)') d1
call checkdp (mpreal (chr1, ndp + 20, nwds), d1, nwds)

write (6, '(a)') 'mpmdi(t1) ='
call mpmdi (t1, d4, i4)
write (6, '(1p,d25.15,i8/d25.15)') d4, i4, d4*2**i4
write (6, '(1p,2d25.15)') d1
call checkdp (t1, d4*2**i4, ndp)

write (6, '(a)') 'mpnrt(t1,i1) ='
call mpwrite (6, ndp + 20, ndp, mpnrt (t1,i1))
write (6, '(1p,2d25.15)') d1 ** (1.d0 / i1)
call checkdp (mpnrt (t1, i1), d1 ** (1.d0 / i1), ndp)

write (6, '(a)') 'mpprodd (t1,d2) ='
call mpwrite (6, ndp + 20, ndp, mpprodd (t1,d2))
write (6, '(1p,2d25.15)') d1*d2
call checkdp (mpprodd (t1,d2), d1 * d2, ndp)

write (6, '(a)') 'mpquotd (t1,d2) ='
call mpwrite (6, ndp + 20, ndp, mpquotd (t1,d2))
write (6, '(1p,2d25.15)') d1/d2
call checkdp (mpquotd (t1,d2), d1 / d2, ndp)

write (6, '(a)') 'mpreal (chr500, nwds) ='
call mpwrite (6, ndp + 20, ndp, mpreal (chr500, nwds))
read (chr500, *) d4
write (6, '(1p,2d25.15)') d4
call checkdp (mpreal (chr500, nwds), d4, ndp)

write (6, '(a)') 'mpreal (chrpi, ndp, nwds) ='
call mpwrite (6, ndp + 20, ndp, mpreal (chrpi, ndp, nwds))
read (chr500, *) d4
write (6, '(1p,2d25.15)') d4
call checkdp (mpreal (chrpi, ndp, nwds), d4, ndp)

write (6, '(a)') 'sign(t1,t2) ='
call mpwrite (6, ndp + 20, ndp, sign (t1,t2))
write (6, '(1p,2d25.15)') sign (d1, d2)
call checkdp (sign (t1, t2), sign (d1, d2), ndp)

write (6, '(a)') 'sin(t2) ='
call mpwrite (6, ndp + 20, ndp, sin (t2))
write (6, '(1p,2d25.15)') sin (d2)
call checkdp (sin (t2), sin (d2), ndp)

write (6, '(a)') 'sinh(t1) ='
call mpwrite (6, ndp + 20, ndp, sinh (t1))
write (6, '(1p,2d25.15)') sinh (d1)
call checkdp (sinh (t1), sinh (d1), ndp)

write (6, '(a)') 'sqrt(t1) ='
call mpwrite (6, ndp + 20, ndp, sqrt (t1))
write (6, '(1p,2d25.15)') sqrt (d1)
call checkdp (sqrt (t1), sqrt (d1), ndp)

write (6, '(a)') 'tan(t1) ='
call mpwrite (6, ndp + 20, ndp, tan (t2))
write (6, '(1p,2d25.15)') tan (d2)
call checkdp (tan (t2), tan (d2), ndp)

write (6, '(a)') 'tanh(t1) ='
call mpwrite (6, ndp + 20, ndp, tanh (t1))
write (6, '(1p,2d25.15)') tanh (d1)
call checkdp (tanh (t1), tanh (d1), ndp)

! write (6, '(a)') 'zeta(t1) ='
! call mpwrite (6, ndp + 20, ndp, zeta (t1))
! write (6, '(1p,2d25.15)') zeta (d1)
! call checkdp (zeta (t1), zeta (d1))

write (6, '(/a/)') 'Complex data operations:'

write (6, '(a)') 'addition: z1+z2 ='
call mpwrite (6, ndp + 20, ndp, z1 + z2)
write (6, 1) dc1 + dc2
call checkdc (z1 + z2, dc1 + dc2, ndp)

write (6, '(a)') 'addition: z1+e2 ='
call mpwrite (6, ndp + 20, ndp, z1 + e2)
write (6, 1) dc1 + e2
call checkdc (z1 + e2, dc1 + e2, ndp)

write (6, '(a)') 'addition: e1+z2 ='
call mpwrite (6, ndp + 20, ndp, e1 + z2)
write (6, 1) e1 + dc2
call checkdc (e1 + z2, e1 + dc2, ndp)

call mpwrite (6, ndp + 20, ndp, z1 + ec2)
write (6, 1) dc1 + ec2
call checkdc (z1 + ec2, dc1 + ec2, ndp)

write (6, '(a)') 'addition: ec1+z2 ='
call mpwrite (6, ndp + 20, ndp, ec1 + z2)
write (6, 1) ec1 + dc2
call checkdc (ec1 + z2, ec1 + dc2, ndp)

write (6, '(a)') 'addition: z1+t2 ='
call mpwrite (6, ndp + 20, ndp, z1 + t2)
write (6, 1) dc1 + d2
call checkdc (z1 + t2, dc1 + d2, ndp)

write (6, '(a)') 'addition: t1+z2 ='
call mpwrite (6, ndp + 20, ndp, t1 + z2)
write (6, 1) d1 + dc2
call checkdc (t1 + z2, d1 + dc2, ndp)

write (6, '(a)') 'subtraction: z1-z2 ='
call mpwrite (6, ndp + 20, ndp, z1 - z2)
write (6, 1) dc1 - dc2
call checkdc (z1 - z2, dc1 - dc2, ndp)

write (6, '(a)') 'subtraction: z1-e2 ='
call mpwrite (6, ndp + 20, ndp, z1 - e2)
write (6, 1) dc1 - e2
call checkdc (z1 - e2, dc1 - e2, ndp)

write (6, '(a)') 'subtraction: e1-z2 ='
call mpwrite (6, ndp + 20, ndp, e1 - z2)
write (6, 1) e1 - dc2
call checkdc (e1 - z2, e1 - dc2, ndp)

write (6, '(a)') 'subtraction: z1-ec2 ='
call mpwrite (6, ndp + 20, ndp, z1 - ec2)
write (6, 1) dc1 - ec2
call checkdc (z1 - ec2, dc1 - ec2, ndp)

write (6, '(a)') 'subtraction: ec1-z2 ='
call mpwrite (6, ndp + 20, ndp, ec1 - z2)
write (6, 1) ec1 - dc2
call checkdc (ec1 - z2, ec1 - dc2, ndp)

write (6, '(a)') 'subtraction: z1-t2 ='
call mpwrite (6, ndp + 20, ndp, z1 - t2)
write (6, 1) dc1 - d2
call checkdc (z1 - t2, dc1 - d2, ndp)

write (6, '(a)') 'subtraction: t1-z2 ='
call mpwrite (6, ndp + 20, ndp, t1 - z2)
write (6, 1) d1 - dc2
call checkdc (t1 - z2, d1 - dc2, ndp)

write (6, '(a)') 'multiplication: z1*z2 ='
call mpwrite (6, ndp + 20, ndp, z1 * z2)
write (6, 1) dc1 * dc2
call checkdc (z1 * z2, dc1 * dc2, ndp)

write (6, '(a)') 'multiplication: z1*e2 ='
call mpwrite (6, ndp + 20, ndp, z1 * e2)
write (6, 1) dc1 * e2
call checkdc (z1 * e2, dc1 * e2, ndp)

write (6, '(a)') 'multiplication: e1*z2 ='
call mpwrite (6, ndp + 20, ndp, e1 * z2)
write (6, 1) e1 * dc2
call checkdc (e1 * z2, e1 * dc2, ndp)

write (6, '(a)') 'multiplication: z1*ec2 ='
call mpwrite (6, ndp + 20, ndp, z1 * ec2)
write (6, 1) dc1 * ec2
call checkdc (z1 * ec2, dc1 * ec2, ndp)

write (6, '(a)') 'multiplication: ec1*z2 ='
call mpwrite (6, ndp + 20, ndp, ec1 * z2)
write (6, 1) ec1 * dc2
call checkdc (ec1 * z2, ec1 * dc2, ndp)

write (6, '(a)') 'multiplication: z1*t2 ='
call mpwrite (6, ndp + 20, ndp, z1 * t2)
write (6, 1) dc1 * d2
call checkdc (z1 * t2, dc1 * d2, ndp)

write (6, '(a)') 'multiplication: t1*z2 ='
call mpwrite (6, ndp + 20, ndp, t1 * z2)
write (6, 1) d1 * dc2
call checkdc (t1 * z2, d1 * dc2, ndp)

write (6, '(a)') 'multiplication: z1*e2 ='
call mpwrite (6, ndp + 20, ndp, z1 * e2)
write (6, 1) dc1 * e2
call checkdc (z1 * e2, dc1 * e2, ndp)

write (6, '(a)') 'multiplication: e1*z2 ='
call mpwrite (6, ndp + 20, ndp, e1 * z2)
write (6, 1) e1 * dc2
call checkdc (e1 * z2, e1 * dc2, ndp)

write (6, '(a)') 'division: z1/z2 ='
call mpwrite (6, ndp + 20, ndp, z1 / z2)
write (6, 1) dc1 / dc2
call checkdc (z1 / z2, dc1 / dc2, ndp)

write (6, '(a)') 'division: z1/e2 ='
call mpwrite (6, ndp + 20, ndp, z1 / e2)
write (6, 1) dc1 / e2
call checkdc (z1 / e2, dc1 / e2, ndp)

write (6, '(a)') 'division: e1/z2 ='
call mpwrite (6, ndp + 20, ndp, e1 / z2)
write (6, 1) e1 / dc2
call checkdc (e1 / z2, e1 / dc2, ndp)

write (6, '(a)') 'division: z1/ec2 ='
call mpwrite (6, ndp + 20, ndp, z1 / ec2)
write (6, 1) dc1 / ec2
call checkdc (z1 / ec2, dc1 / ec2, ndp)

write (6, '(a)') 'division: ec1/z2 ='
call mpwrite (6, ndp + 20, ndp, ec1 / z2)
write (6, 1) ec1 / dc2
call checkdc (ec1 / z2, ec1 / dc2, ndp)

write (6, '(a)') 'division: z1/t2 ='
call mpwrite (6, ndp + 20, ndp, z1 / t2)
write (6, 1) dc1 / d2
call checkdc (z1 / t2, dc1 / d2, ndp)

write (6, '(a)') 'division: t1/z2 ='
call mpwrite (6, ndp + 20, ndp, t1 / z2)
write (6, 1) d1 / dc2
call checkdc (t1 / z2, d1 / dc2, ndp)

write (6, '(a)') 'exponentiation: z1**i1 ='
call mpwrite (6, ndp + 20, ndp, z1 ** i1)
write (6, 1) dc1 ** i1
call checkdc (z1 ** i1, dc1 ** i1, ndp)

write (6, '(a)') 'exponentiation: z1**z2 ='
call mpwrite (6, ndp + 20, ndp, z1 ** z2)
write (6, 1) dc1 ** dc2
call checkdc (z1 ** z2, dc1 ** dc2, ndp)

write (6, '(a)') 'exponentiation: t1**z2 ='
call mpwrite (6, ndp + 20, ndp, t1 ** z2)
write (6, 1) d1 ** dc2
call checkdc (t1 ** z2, d1 ** dc2, ndp)

write (6, '(a)') 'exponentiation: z1**t2 ='
call mpwrite (6, ndp + 20, ndp, z1 ** t2)
write (6, 1) dc1 ** d2
call checkdc (z1 ** t2, dc1 ** d2, ndp)

write (6, '(a)') 'equal test: z1 == z2, e1 == z2, z1 == e2, ec1 == z2, z1 == ec2'
l1 = z1 == z2; l2 = dc1 == dc2; l3 = e1 == z2; l4 = e1 == dc2; l5 = z1 == e2
l6 = ec1 == e2; l7 = ec1 == z2; l8 = ec1 == dc2; l9 = z1 == ec2; l10 = dc1 == ec2
write (6, '(10l4)') l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)
call checkl (l7, l8); call checkl (l9, l10)

write (6, '(a)') 'not-equal test: z1 /= z2, e1 /= z2, z1 /= e2, ec1 /= z2, z1 /= ec2'
l1 = z1 /= z2; l2 = dc1 /= dc2; l3 = e1 /= z2; l4 = e1 /= dc2; l5 = z1 /= e2
l6 = ec1 /= e2; l7 = ec1 /= z2; l8 = ec1 /= dc2; l9 = z1 /= ec2; l10 = dc1 /= ec2
write (6, '(10l4)') l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
call checkl (l1, l2); call checkl (l3, l4); call checkl (l5, l6)
call checkl (l7, l8); call checkl (l9, l10)

write (6, '(a)') 'abs(z2) ='
call mpwrite (6, ndp + 20, ndp, abs (z2))
write (6, '(1p,2d25.15)') abs (dc2)
call checkdp (abs (z2), abs (dc2), ndp)

write (6, '(a)') 'aimag(z1) ='
call mpwrite (6, ndp + 20, ndp, aimag (z1))
write (6, '(1p,2d25.15)') aimag (dc1)
call checkdp (aimag (z1), aimag (dc1), ndp)

write (6, '(a)') 'conjg(z1) ='
call mpwrite (6, ndp + 20, ndp, conjg (z1))
write (6, '(1p,2d25.15)') conjg (dc1)
call checkdc (conjg (z1), conjg (dc1), ndp)

write (6, '(a)') 'cos(z2) ='
call mpwrite (6, ndp + 20, ndp, cos (z2))
write (6, '(1p,2d25.15)') cos (dc2)
call checkdc (cos (z2), cos (dc2), ndp)

write (6, '(a)') 'dcmplx(z1) ='
write (6, '(1p,2d25.15)') cmplx (dc1, kind=mprknd)
write (6, '(1p,2d25.15)') dc1
call checkdc (z1, dc1, ndp)

write (6, '(a)') 'exp(z1) ='
call mpwrite (6, ndp + 20, ndp, exp (z1))
write (6, '(1p,2d25.15)') exp (dc1)
call checkdc (exp (z1), exp (dc1), ndp)

write (6, '(a)') 'log(z1) ='
call mpwrite (6, ndp + 20, ndp, log (z1))
write (6, '(1p,2d25.15)') log (dc1)
call checkdc (log (z1), log (dc1), ndp)

write (6, '(a)') 'mpcmplx(t1,t2) ='
call mpwrite (6, ndp + 20, ndp, mpcmplx (t1, t2, nwds))
write (6, '(1p,2d25.15)') cmplx (d1, d2, mprknd)
call checkdc (mpcmplx (t1, t2, nwds), cmplx (d1, d2, mprknd), ndp)

write (6, '(a)') 'mpreal(z1) ='
call mpwrite (6, ndp + 20, ndp, mpreal (z1, nwds))
write (6, '(1p,2d25.15)') dble (dc1)
call checkdp (mpreal (z1, nwds), dble (dc1), ndp)

write (6, '(a)') 'sin(z2) ='
call mpwrite (6, ndp + 20, ndp, sin (z2))
write (6, '(1p,2d25.15)') sin (dc2)
call checkdc (sin (z2), sin (dc2), ndp)

write (6, '(a)') 'sqrt(z1) ='
call mpwrite (6, ndp + 20, ndp, sqrt (z1))
write (6, '(1p,2d25.15)') sqrt (dc1)
call checkdc (sqrt (z1), sqrt (dc1), ndp)

write (6, '(a)') 'sqrt(z2) ='
call mpwrite (6, ndp + 20, ndp, sqrt (z2))
write (6, '(1p,2d25.15)') sqrt (dc2)
call checkdc (sqrt (z2), sqrt (dc2), ndp)

stop
end

subroutine checkdp (t1, d1, ndp)
use mpmodule
implicit none
type (mp_real) t1
real (8) d1, dtol
parameter (dtol = 1d-14)
integer ndp

if (abs ((dble (t1) - d1) / d1) > dtol) then
  write (6, '(a)') 'error:'
  call mpwrite (6, ndp + 20, ndp, t1)
  write (6, '(1p,d30.16)') d1
  stop
endif

return
end

subroutine checkdc (z1, dc1, ndp)
use mpmodule
implicit none
type (mp_complex) z1
complex (8) dc1
real (8) dtol
parameter (dtol = 1d-14)
integer ndp

if (abs ((dcmplx (z1) - dc1) / dc1) > dtol) then
  write (6, '("error")')
  call mpwrite (6, ndp + 20, ndp, z1)
  write (6, '(1p,2d30.16)') dc1
  stop
endif

return
end

subroutine checkl (l1, l2)
logical l1, l2

if (l1 .neqv. l2) then
  write (6, '("error:",2l4)') l1, l2
  stop
endif 

return
end
