module chisqr

    use wavepack_kinds, only: rk
    implicit none

contains

!**************************************************************************
!      dcdflib: available at http://www.netlib.org/random/
!
!   ***note (c.t.) *** this is a subset of dcdflib, which only includes
!            the routines used by the chi square cumulative distribution.
!
! modified: november 1999 by arjan van dijk to include implicit none and
!           to convert all routines to double precision.
!**************************************************************************

!**************************************************************************
!                                           dcdflib
!                  library of fortran routines for cumulative distribution
!                       functions, inverses, and other parameters
!                                    (february, 1994)
!                            compiled and written by:
!                                     barry w. brown
!                                      james lovato
!                                      kathy russell
!
!
!                              legalities
!
! code that appeared  in an    acm  publication  is subject  to    their
! algorithms policy:
!
!        submittal of  an  algorithm  for publication  in   one of   the  acm
!        transactions implies that unrestricted use  of the algorithm within  a
!        computer is permissible.     general permission  to copy and  distribute
!        the algorithm without fee is granted provided that the copies        are not
!        made  or     distributed for  direct   commercial  advantage.    the acm
!        copyright notice and the title of the publication and its date appear,
!        and  notice is given that copying  is by permission of the association
!        for computing machinery.  to copy otherwise, or to republish, requires
!        a fee and/or specific permission.
!
!        krogh, f.  algorithms        policy.  acm  tran.   math.  softw.   13(1987),
!        183-186.
!
! we place the dcdflib code that we have written in the public domain.
!
!                                     no warranty
!
!        we provide absolutely        no warranty  of any  kind  either  expressed or
!        implied,  including but   not limited to,  the  implied  warranties of
!        merchantability and fitness for a particular purpose.        the entire risk
!        as to the quality and performance of the program is  with you.  should
!        this program prove  defective, you assume  the cost  of  all necessary
!        servicing, repair or correction.
!
!        in no        event  shall the university  of texas or  any  of its component
!        institutions including m. d. anderson hospital be liable  to you for
!        damages, including any  lost profits, lost monies,   or other special,
!        incidental   or  consequential damages   arising     out  of  the use or
!        inability to use (including but not limited to loss of data or data or
!        its analysis being  rendered inaccurate or  losses sustained  by third
!        parties) the program.
!
!        (above no warranty modified from the gnu no warranty statement.)
!

    subroutine cdfchi(which, p, q, x, df, status, bound)
!**********************************************************************
!
!      subroutine cdfchi( which, p, q, x, df, status, bound )
!             cumulative distribution function
!             chi-square distribution
!
!
!                         function
!
!
!     calculates any one parameter of the chi-square
!     distribution given values for the others.
!
!
!                         arguments
!
!
!     which --> integer indicating which of the next three argument
!             values is to be calculated from the others.
!             legal range: 1..3
!             iwhich = 1 : calculate p and q from x and df
!             iwhich = 2 : calculate x from p,q and df
!             iwhich = 3 : calculate df from p,q and x
!                  integer which
!
!     p <--> the integral from 0 to x of the chi-square
!            distribution.
!            input range: [0, 1].
!                  double precision p
!
!     q <--> 1-p.
!            input range: (0, 1].
!            p + q = 1.0.
!                  double precision q
!
!     x <--> upper limit of integration of the non-central
!            chi-square distribution.
!            input range: [0, +infinity).
!            search range: [0,1e300]
!                  double precision x
!
!     df <--> degrees of freedom of the
!             chi-square distribution.
!             input range: (0, +infinity).
!             search range: [ 1e-300, 1e300]
!                  double precision df
!
!     status <-- 0 if calculation completed correctly
!             -i if input parameter number i is out of range
!              1 if answer appears to be lower than lowest
!                search bound
!              2 if answer appears to be higher than greatest
!                search bound
!              3 if p + q .ne. 1
!             10 indicates error returned from cumgam.  see
!                references in cdfgam
!                  integer status
!
!     bound <-- undefined if status is 0
!
!             bound exceeded by parameter number i if status
!             is negative.
!
!             lower search bound if status is 1.
!
!             upper search bound if status is 2.
!
!
!                         method
!
!
!     formula    26.4.19   of abramowitz  and       stegun, handbook  of
!     mathematical functions   (1966) is used   to reduce the chisqure
!     distribution to the incomplete distribution.
!
!     computation of other parameters involve a seach for a value that
!     produces  the desired  value  of p.   the search relies  on  the
!     monotinicity of p with the other parameter.
!
!**********************************************************************
        implicit none
!     .. parameters ..
        real(rk) :: tol
        parameter(tol=1.0d-8)
        real(rk) :: atol
        parameter(atol=1.0d-50)
        real(rk) :: zero, inf
        parameter(zero=1.0d-300, inf=1.0d300)
!     ..
!     .. scalar arguments ..
        real(rk) :: bound, df, p, q, x
        integer status, which
!     ..
!     .. local scalars ..
        real(rk) :: fx, cum, ccum, pq, porq
        logical qhi, qleft, qporq
!     ..
!     .. external functions ..
        real(rk) :: spmpar
        external spmpar
!     ..
!     .. external subroutines ..
        external dinvr, dstinv, cumchi
!     ..
!     .. executable statements ..
!
!     check arguments
!
        if (.not. ((which < 1) .or. (which > 3))) goto 30
        if (.not. (which < 1)) goto 10
        bound = 1.0d0
        goto 20

10      bound = 3.0d0
20      status = -1
        return

30      if (which == 1) goto 70
!
!     p
!
        if (.not. ((p < 0.0d0) .or. (p > 1.0d0))) goto 60
        if (.not. (p < 0.0d0)) goto 40
        bound = 0.0d0
        goto 50

40      bound = 1.0d0
50      status = -2
        return

60      continue
70      if (which == 1) goto 110
!
!     q
!
        if (.not. ((q <= 0.0d0) .or. (q > 1.0d0))) goto 100
        if (.not. (q <= 0.0d0)) goto 80
        bound = 0.0d0
        goto 90

80      bound = 1.0d0
90      status = -3
        return

100     continue
110     if (which == 2) goto 130
!
!     x
!
        if (.not. (x < 0.0d0)) goto 120
        bound = 0.0d0
        status = -4
        return

120     continue
130     if (which == 3) goto 150
!
!     df
!
        if (.not. (df <= 0.0d0)) goto 140
        bound = 0.0d0
        status = -5
        return

140     continue
150     if (which == 1) goto 190
!
!     p + q
!
        pq = p + q
        if (.not. (abs(((pq) - 0.5d0) - 0.5d0) > (3.0d0*spmpar(1)))) goto 180
        if (.not. (pq < 0.0d0)) goto 160
        bound = 0.0d0
        goto 170

160     bound = 1.0d0
170     status = 3
        return

180     continue
190     if (which == 1) goto 220
!
!     select the minimum of p or q
!
        qporq = p <= q
        if (.not. (qporq)) goto 200
        porq = p
        goto 210

200     porq = q
210     continue
!
!     calculate answers
!
220     if ((1) == (which)) then
!
!     calculating p and q
!
            status = 0
            call cumchi(x, df, p, q)
            if (porq > 1.5d0) then
                status = 10
                return

            end if

        else if ((2) == (which)) then
!
!     calculating x
!
            x = 5.0d0
            call dstinv(0.0d0, inf, 0.5d0, 0.5d0, 5.0d0, atol, tol)
            status = 0
            call dinvr(status, x, fx, qleft, qhi)
230         if (.not. (status == 1)) goto 270
            call cumchi(x, df, cum, ccum)
            if (.not. (qporq)) goto 240
            fx = cum - p
            goto 250

240         fx = ccum - q
250         if (.not. ((fx + porq) > 1.5d0)) goto 260
            status = 10
            return

260         call dinvr(status, x, fx, qleft, qhi)
            goto 230

270         if (.not. (status == -1)) goto 300
            if (.not. (qleft)) goto 280
            status = 1
            bound = 0.0d0
            goto 290

280         status = 2
            bound = inf
290         continue
300         continue

        else if ((3) == (which)) then
!
!     calculating df
!
            df = 5.0d0
            call dstinv(zero, inf, 0.5d0, 0.5d0, 5.0d0, atol, tol)
            status = 0
            call dinvr(status, df, fx, qleft, qhi)
310         if (.not. (status == 1)) goto 350
            call cumchi(x, df, cum, ccum)
            if (.not. (qporq)) goto 320
            fx = cum - p
            goto 330

320         fx = ccum - q
330         if (.not. ((fx + porq) > 1.5d0)) goto 340
            status = 10
            return

340         call dinvr(status, df, fx, qleft, qhi)
            goto 310

350         if (.not. (status == -1)) goto 380
            if (.not. (qleft)) goto 360
            status = 1
            bound = zero
            goto 370

360         status = 2
            bound = inf
370         continue
380     end if

        return

    end subroutine cdfchi

    subroutine cumchi(x, df, cum, ccum)
!**********************************************************************
!
!     subroutine function cumchi(x,df,cum,ccum)
!             cumulative of the chi-square distribution
!
!
!                         function
!
!
!     calculates the cumulative chi-square distribution.
!
!
!                         arguments
!
!
!     x       --> upper limit of integration of the
!               chi-square distribution.
!                             x is double precision
!
!     df      --> degrees of freedom of the
!               chi-square distribution.
!                             df is double precision
!
!     cum <-- cumulative chi-square distribution.
!                             cum is double precisio
!
!     ccum <-- compliment of cumulative chi-square distribution.
!                             ccum is double precisi
!
!
!                         method
!
!
!     calls incomplete gamma function (cumgam)
!
!**********************************************************************
        implicit none
!     .. scalar arguments ..
        real(rk) :: df, x, cum, ccum
!     ..
!     .. local scalars ..
        real(rk) :: a, xx
!     ..
!     .. external subroutines ..
        external cumgam
!     ..
!     .. executable statements ..
        a = df*0.5d0
        xx = x*0.5d0
        call cumgam(xx, a, cum, ccum)
        return

    end subroutine cumchi

    subroutine cumgam(x, a, cum, ccum)
!**********************************************************************
!
!     subroutine cumgam(x,a,cum,ccum)
!           double precision cumulative incomplete gamma distribution
!
!
!                         function
!
!
!     computes   the  cumulative of    the     incomplete   gamma
!     distribution, i.e., the integral from 0 to x of
!          (1/gam(a))*exp(-t)*t**(a-1) dt
!     where gam(a) is the complete gamma function of a, i.e.,
!          gam(a) = integral from 0 to infinity of
!                  exp(-t)*t**(a-1) dt
!
!
!                         arguments
!
!
!     x --> the upper limit of integration of the incomplete gamma.
!                            x is double precision
!
!     a --> the shape parameter of the incomplete gamma.
!                            a is double precision
!
!     cum <-- cumulative incomplete gamma distribution.
!                          cum is double precision
!
!     ccum <-- compliment of cumulative incomplete gamma distribution.
!                            ccum is double precisio
!
!
!                         method
!
!
!     calls the routine gratio.
!
!**********************************************************************
        implicit none
!
!     ..
!     .. scalar arguments ..
        real(rk) :: a, x, cum, ccum
!     ..
!     .. external routines ..
        external gratio
!     ..
!     .. executable statements ..
        if (.not. (x <= 0.0d0)) goto 10
        cum = 0.0d0
        ccum = 1.0d0
        return

10      call gratio(a, x, cum, ccum, 0)

!     call gratio routine

        return

    end subroutine cumgam
    
    subroutine dinvr(status, x, fx, qleft, qhi)
!**********************************************************************
!
!     subroutine dinvr(status, x, fx, qleft, qhi)
!          double precision
!          bounds the zero of the function and invokes zror
!                  reverse communication
!
!
!                         function
!
!
!     bounds the    function  and  invokes  zror   to perform the   zero
!     finding.  stinvr   must  have   been  called  before this     routine
!     in order to set its parameters.
!
!
!                         arguments
!
!
!     status <--> at the beginning of a zero finding problem, status
!               should be set to 0 and invr invoked.       (the value
!               of parameters other than x will be ignored on this cal
!
!               when invr needs the function evaluated, it will set
!               status to 1 and return.  the value of the function
!               should be set in fx and invr again called without
!               changing any of its other parameters.
!
!               when invr has finished without error, it will return
!               with status 0.  in that case x is approximately a root
!               of f(x).
!
!               if invr cannot bound the function, it returns status
!               -1 and sets qleft and qhi.
!                    integer status
!
!     x <-- the value of x at which f(x) is to be evaluated.
!                    double precision x
!
!     fx --> the value of f(x) calculated when invr returns with
!            status = 1.
!                    double precision fx
!
!     qleft <-- defined only if qmfinv returns .false.  in that
!          case it is .true. if the stepping search terminated
!          unsucessfully at small.  if it is .false. the search
!          terminated unsucessfully at big.
!                  qleft is logical
!
!     qhi <-- defined only if qmfinv returns .false.  in that
!          case it is .true. if f(x) .gt. y at the termination
!          of the search and .false. if f(x) .lt. y at the
!          termination of the search.
!                  qhi is logical

!
!**********************************************************************
        implicit none
!     .. scalar arguments ..
        real(rk) :: fx, x, zabsst, zabsto, zbig, zrelst, zrelto, zsmall, zstpmu
        integer status
        logical qhi, qleft
!     ..
!     .. local scalars ..
        real(rk) :: absstp, abstol, big, fbig, fsmall, relstp, reltol, small, step, &
            stpmul, xhi, xlb, xlo, xsave, xub, yy, zx, zy, zz
        integer i99999
        logical qbdd, qcond, qdum1, qdum2, qincr, qlim, qup
!     ..
!     .. external subroutines ..
        external dstzr, dzror
!     ..
!     .. intrinsic functions ..
        intrinsic abs, max, min
!     ..
!     .. statement functions ..
        logical qxmon
!     ..
!     .. save statement ..
        save
!     ..
!     .. statement function definitions ..
        qxmon(zx, zy, zz) = zx <= zy .and. zy <= zz
!     ..
!     .. executable statements ..

        if (status > 0) goto 310

        qcond = .not. qxmon(small, x, big)
        if (qcond) stop ' small, x, big not monotone in invr'
        xsave = x
!
!     see that small and big bound the zero and set qincr
!
        x = small
!     get-function-value
        assign 10 to i99999
        goto 300

10      fsmall = fx
        x = big
!     get-function-value
        assign 20 to i99999
        goto 300

20      fbig = fx
        qincr = fbig > fsmall
        if (.not. (qincr)) goto 50
        if (fsmall <= 0.0d0) goto 30
        status = -1
        qleft = .true.
        qhi = .true.
        return

30      if (fbig >= 0.0d0) goto 40
        status = -1
        qleft = .false.
        qhi = .false.
        return

40      goto 80

50      if (fsmall >= 0.0d0) goto 60
        status = -1
        qleft = .true.
        qhi = .false.
        return

60      if (fbig <= 0.0d0) goto 70
        status = -1
        qleft = .false.
        qhi = .true.
        return

70      continue
80      x = xsave
        step = max(absstp, relstp*abs(x))
!      yy = f(x) - y
!     get-function-value
        assign 90 to i99999
        goto 300

90      yy = fx
        if (.not. (yy == 0.0d0)) goto 100
        status = 0
        return

100     qup = (qincr .and. (yy < 0.0d0)) .or. (.not. qincr .and. (yy > 0.0d0))
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     handle case in which we must step higher
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (.not. (qup)) goto 170
        xlb = xsave
        xub = min(xlb + step, big)
        goto 120

110     if (qcond) goto 150
!      yy = f(xub) - y
120     x = xub
!     get-function-value
        assign 130 to i99999
        goto 300

130     yy = fx
        qbdd = (qincr .and. (yy >= 0.0d0)) .or. (.not. qincr .and. (yy <= 0.0d0))
        qlim = xub >= big
        qcond = qbdd .or. qlim
        if (qcond) goto 140
        step = stpmul*step
        xlb = xub
        xub = min(xlb + step, big)
140     goto 110

150     if (.not. (qlim .and. .not. qbdd)) goto 160
        status = -1
        qleft = .false.
        qhi = .not. qincr
        x = big
        return

160     goto 240
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     handle case in which we must step lower
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
170     xub = xsave
        xlb = max(xub - step, small)
        goto 190

180     if (qcond) goto 220
!      yy = f(xlb) - y
190     x = xlb
!     get-function-value
        assign 200 to i99999
        goto 300

200     yy = fx
        qbdd = (qincr .and. (yy <= 0.0d0)) .or. (.not. qincr .and. (yy >= 0.0d0))
        qlim = xlb <= small
        qcond = qbdd .or. qlim
        if (qcond) goto 210
        step = stpmul*step
        xub = xlb
        xlb = max(xub - step, small)
210     goto 180

220     if (.not. (qlim .and. .not. qbdd)) goto 230
        status = -1
        qleft = .true.
        qhi = qincr
        x = small
        return

230     continue
240     call dstzr(xlb, xub, abstol, reltol)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     if we reach here, xlb and xub bound the zero of f.
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        status = 0
        goto 260

250     if (.not. (status == 1)) goto 290
260     call dzror(status, x, fx, xlo, xhi, qdum1, qdum2)
        if (.not. (status == 1)) goto 280
!     get-function-value
        assign 270 to i99999
        goto 300

270     continue
280     goto 250

290     x = xlo
        status = 0
        return

        entry dstinv(zsmall, zbig, zabsst, zrelst, zstpmu, zabsto, zrelto)
!**********************************************************************
!
!      subroutine dstinv( small, big, absstp, relstp, stpmul,
!     +                  abstol, reltol )
!      double precision - set inverse finder - reverse communication
!
!
!                         function
!
!
!     concise description - given a monotone function f finds x
!     such that f(x) = y.  uses reverse communication -- see invr.
!     this routine sets quantities needed by invr.
!
!          more precise description of invr -
!
!     f must be a monotone function, the results of qmfinv are
!     otherwise undefined.  qincr must be .true. if f is non-
!     decreasing and .false. if f is non-increasing.
!
!     qmfinv will return .true. if and only if f(small) and
!     f(big) bracket y, i. e.,
!          qincr is .true. and f(small).le.y.le.f(big) or
!          qincr is .false. and f(big).le.y.le.f(small)
!
!     if qmfinv returns .true., then the x returned satisfies
!     the following condition.  let
!             tol(x) = max(abstol,reltol*abs(x))
!     then if qincr is .true.,
!          f(x-tol(x)) .le. y .le. f(x+tol(x))
!     and if qincr is .false.
!          f(x-tol(x)) .ge. y .ge. f(x+tol(x))
!
!
!                         arguments
!
!
!     small --> the left endpoint of the interval to be
!          searched for a solution.
!                  small is double precision
!
!     big --> the right endpoint of the interval to be
!          searched for a solution.
!                  big is double precision
!
!     absstp, relstp --> the initial step size in the search
!          is max(absstp,relstp*abs(x)). see algorithm.
!                  absstp is double precision
!                  relstp is double precision
!
!     stpmul --> when a step doesn't bound the zero, the step
!              size is multiplied by stpmul and another step
!              taken.  a popular value is 2.0
!                  double precision stpmul
!
!     abstol, reltol --> two numbers that determine the accuracy
!          of the solution.  see function for a precise definition.
!                  abstol is double precision
!                  reltol is double precision
!
!
!                         method
!
!
!     compares f(x) with y for the input value of x then uses qincr
!     to determine whether to step left or right to bound the
!     desired x.  the initial step size is
!          max(absstp,relstp*abs(s)) for the input value of x.
!     iteratively steps right or left until it bounds x.
!     at each step which doesn't bound x, the step size is doubled.
!     the routine is careful never to step beyond small or big.  if
!     it hasn't bounded x at small or big, qmfinv returns .false.
!     after setting qleft and qhi.
!
!     if x is successfully bounded then algorithm r of the paper
!     'two efficient algorithms with guaranteed convergence for
!     finding a zero of a function' by j. c. p. bus and
!     t. j. dekker in acm transactions on mathematical
!     software, volume 1, no. 4 page 330 (dec. '75) is employed
!     to find the zero of the function f(x)-y. this is routine
!     qrzero.
!
!**********************************************************************
        small = zsmall
        big = zbig
        absstp = zabsst
        relstp = zrelst
        stpmul = zstpmu
        abstol = zabsto
        reltol = zrelto
        return

!      stop '*** execution flowing into flecs procedures ***'
!     to get-function-value
300     status = 1
        return

310     continue
        goto i99999

    end subroutine dinvr

    subroutine dzror(status, x, fx, xlo, xhi, qleft, qhi)
!**********************************************************************
!
!     subroutine dzror(status, x, fx, xlo, xhi, qleft, qhi)
!     double precision zero of a function -- reverse communication
!
!
!                         function
!
!
!     performs the zero finding.  stzror must have been called before
!     this routine in order to set its parameters.
!
!
!                         arguments
!
!
!     status <--> at the beginning of a zero finding problem, status
!               should be set to 0 and zror invoked.       (the value
!               of other parameters will be ignored on this call.)
!
!               when zror needs the function evaluated, it will set
!               status to 1 and return.  the value of the function
!               should be set in fx and zror again called without
!               changing any of its other parameters.
!
!               when zror has finished without error, it will return
!               with status 0.  in that case (xlo,xhi) bound the answe
!
!               if zror finds an error (which implies that f(xlo)-y an
!               f(xhi)-y have the same sign, it returns status -1.  in
!               this case, xlo and xhi are undefined.
!                    integer status
!
!     x <-- the value of x at which f(x) is to be evaluated.
!                    double precision x
!
!     fx --> the value of f(x) calculated when zror returns with
!            status = 1.
!                    double precision fx
!
!     xlo <-- when zror returns with status = 0, xlo bounds the
!             inverval in x containing the solution below.
!                    double precision xlo
!
!     xhi <-- when zror returns with status = 0, xhi bounds the
!             inverval in x containing the solution above.
!                    double precision xhi
!
!     qleft <-- .true. if the stepping search terminated unsucessfully
!              at xlo.  if it is .false. the search terminated
!              unsucessfully at xhi.
!                  qleft is logical
!
!     qhi <-- .true. if f(x) .gt. y at the termination of the
!              search and .false. if f(x) .lt. y at the
!              termination of the search.
!                  qhi is logical
!
!**********************************************************************
        implicit none
!     .. scalar arguments ..
        real(rk) :: fx, x, xhi, xlo, zabstl, zreltl, zxhi, zxlo
        integer status
        logical qhi, qleft
!     ..
!     .. save statement ..
        save
!     ..
!     .. local scalars ..
        real(rk) :: a, abstol, b, c, d, fa, fb, fc, fd, fda, fdb, m, mb, p, q, &
            reltol, tol, w, xxhi, xxlo, zx
        integer ext, i99999
        logical first, qrzero
!     ..
!     .. intrinsic functions ..
        intrinsic abs, max, sign
!     ..
!     .. statement functions ..
        real(rk) :: ftol
!     ..
!     .. statement function definitions ..
        ftol(zx) = 0.5d0*max(abstol, reltol*abs(zx))
!     ..
!     .. executable statements ..

        if (status > 0) goto 280
        xlo = xxlo
        xhi = xxhi
        b = xlo
        x = xlo
!     get-function-value
        assign 10 to i99999
        goto 270

10      fb = fx
        xlo = xhi
        a = xlo
        x = xlo
!     get-function-value
        assign 20 to i99999
        goto 270
!
!     check that f(zxlo) < 0 < f(zxhi)  or
!              f(zxlo) > 0 > f(zxhi)
!
20      if (.not. (fb < 0.0d0)) goto 40
        if (.not. (fx < 0.0d0)) goto 30
        status = -1
        qleft = fx < fb
        qhi = .false.
        return

30      continue
40      if (.not. (fb > 0.0d0)) goto 60
        if (.not. (fx > 0.0d0)) goto 50
        status = -1
        qleft = fx > fb
        qhi = .true.
        return

50      continue
60      fa = fx
!
        first = .true.
70      c = a
        fc = fa
        ext = 0
80      if (.not. (abs(fc) < abs(fb))) goto 100
        if (.not. (c /= a)) goto 90
        d = a
        fd = fa
90      a = b
        fa = fb
        xlo = c
        b = xlo
        fb = fc
        c = a
        fc = fa
100     tol = ftol(xlo)
        m = (c + b)*.5d0
        mb = m - b
        if (.not. (abs(mb) > tol)) goto 240
        if (.not. (ext > 3)) goto 110
        w = mb
        goto 190

110     tol = sign(tol, mb)
        p = (b - a)*fb
        if (.not. (first)) goto 120
        q = fa - fb
        first = .false.
        goto 130

120     fdb = (fd - fb)/(d - b)
        fda = (fd - fa)/(d - a)
        p = fda*p
        q = fdb*fa - fda*fb
130     if (.not. (p < 0.0d0)) goto 140
        p = -p
        q = -q
140     if (ext == 3) p = p*2.0d0
        if (.not. ((p*1.0d0) == 0.0d0 .or. p <= (q*tol))) goto 150
        w = tol
        goto 180

150     if (.not. (p < (mb*q))) goto 160
        w = p/q
        goto 170

160     w = mb
170     continue
180     continue
190     d = a
        fd = fa
        a = b
        fa = fb
        b = b + w
        xlo = b
        x = xlo
!     get-function-value
        assign 200 to i99999
        goto 270

200     fb = fx
        if (.not. ((fc*fb) >= 0.0d0)) goto 210
        goto 70

210     if (.not. (w == mb)) goto 220
        ext = 0
        goto 230

220     ext = ext + 1
230     goto 80

240     xhi = c
        qrzero = (fc >= 0.0d0 .and. fb <= 0.0d0) .or. (fc < 0.0d0 .and. fb >= 0.0d0)
        if (.not. (qrzero)) goto 250
        status = 0
        goto 260

250     status = -1
260     return

        entry dstzr(zxlo, zxhi, zabstl, zreltl)
!**********************************************************************
!
!     subroutine dstzr( xlo, xhi, abstol, reltol )
!     double precision set zero finder - reverse communication version
!
!
!                         function
!
!
!
!     sets quantities needed by zror.  the function of zror
!     and the quantities set is given here.
!
!     concise description - given a function f
!     find xlo such that f(xlo) = 0.
!
!          more precise description -
!
!     input condition. f is a double precision function of a single
!     double precision argument and xlo and xhi are such that
!          f(xlo)*f(xhi)  .le.        0.0
!
!     if the input condition is met, qrzero returns .true.
!     and output values of xlo and xhi satisfy the following
!          f(xlo)*f(xhi)  .le. 0.
!          abs(f(xlo)  .le. abs(f(xhi)
!          abs(xlo-xhi)  .le. tol(x)
!     where
!          tol(x) = max(abstol,reltol*abs(x))
!
!     if this algorithm does not find xlo and xhi satisfying
!     these conditions then qrzero returns .false.  this
!     implies that the input condition was not met.
!
!
!                         arguments
!
!
!     xlo --> the left endpoint of the interval to be
!           searched for a solution.
!                  xlo is double precision
!
!     xhi --> the right endpoint of the interval to be
!           for a solution.
!                  xhi is double precision
!
!     abstol, reltol --> two numbers that determine the accuracy
!                    of the solution.  see function for a
!                    precise definition.
!                  abstol is double precision
!                  reltol is double precision
!
!
!                         method
!
!
!     algorithm r of the paper 'two efficient algorithms with
!     guaranteed convergence for finding a zero of a function'
!     by j. c. p. bus and t. j. dekker in acm transactions on
!     mathematical software, volume 1, no. 4 page 330
!     (dec. '75) is employed to find the zero of f(x)-y.
!
!**********************************************************************
        xxlo = zxlo
        xxhi = zxhi
        abstol = zabstl
        reltol = zreltl
        return

!      stop '*** execution flowing into flecs procedures ***'
!     to get-function-value
270     status = 1
        return

280     continue
        goto i99999

    end subroutine dzror
    real(rk) function erf(x)
!-----------------------------------------------------------------------
!             evaluation of the real error function
!-----------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: x
!     ..
!     .. local scalars ..
        real(rk) :: ax, bot, c, t, top, x2
!     ..
!     .. local arrays ..
        real(rk) :: a(5), b(3), p(8), q(8), r(5), s(4)
!     ..
!     .. intrinsic functions ..
        intrinsic abs, exp, sign
!     ..
!     .. data statements ..
!-------------------------
!-------------------------
!-------------------------
!-------------------------
        data c/.564189583547756d0/
        data a(1)/.771058495001320d-04/, a(2)/-.133733772997339d-02/, &
            a(3)/.323076579225834d-01/, a(4)/.479137145607681d-01/, a(5)/.128379167095513d+00/
        data b(1)/.301048631703895d-02/, b(2)/.538971687740286d-01/, &
            b(3)/.375795757275549d+00/
        data p(1)/-1.36864857382717d-07/, p(2)/5.64195517478974d-01/, &
            p(3)/7.21175825088309d+00/, p(4)/4.31622272220567d+01/, &
            p(5)/1.52989285046940d+02/, p(6)/3.39320816734344d+02/, &
            p(7)/4.51918953711873d+02/, p(8)/3.00459261020162d+02/
        data q(1)/1.00000000000000d+00/, q(2)/1.27827273196294d+01/, &
            q(3)/7.70001529352295d+01/, q(4)/2.77585444743988d+02/, &
            q(5)/6.38980264465631d+02/, q(6)/9.31354094850610d+02/, &
            q(7)/7.90950925327898d+02/, q(8)/3.00459260956983d+02/
        data r(1)/2.10144126479064d+00/, r(2)/2.62370141675169d+01/, &
            r(3)/2.13688200555087d+01/, r(4)/4.65807828718470d+00/, r(5)/2.82094791773523d-01/
        data s(1)/9.41537750555460d+01/, s(2)/1.87114811799590d+02/, &
            s(3)/9.90191814623914d+01/, s(4)/1.80124575948747d+01/
!     ..
!     .. executable statements ..
!-------------------------
        ax = abs(x)
        if (ax > 0.5d0) goto 10
        t = x*x
        top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0d0
        bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0d0
        erf = x*(top/bot)
        return
!
10      if (ax > 4.0d0) goto 20
        top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax + p(7))*ax + p(8)
        bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax + q(7))*ax + q(8)
        erf = 0.5d0 + (0.5d0 - exp(-x*x)*top/bot)
        if (x < 0.0d0) erf = -erf
        return
!
20      if (ax >= 5.8d0) goto 30
        x2 = x*x
        t = 1.0d0/x2
        top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
        bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0d0
        erf = (c - top/(x2*bot))/ax
        erf = 0.5d0 + (0.5d0 - exp(-x2)*erf)
        if (x < 0.0d0) erf = -erf
        return
!
30      erf = sign(1.0d0, x)
        return

    end function erf
    real(rk) function erfc1(ind, x)
!-----------------------------------------------------------------------
!         evaluation of the complementary error function
!
!          erfc1(ind,x) = erfc(x)          if ind = 0
!          erfc1(ind,x) = exp(x*x)*erfc(x)   otherwise
!-----------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: x
        integer ind
!     ..
!     .. local scalars ..
        real(rk) :: ax, bot, c, e, t, top, w
!     ..
!     .. local arrays ..
        real(rk) :: a(5), b(3), p(8), q(8), r(5), s(4)
!     ..
!     .. external functions ..
        real(rk) :: exparg
        external exparg
!     ..
!     .. intrinsic functions ..
        intrinsic abs, dble, exp
!     ..
!     .. data statements ..
!-------------------------
!-------------------------
!-------------------------
!-------------------------
        data c/.564189583547756d0/
        data a(1)/.771058495001320d-04/, a(2)/-.133733772997339d-02/, &
            a(3)/.323076579225834d-01/, a(4)/.479137145607681d-01/, a(5)/.128379167095513d+00/
        data b(1)/.301048631703895d-02/, b(2)/.538971687740286d-01/, &
            b(3)/.375795757275549d+00/
        data p(1)/-1.36864857382717d-07/, p(2)/5.64195517478974d-01/, &
            p(3)/7.21175825088309d+00/, p(4)/4.31622272220567d+01/, &
            p(5)/1.52989285046940d+02/, p(6)/3.39320816734344d+02/, &
            p(7)/4.51918953711873d+02/, p(8)/3.00459261020162d+02/
        data q(1)/1.00000000000000d+00/, q(2)/1.27827273196294d+01/, &
            q(3)/7.70001529352295d+01/, q(4)/2.77585444743988d+02/, &
            q(5)/6.38980264465631d+02/, q(6)/9.31354094850610d+02/, &
            q(7)/7.90950925327898d+02/, q(8)/3.00459260956983d+02/
        data r(1)/2.10144126479064d+00/, r(2)/2.62370141675169d+01/, &
            r(3)/2.13688200555087d+01/, r(4)/4.65807828718470d+00/, r(5)/2.82094791773523d-01/
        data s(1)/9.41537750555460d+01/, s(2)/1.87114811799590d+02/, &
            s(3)/9.90191814623914d+01/, s(4)/1.80124575948747d+01/
!     ..
!     .. executable statements ..
!-------------------------
!
!                   abs(x) .le. 0.5
!
        ax = abs(x)
        if (ax > 0.5d0) goto 10
        t = x*x
        top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0d0
        bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0d0
        erfc1 = 0.5d0 + (0.5d0 - x*(top/bot))
        if (ind /= 0) erfc1 = exp(t)*erfc1
        return
!
!                0.5 .lt. abs(x) .le. 4
!
10      if (ax > 4.0d0) goto 20
        top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax + p(7))*ax + p(8)
        bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax + q(7))*ax + q(8)
        erfc1 = top/bot
        goto 40
!
!                    abs(x) .gt. 4
!
20      if (x <= -5.6d0) goto 60
        if (ind /= 0) goto 30
        if (x > 100.0d0) goto 70
        if (x*x > -exparg(1)) goto 70
!
30      t = (1.0d0/x)**2
        top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
        bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0d0
        erfc1 = (c - t*top/bot)/ax
!
!                    final assembly
!
40      if (ind == 0) goto 50
        if (x < 0.0d0) erfc1 = 2.0d0*exp(x*x) - erfc1
        return

50      w = dble(x)*dble(x)
        t = w
        e = w - dble(t)
        erfc1 = ((0.5d0 + (0.5d0 - e))*exp(-t))*erfc1
        if (x < 0.0d0) erfc1 = 2.0d0 - erfc1
        return
!
!             limit value for large negative x
!
60      erfc1 = 2.0d0
        if (ind /= 0) erfc1 = 2.0d0*exp(x*x)
        return
!
!             limit value for large positive x
!                  when ind = 0
!
70      erfc1 = 0.0d0
        return

    end function erfc1
    real(rk) function exparg(l)
!--------------------------------------------------------------------
!     if l = 0 then  exparg(l) = the largest positive w for which
!     exp(w) can be computed.
!
!     if l is nonzero then  exparg(l) = the largest negative w for
!     which the computed value of exp(w) is nonzero.
!
!     note... only an approximate value for exparg(l) is needed.
!--------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        integer l
!     ..
!     .. local scalars ..
        real(rk) :: lnb
        integer b, m
!     ..
!     .. external functions ..
        integer ipmpar
        external ipmpar
!     ..
!     .. intrinsic functions ..
        intrinsic dble, log
!     ..
!     .. executable statements ..
!
        b = ipmpar(4)
        if (b /= 2) goto 10
        lnb = .69314718055995d0
        goto 40

10      if (b /= 8) goto 20
        lnb = 2.0794415416798d0
        goto 40

20      if (b /= 16) goto 30
        lnb = 2.7725887222398d0
        goto 40

30      lnb = log(dble(b))
!
40      if (l == 0) goto 50
        m = ipmpar(9) - 1
        exparg = 0.99999d0*(m*lnb)
        return

50      m = ipmpar(10)
        exparg = 0.99999d0*(m*lnb)
        return

    end function exparg

    real(rk) function gam1(a)
!     ------------------------------------------------------------------
!     computation of 1/gamma(a+1) - 1  for -0.5 .le. a .le. 1.5
!     ------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: a
!     ..
!     .. local scalars ..
        real(rk) :: bot, d, s1, s2, t, top, w
!     ..
!     .. local arrays ..
        real(rk) :: p(7), q(5), r(9)
!     ..
!     .. data statements ..
!     -------------------
!     -------------------
!     -------------------
!     -------------------
        data p(1)/.577215664901533d+00/, p(2)/-.409078193005776d+00/, &
            p(3)/-.230975380857675d+00/, p(4)/.597275330452234d-01/, &
            p(5)/.766968181649490d-02/, p(6)/-.514889771323592d-02/, p(7)/.589597428611429d-03/
        data q(1)/.100000000000000d+01/, q(2)/.427569613095214d+00/, &
            q(3)/.158451672430138d+00/, q(4)/.261132021441447d-01/, q(5)/.423244297896961d-02/
        data r(1)/-.422784335098468d+00/, r(2)/-.771330383816272d+00/, &
            r(3)/-.244757765222226d+00/, r(4)/.118378989872749d+00/, &
            r(5)/.930357293360349d-03/, r(6)/-.118290993445146d-01/, &
            r(7)/.223047661158249d-02/, r(8)/.266505979058923d-03/, r(9)/-.132674909766242d-03/
        data s1/.273076135303957d+00/, s2/.559398236957378d-01/
!     ..
!     .. executable statements ..
!     -------------------
        t = a
        d = a - 0.5d0
        if (d > 0.0d0) t = d - 0.5d0
        if (t) 40, 10, 20
!
10      gam1 = 0.0d0
        return
!
20      top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t + p(2))*t + p(1)
        bot = (((q(5)*t + q(4))*t + q(3))*t + q(2))*t + 1.0d0
        w = top/bot
        if (d > 0.0d0) goto 30
        gam1 = a*w
        return

30      gam1 = (t/a)*((w - 0.5d0) - 0.5d0)
        return
!
40      top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t + r(4))*t + r(3))*t + r(2))*t + r(1)
        bot = (s2*t + s1)*t + 1.0d0
        w = top/bot
        if (d > 0.0d0) goto 50
        gam1 = a*((w + 0.5d0) + 0.5d0)
        return

50      gam1 = t*w/a
        return

    end function gam1
    real(rk) function gamma(a)
!-----------------------------------------------------------------------
!
!         evaluation of the gamma function for real arguments
!
!                      -----------
!
!     gamma(a) is assigned the value 0 when the gamma function cannot
!     be computed.
!
!-----------------------------------------------------------------------
!     written by alfred h. morris, jr.
!          naval surface weapons center
!          dahlgren, virginia
!-----------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: a
!     ..
!     .. local scalars ..
        real(rk) :: bot, d, g, lnx, pi, r1, r2, r3, r4, r5, s, t, top, w, x, z
        integer i, j, m, n
!     ..
!     .. local arrays ..
        real(rk) :: p(7), q(7)
!     ..
!     .. external functions ..
        real(rk) :: exparg, spmpar
        external exparg, spmpar
!     ..
!     .. intrinsic functions ..
        intrinsic abs, dble, log, exp, int, mod, sin
!     ..
!     .. data statements ..
!--------------------------
!     d = 0.5*(ln(2*pi) - 1)
!--------------------------
!--------------------------
!--------------------------
        data pi/3.1415926535898d0/
        data d/.41893853320467274178d0/
        data p(1)/.539637273585445d-03/, p(2)/.261939260042690d-02/, &
            p(3)/.204493667594920d-01/, p(4)/.730981088720487d-01/, &
            p(5)/.279648642639792d+00/, p(6)/.553413866010467d+00/, p(7)/1.0d0/
        data q(1)/-.832979206704073d-03/, q(2)/.470059485860584d-02/, &
            q(3)/.225211131035340d-01/, q(4)/-.170458969313360d+00/, &
            q(5)/-.567902761974940d-01/, q(6)/.113062953091122d+01/, q(7)/1.0d0/
        data r1/.820756370353826d-03/, r2/-.595156336428591d-03/, &
            r3/.793650663183693d-03/, r4/-.277777777770481d-02/, r5/.833333333333333d-01/
!     ..
!     .. executable statements ..
!--------------------------
        gamma = 0.0d0
        x = a
        if (abs(a) >= 15.0d0) goto 110
!-----------------------------------------------------------------------
!            evaluation of gamma(a) for abs(a) .lt. 15
!-----------------------------------------------------------------------
        t = 1.0d0
        m = int(a) - 1
!
!     let t be the product of a-j when a .ge. 2
!
        if (m) 40, 30, 10
10      do j = 1, m
            x = x - 1.0d0
            t = x*t
        end do
30      x = x - 1.0d0
        goto 80
!
!     let t be the product of a+j when a .lt. 1
!
40      t = a
        if (a > 0.0d0) goto 70
        m = -m - 1
        if (m == 0) goto 60
        do j = 1, m
            x = x + 1.0d0
            t = x*t
        end do
60      x = (x + 0.5d0) + 0.5d0
        t = x*t
        if (t == 0.0d0) return
!
70      continue
!
!     the following code checks if 1/t can overflow. this
!     code may be omitted if desired.
!
        if (abs(t) >= 1.d-30) goto 80
        if (abs(t)*spmpar(3) <= 1.0001d0) return
        gamma = 1.0d0/t
        return
!
!     compute gamma(1 + x) for  0 .le. x .lt. 1
!
80      top = p(1)
        bot = q(1)
        do i = 2, 7
            top = p(i) + x*top
            bot = q(i) + x*bot
        end do
        gamma = top/bot
!
!     termination
!
        if (a < 1.0d0) goto 100
        gamma = gamma*t
        return

100     gamma = gamma/t
        return
!-----------------------------------------------------------------------
!            evaluation of gamma(a) for abs(a) .ge. 15
!-----------------------------------------------------------------------
110     if (abs(a) >= 1.d3) return
        if (a > 0.0d0) goto 120
        x = -a
        n = x
        t = x - n
        if (t > 0.9d0) t = 1.0d0 - t
        s = sin(pi*t)/pi
        if (mod(n, 2) == 0) s = -s
        if (s == 0.0d0) return
!
!     compute the modified asymptotic sum
!
120     t = 1.0d0/(x*x)
        g = ((((r1*t + r2)*t + r3)*t + r4)*t + r5)/x
!
!     one may replace the next statement with  lnx = alog(x)
!     but less accuracy will normally be obtained.
!
        lnx = log(x)
!
!     final assembly
!
        z = x
        g = (d + g) + (z - 0.5d0)*(lnx - 1.d0)
        w = g
        t = g - dble(w)
        if (w > 0.99999d0*exparg(0)) return
        gamma = exp(w)*(1.0d0 + t)
        if (a < 0.0d0) gamma = (1.0d0/(gamma*s))/x
        return

    end function gamma

    subroutine gratio(a, x, ans, qans, ind)
! ----------------------------------------------------------------------
!        evaluation of the incomplete gamma ratio functions
!                    p(a,x) and q(a,x)
!
!                   ----------
!
!     it is assumed that a and x are nonnegative, where a and x
!     are not both 0.
!
!     ans and qans are variables. gratio assigns ans the value
!     p(a,x) and qans the value q(a,x). ind may be any integer.
!     if ind = 0 then the user is requesting as much accuracy as
!     possible (up to 14 significant digits). otherwise, if
!     ind = 1 then accuracy is requested to within 1 unit of the
!     6-th significant digit, and if ind .ne. 0,1 then accuracy
!     is requested to within 1 unit of the 3rd significant digit.
!
!     error return ...
!        ans is assigned the value 2 when a or x is negative,
!     when a*x = 0, or when p(a,x) and q(a,x) are indeterminant.
!     p(a,x) and q(a,x) are computationally indeterminant when
!     x is exceedingly close to a and a is extremely large.
! ----------------------------------------------------------------------
!     written by alfred h. morris, jr.
!        naval surface weapons center
!        dahlgren, virginia
!     --------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: a, ans, qans, x
        integer ind
!     ..
!     .. local scalars ..
        real(rk) :: a2n, a2nm1, acc, alog10, am0, amn, an, an0, apn, b2n, b2nm1, &
            c, c0, c1, c2, c3, c4, c5, c6, cma, d10, d20, d30, d40, d50, d60, d70, e, e0, &
            g, h, j, l, r, rt2pin, rta, rtpi, rtx, s, sum, t, t1, third, tol, twoa, u, w, x0, y, z
        integer i, iop, m, max, n
!     ..
!     .. local arrays ..
        real(rk) :: acc0(3), big(3), d0(13), d1(12), d2(10), d3(8), d4(6), d5(4), d6(2), e00(3), wk(20), x00(3)
!     ..
!     .. external functions ..
        real(rk) :: erf, erfc1, gam1, gamma, rexp, rlog, spmpar
        external erf, erfc1, gam1, gamma, rexp, rlog, spmpar

!     ..
!     .. data statements ..
!     --------------------
!     --------------------
!     alog10 = ln(10)
!     rt2pin = 1/sqrt(2*pi)
!     rtpi   = sqrt(pi)
!     --------------------
!     --------------------
!     --------------------
!     --------------------
!     --------------------
!     --------------------
!     --------------------
!     --------------------
!     --------------------
        data acc0(1)/5.d-15/, acc0(2)/5.d-7/, acc0(3)/5.d-4/
        data big(1)/20.0d0/, big(2)/14.0d0/, big(3)/10.0d0/
        data e00(1)/.25d-3/, e00(2)/.25d-1/, e00(3)/.14d0/
        data x00(1)/31.0d0/, x00(2)/17.0d0/, x00(3)/9.7d0/
        data alog10/2.30258509299405d0/
        data rt2pin/.398942280401433d0/
        data rtpi/1.77245385090552d0/
        data third/.333333333333333d0/
        data d0(1)/.833333333333333d-01/, d0(2)/-.148148148148148d-01/, &
            d0(3)/.115740740740741d-02/, d0(4)/.352733686067019d-03/, &
            d0(5)/-.178755144032922d-03/, d0(6)/.391926317852244d-04/, &
            d0(7)/-.218544851067999d-05/, d0(8)/-.185406221071516d-05/, &
            d0(9)/.829671134095309d-06/, d0(10)/-.176659527368261d-06/, &
            d0(11)/.670785354340150d-08/, d0(12)/.102618097842403d-07/, d0(13)/-.438203601845335d-08/
        data d10/-.185185185185185d-02/, d1(1)/-.347222222222222d-02/, &
            d1(2)/.264550264550265d-02/, d1(3)/-.990226337448560d-03/, &
            d1(4)/.205761316872428d-03/, d1(5)/-.401877572016461d-06/, &
            d1(6)/-.180985503344900d-04/, d1(7)/.764916091608111d-05/, &
            d1(8)/-.161209008945634d-05/, d1(9)/.464712780280743d-08/, &
            d1(10)/.137863344691572d-06/, d1(11)/-.575254560351770d-07/, d1(12)/.119516285997781d-07/
        data d20/.413359788359788d-02/, d2(1)/-.268132716049383d-02/, &
            d2(2)/.771604938271605d-03/, d2(3)/.200938786008230d-05/, &
            d2(4)/-.107366532263652d-03/, d2(5)/.529234488291201d-04/, &
            d2(6)/-.127606351886187d-04/, d2(7)/.342357873409614d-07/, &
            d2(8)/.137219573090629d-05/, d2(9)/-.629899213838006d-06/, &
            d2(10)/.142806142060642d-06/
        data d30/.649434156378601d-03/, d3(1)/.229472093621399d-03/, &
            d3(2)/-.469189494395256d-03/, d3(3)/.267720632062839d-03/, &
            d3(4)/-.756180167188398d-04/, d3(5)/-.239650511386730d-06/, &
            d3(6)/.110826541153473d-04/, d3(7)/-.567495282699160d-05/, d3(8)/.142309007324359d-05/
        data d40/-.861888290916712d-03/, d4(1)/.784039221720067d-03/, &
            d4(2)/-.299072480303190d-03/, d4(3)/-.146384525788434d-05/, &
            d4(4)/.664149821546512d-04/, d4(5)/-.396836504717943d-04/, d4(6)/.113757269706784d-04/
        data d50/-.336798553366358d-03/, d5(1)/-.697281375836586d-04/, &
            d5(2)/.277275324495939d-03/, d5(3)/-.199325705161888d-03/, d5(4)/.679778047793721d-04/
        data d60/.531307936463992d-03/, d6(1)/-.592166437353694d-03/, d6(2)/.270878209671804d-03/
        data d70/.344367606892378d-03/
!     ..
!     .. executable statements ..
!     --------------------
!     ****** e is a machine dependent constant. e is the smallest
!            floating point number for which 1.0 + e .gt. 1.0 .
!
        e = spmpar(1)
!
!     --------------------
        if (a < 0.0d0 .or. x < 0.0d0) goto 430
        if (a == 0.0d0 .and. x == 0.0d0) goto 430
        if (a*x == 0.0d0) goto 420
!
        iop = ind + 1
        if (iop /= 1 .and. iop /= 2) iop = 3
        acc = dmax1(acc0(iop), e)
        e0 = e00(iop)
        x0 = x00(iop)
!
!            select the appropriate algorithm
!
        if (a >= 1.0d0) goto 10
        if (a == 0.5d0) goto 390
        if (x < 1.1d0) goto 160
        t1 = a*log(x) - x
        u = a*exp(t1)
        if (u == 0.0d0) goto 380
        r = u*(1.0d0 + gam1(a))
        goto 250
!
10      if (a >= big(iop)) goto 30
        if (a > x .or. x >= x0) goto 20
        twoa = a + a
        m = int(twoa)
        if (twoa /= dble(m)) goto 20
        i = m/2
        if (a == dble(i)) goto 210
        goto 220

20      t1 = a*log(x) - x
        r = exp(t1)/gamma(a)
        goto 40
!
30      l = x/a
        if (l == 0.0d0) goto 370
        s = 0.5d0 + (0.5d0 - l)
        z = rlog(l)
        if (z >= 700.0d0/a) goto 410
        y = a*z
        rta = sqrt(a)
        if (abs(s) <= e0/rta) goto 330
        if (abs(s) <= 0.4d0) goto 270
!
        t = (1.0d0/a)**2
        t1 = (((0.75d0*t - 1.0d0)*t + 3.5d0)*t - 105.0d0)/(a*1260.0d0)
        t1 = t1 - y
        r = rt2pin*rta*exp(t1)
!
40      if (r == 0.0d0) goto 420
        if (x <= dmax1(a, alog10)) goto 50
        if (x < x0) goto 250
        goto 100
!
!               taylor series for p/r
!
50      apn = a + 1.0d0
        t = x/apn
        wk(1) = t
        do n = 2, 20
            apn = apn + 1.0d0
            t = t*(x/apn)
            if (t <= 1.d-3) goto 70
            wk(n) = t
        end do
        n = 20
!
70      sum = t
        tol = 0.5d0*acc
80      apn = apn + 1.0d0
        t = t*(x/apn)
        sum = sum + t
        if (t > tol) goto 80
!
        max = n - 1
        do m = 1, max
            n = n - 1
            sum = sum + wk(n)
        end do
        ans = (r/a)*(1.0d0 + sum)
        qans = 0.5d0 + (0.5d0 - ans)
        return
!
!               asymptotic expansion
!
100     amn = a - 1.0d0
        t = amn/x
        wk(1) = t
        do n = 2, 20
            amn = amn - 1.0d0
            t = t*(amn/x)
            if (abs(t) <= 1.d-3) goto 120
            wk(n) = t
        end do
        n = 20
!
120     sum = t
130     if (abs(t) <= acc) goto 140
        amn = amn - 1.0d0
        t = t*(amn/x)
        sum = sum + t
        goto 130
!
140     max = n - 1
        do m = 1, max
            n = n - 1
            sum = sum + wk(n)
        end do
        qans = (r/x)*(1.0d0 + sum)
        ans = 0.5d0 + (0.5d0 - qans)
        return
!
!             taylor series for p(a,x)/x**a
!
160     an = 3.0d0
        c = x
        sum = x/(a + 3.0d0)
        tol = 3.0d0*acc/(a + 1.0d0)
170     an = an + 1.0d0
        c = -c*(x/an)
        t = c/(a + an)
        sum = sum + t
        if (abs(t) > tol) goto 170
        j = a*x*((sum/6.0d0 - 0.5d0/(a + 2.0d0))*x + 1.0d0/(a + 1.0d0))
!
        z = a*log(x)
        h = gam1(a)
        g = 1.0d0 + h
        if (x < 0.25d0) goto 180
        if (a < x/2.59d0) goto 200
        goto 190

180     if (z > -.13394d0) goto 200
!
190     w = exp(z)
        ans = w*g*(0.5d0 + (0.5d0 - j))
        qans = 0.5d0 + (0.5d0 - ans)
        return
!
200     l = rexp(z)
        w = 0.5d0 + (0.5d0 + l)
        qans = (w*j - l)*g - h
        if (qans < 0.0d0) goto 380
        ans = 0.5d0 + (0.5d0 - qans)
        return
!
!             finite sums for q when a .ge. 1
!               and 2*a is an integer
!
210     sum = exp(-x)
        t = sum
        n = 1
        c = 0.0d0
        goto 230
!
220     rtx = sqrt(x)
        sum = erfc1(0, rtx)
        t = exp(-x)/(rtpi*rtx)
        n = 0
        c = -0.5d0
!
230     if (n == i) goto 240
        n = n + 1
        c = c + 1.0d0
        t = (x*t)/c
        sum = sum + t
        goto 230

240     qans = sum
        ans = 0.5d0 + (0.5d0 - qans)
        return
!
!              continued fraction expansion
!
250     tol = dmax1(5.0d0*e, acc)
        a2nm1 = 1.0d0
        a2n = 1.0d0
        b2nm1 = x
        b2n = x + (1.0d0 - a)
        c = 1.0d0
260     a2nm1 = x*a2n + c*a2nm1
        b2nm1 = x*b2n + c*b2nm1
        am0 = a2nm1/b2nm1
        c = c + 1.0d0
        cma = c - a
        a2n = a2nm1 + cma*a2n
        b2n = b2nm1 + cma*b2n
        an0 = a2n/b2n
        if (abs(an0 - am0) >= tol*an0) goto 260
!
        qans = r*an0
        ans = 0.5d0 + (0.5d0 - qans)
        return
!
!              general temme expansion
!
270     if (abs(s) <= 2.0d0*e .and. a*e*e > 3.28d-3) goto 430
        c = exp(-y)
        w = 0.5d0*erfc1(1, sqrt(y))
        u = 1.0d0/a
        z = sqrt(z + z)
        if (l < 1.0d0) z = -z
        if (iop - 2) 280, 290, 300
!
280     if (abs(s) <= 1.d-3) goto 340
        c0 = ((((((((((((d0(13)*z + d0(12))*z + d0(11))*z + d0(10))*z + d0(9))*z + d0(8))*z + d0(7))*z + d0(6))*z &
                  + d0(5))*z + d0(4))*z + d0(3))*z + d0(2))*z + d0(1))*z - third
        c1 = (((((((((((d1(12)*z + d1(11))*z + d1(10))*z + d1(9))*z + d1(8))*z + d1(7))*z + d1(6))*z + d1(5))*z &
                 + d1(4))*z + d1(3))*z + d1(2))*z + d1(1))*z + d10
        c2 = (((((((((d2(10)*z + d2(9))*z + d2(8))*z + d2(7))*z + d2(6))*z + d2(5))*z + d2(4))*z + &
                d2(3))*z + d2(2))*z + d2(1))*z + d20
        c3 = (((((((d3(8)*z + d3(7))*z + d3(6))*z + d3(5))*z + d3(4))*z + d3(3))*z + d3(2))*z + d3(1))*z + d30
        c4 = (((((d4(6)*z + d4(5))*z + d4(4))*z + d4(3))*z + d4(2))*z + d4(1))*z + d40
        c5 = (((d5(4)*z + d5(3))*z + d5(2))*z + d5(1))*z + d50
        c6 = (d6(2)*z + d6(1))*z + d60
        t = ((((((d70*u + c6)*u + c5)*u + c4)*u + c3)*u + c2)*u + c1)*u + c0
        goto 310
!
290     c0 = (((((d0(6)*z + d0(5))*z + d0(4))*z + d0(3))*z + d0(2))*z + d0(1))*z - third
        c1 = (((d1(4)*z + d1(3))*z + d1(2))*z + d1(1))*z + d10
        c2 = d2(1)*z + d20
        t = (c2*u + c1)*u + c0
        goto 310
!
300     t = ((d0(3)*z + d0(2))*z + d0(1))*z - third
!
310     if (l < 1.0d0) goto 320
        qans = c*(w + rt2pin*t/rta)
        ans = 0.5d0 + (0.5d0 - qans)
        return

320     ans = c*(w - rt2pin*t/rta)
        qans = 0.5d0 + (0.5d0 - ans)
        return
!
!             temme expansion for l = 1
!
330     if (a*e*e > 3.28d-3) goto 430
        c = 0.5d0 + (0.5d0 - y)
        w = (0.5d0 - sqrt(y)*(0.5d0 + (0.5d0 - y/3.0d0))/rtpi)/c
        u = 1.0d0/a
        z = sqrt(z + z)
        if (l < 1.0d0) z = -z
        if (iop - 2) 340, 350, 360
!
340     c0 = ((((((d0(7)*z + d0(6))*z + d0(5))*z + d0(4))*z + d0(3))*z + d0(2))*z + d0(1))*z - third
        c1 = (((((d1(6)*z + d1(5))*z + d1(4))*z + d1(3))*z + d1(2))*z + d1(1))*z + d10
        c2 = ((((d2(5)*z + d2(4))*z + d2(3))*z + d2(2))*z + d2(1))*z + d20
        c3 = (((d3(4)*z + d3(3))*z + d3(2))*z + d3(1))*z + d30
        c4 = (d4(2)*z + d4(1))*z + d40
        c5 = (d5(2)*z + d5(1))*z + d50
        c6 = d6(1)*z + d60
        t = ((((((d70*u + c6)*u + c5)*u + c4)*u + c3)*u + c2)*u + c1)*u + c0
        goto 310
!
350     c0 = (d0(2)*z + d0(1))*z - third
        c1 = d1(1)*z + d10
        t = (d20*u + c1)*u + c0
        goto 310
!
360     t = d0(1)*z - third
        goto 310
!
!                   special cases
!
370     ans = 0.0d0
        qans = 1.0d0
        return
!
380     ans = 1.0d0
        qans = 0.0d0
        return
!
390     if (x >= 0.25d0) goto 400
        ans = erf(sqrt(x))
        qans = 0.5d0 + (0.5d0 - ans)
        return

400     qans = erfc1(0, sqrt(x))
        ans = 0.5d0 + (0.5d0 - qans)
        return
!
410     if (abs(s) <= 2.0d0*e) goto 430
420     if (x <= a) goto 370
        goto 380
!
!                   error return
!
430     ans = 2.0d0
        return

    end subroutine gratio

    integer function ipmpar(i)
!-----------------------------------------------------------------------
!
!     ipmpar provides the integer machine constants for the computer
!     that is used. it is assumed that the argument i is an integer
!     having one of the values 1-10. ipmpar(i) has the value ...
!
!  integers.
!
!     assume integers are represented in the n-digit, base-a form
!
!             sign ( x(n-1)*a**(n-1) + ... + x(1)*a + x(0) )
!
!             where 0 .le. x(i) .lt. a for i=0,...,n-1.
!
!     ipmpar(1) = a, the base.
!
!     ipmpar(2) = n, the number of base-a digits.
!
!     ipmpar(3) = a**n - 1, the largest magnitude.
!
!  floating-point numbers.
!
!     it is assumed that the single and double precision floating
!     point arithmetics have the same base, say b, and that the
!     nonzero numbers are represented in the form
!
!             sign (b**e) * (x(1)/b + ... + x(m)/b**m)
!
!             where x(i) = 0,1,...,b-1 for i=1,...,m,
!             x(1) .ge. 1, and emin .le. e .le. emax.
!
!     ipmpar(4) = b, the base.
!
!  single-precision
!
!     ipmpar(5) = m, the number of base-b digits.
!
!     ipmpar(6) = emin, the smallest exponent e.
!
!     ipmpar(7) = emax, the largest exponent e.
!
!  double-precision
!
!     ipmpar(8) = m, the number of base-b digits.
!
!     ipmpar(9) = emin, the smallest exponent e.
!
!     ipmpar(10) = emax, the largest exponent e.
!
!-----------------------------------------------------------------------
!
!     to define this function for the computer being used, activate
!     the data statments for the computer by removing the c from
!     column 1. (all the other data statements should have c in
!     column 1.)
!
!-----------------------------------------------------------------------
!
!     ipmpar is an adaptation of the function i1mach, written by
!     p.a. fox, a.d. hall, and n.l. schryer (bell laboratories).
!     ipmpar was formed by a.h. morris (nswc). the constants are
!     from bell laboratories, nswc, and other sources.
!
!-----------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        integer i
!     ..
!     .. local arrays ..
        integer imach(10)
!     ..
!     .. data statements ..
!
!     machine constants for amdahl machines.
!
!     data imach( 1) /   2 /
!     data imach( 2) /  31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /  16 /
!     data imach( 5) /   6 /
!     data imach( 6) / -64 /
!     data imach( 7) /  63 /
!     data imach( 8) /  14 /
!     data imach( 9) / -64 /
!     data imach(10) /  63 /
!
!     machine constants for the at&t 3b series, at&t
!     pc 7300, and at&t 6300.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the burroughs 1700 system.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   33 /
!     data imach( 3) / 8589934591 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -256 /
!     data imach( 7) /  255 /
!     data imach( 8) /   60 /
!     data imach( 9) / -256 /
!     data imach(10) /  255 /
!
!     machine constants for the burroughs 5700 system.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /    8 /
!     data imach( 5) /   13 /
!     data imach( 6) /  -50 /
!     data imach( 7) /   76 /
!     data imach( 8) /   26 /
!     data imach( 9) /  -50 /
!     data imach(10) /   76 /
!
!     machine constants for the burroughs 6700/7700 systems.
!
!     data imach( 1) /      2 /
!     data imach( 2) /     39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /      8 /
!     data imach( 5) /     13 /
!     data imach( 6) /    -50 /
!     data imach( 7) /     76 /
!     data imach( 8) /     26 /
!     data imach( 9) / -32754 /
!     data imach(10) /  32780 /
!
!     machine constants for the cdc 6000/7000 series
!     60 bit arithmetic, and the cdc cyber 995 64 bit
!     arithmetic (nos operating system).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   48 /
!     data imach( 3) / 281474976710655 /
!     data imach( 4) /    2 /
!     data imach( 5) /   48 /
!     data imach( 6) / -974 /
!     data imach( 7) / 1070 /
!     data imach( 8) /   95 /
!     data imach( 9) / -926 /
!     data imach(10) / 1070 /
!
!     machine constants for the cdc cyber 995 64 bit
!     arithmetic (nos/ve operating system).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    48 /
!     data imach( 6) / -4096 /
!     data imach( 7) /  4095 /
!     data imach( 8) /    96 /
!     data imach( 9) / -4096 /
!     data imach(10) /  4095 /
!
!     machine constants for the cray 1, xmp, 2, and 3.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    47 /
!     data imach( 6) / -8189 /
!     data imach( 7) /  8190 /
!     data imach( 8) /    94 /
!     data imach( 9) / -8099 /
!     data imach(10) /  8190 /
!
!     machine constants for the data general eclipse s/200.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     machine constants for the harris 220.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   23 /
!     data imach( 3) / 8388607 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   38 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the honeywell 600/6000
!     and dps 8/70 series.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   63 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the hp 2100
!     3 word double precision option with ftn4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   39 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     machine constants for the hp 2100
!     4 word double precision option with ftn4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   55 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     machine constants for the hp 9000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -126 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the ibm 360/370 series,
!     the icl 2900, the itel as/6, the xerox sigma
!     5/7/9 and the sel systems 85/86.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     machine constants for the ibm pc.
!
!      data imach(1)/2/
!      data imach(2)/31/
!      data imach(3)/2147483647/
!      data imach(4)/2/
!      data imach(5)/24/
!      data imach(6)/-125/
!      data imach(7)/128/
!      data imach(8)/53/
!      data imach(9)/-1021/
!      data imach(10)/1024/
!
!     machine constants for the macintosh ii - absoft
!     macfortran ii.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the microvax - vms fortran.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the pdp-10 (ka processor).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   54 /
!     data imach( 9) / -101 /
!     data imach(10) /  127 /
!
!     machine constants for the pdp-10 (ki processor).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   62 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     machine constants for the pdp-11 fortran supporting
!     32-bit integer arithmetic.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     machine constants for the sequent balance 8000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for the silicon graphics iris-4d
!     series (mips r3000 processor).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     machine constants for ieee arithmetic machines, such as the at&t
!     3b series, motorola 68000 based machines (e.g. sun 3 and at&t
!     pc 7300), and 8087 based micros (e.g. ibm pc and at&t 6300).
!
        data imach(1)/2/
        data imach(2)/31/
        data imach(3)/2147483647/
        data imach(4)/2/
        data imach(5)/24/
        data imach(6)/-125/
        data imach(7)/128/
        data imach(8)/53/
        data imach(9)/-1021/
        data imach(10)/1024/
!
!     machine constants for the univac 1100 series.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   60 /
!     data imach( 9) /-1024 /
!     data imach(10) / 1023 /
!
!     machine constants for the vax 11/780.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
        ipmpar = imach(i)
        return

    end function ipmpar
    
    real(rk) function rexp(x)
!-----------------------------------------------------------------------
!            evaluation of the function exp(x) - 1
!-----------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: x
!     ..
!     .. local scalars ..
        real(rk) :: p1, p2, q1, q2, q3, q4, w
!     ..
!     .. intrinsic functions ..
        intrinsic abs, exp
!     ..
!     .. data statements ..
        data p1/.914041914819518d-09/, p2/.238082361044469d-01/, &
            q1/-.499999999085958d+00/, q2/.107141568980644d+00/, &
            q3/-.119041179760821d-01/, q4/.595130811860248d-03/
!     ..
!     .. executable statements ..
!-----------------------
        if (abs(x) > 0.15d0) goto 10
        rexp = x*(((p2*x + p1)*x + 1.0d0)/((((q4*x + q3)*x + q2)*x + q1)*x + 1.0d0))
        return
!
10      w = exp(x)
        if (x > 0.0d0) goto 20
        rexp = (w - 0.5d0) - 0.5d0
        return

20      rexp = w*(0.5d0 + (0.5d0 - 1.0d0/w))
        return

    end function rexp
    
    real(rk) function rlog(x)
!     -------------------
!     computation of  x - 1 - ln(x)
!     -------------------
        implicit none
!     .. scalar arguments ..
        real(rk) :: x
!     ..
!     .. local scalars ..
        real(rk) :: a, b, p0, p1, p2, q1, q2, r, t, u, w, w1
!     ..
!     .. intrinsic functions ..
        intrinsic dble, log
!     ..
!     .. data statements ..
!     -------------------
        data a/.566749439387324d-01/
        data b/.456512608815524d-01/
        data p0/.333333333333333d+00/, p1/-.224696413112536d+00/, p2/.620886815375787d-02/
        data q1/-.127408923933623d+01/, q2/.354508718369557d+00/
!     ..
!     .. executable statements ..
!     -------------------
        if (x < 0.61d0 .or. x > 1.57d0) goto 40
        if (x < 0.82d0) goto 10
        if (x > 1.18d0) goto 20
!
!              argument reduction
!
        u = (x - 0.5d0) - 0.5d0
        w1 = 0.0d0
        goto 30
!
10      u = dble(x) - 0.7d0
        u = u/0.7d0
        w1 = a - u*0.3d0
        goto 30
!
20      u = 0.75d0*dble(x) - 1.d0
        w1 = b + u/3.0d0
!
!             series expansion
!
30      r = u/(u + 2.0d0)
        t = r*r
        w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0d0)
        rlog = 2.0d0*t*(1.0d0/(1.0d0 - r) - r*w) + w1
        return
!
!
40      r = (x - 0.5d0) - 0.5d0
        rlog = r - log(x)
        return

    end function rlog
    
    real(rk) function spmpar(i)
!-----------------------------------------------------------------------
!
!     spmpar provides the single precision machine constants for
!     the computer being used. it is assumed that the argument
!     i is an integer having one of the values 1, 2, or 3. if the
!     single precision arithmetic being used has m base b digits and
!     its smallest and largest exponents are emin and emax, then
!
!        spmpar(1) = b**(1 - m), the machine precision,
!
!        spmpar(2) = b**(emin - 1), the smallest magnitude,
!
!        spmpar(3) = b**emax*(1 - b**(-m)), the largest magnitude.
!
!-----------------------------------------------------------------------
!     written by
!        alfred h. morris, jr.
!        naval surface warfare center
!        dahlgren virginia
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     modified by barry w. brown to return double precision machine
!     constants for the computer being used.  this modification was
!     made as part of converting bratio to double precision
!-----------------------------------------------------------------------
        implicit none
!     .. scalar arguments ..
        integer i
!     ..
!     .. local scalars ..
        real(rk) :: b, binv, bm1, one, w, z
        integer emax, emin, ibeta, m
!     ..
!     .. external functions ..
        integer ipmpar
        external ipmpar
!     ..
!     .. intrinsic functions ..
        intrinsic dble
!     ..
!     .. executable statements ..
!
        if (i > 1) goto 10
        b = ipmpar(4)
        m = ipmpar(8)
        spmpar = b**(1 - m)
        return
!
10      if (i > 2) goto 20
        b = ipmpar(4)
        emin = ipmpar(9)
        one = dble(1)
        binv = one/b
        w = b**(emin + 2)
        spmpar = ((w*binv)*binv)*binv
        return
!
20      ibeta = ipmpar(4)
        m = ipmpar(8)
        emax = ipmpar(10)
!
        b = ibeta
        bm1 = ibeta - 1
        one = dble(1)
        z = b**(m - 1)
        w = ((z - one)*b + bm1)/(b*z)
!
        z = b**(emax - 2)
        spmpar = ((w*z)*b)*b
        return

    end function spmpar

end module chisqr
