module wavepack

    use wavepack_kinds, only: rk
    use chisqr, only: cdfchi
    implicit none
    private

    public :: wavelet, wave_signif, rk

contains

    !****************************************************************************
    !
    ! wavepack:  routines to compute the wavelet transform of a time series,
    !            and significance levels.
    !
    ! written by:   christopher torrence and gilbert p. compo
    !
    ! available from:  http://paos.colorado.edu/research/wavelets/
    !
    ! requires the following packages:   cfftpack, chisqr
    !
    ! notes:
    !
    !  (1) all routines are written in single precision (double precision),
    !      except for chisqr, which requires double precision input.
    !      single precision should be sufficient for most applications.
    !
    !  (2) the cfftpack and chisqr routines were not written by us,
    !      and no guarentees are made as to their reliability or efficiency.
    !
    !  (3) no provision is made for output to files or to graphics packages.
    !      the user is expected to call these routines from within
    !      their own programs. see sample program "wavetest.f".
    !
    !  (4) little error checking is done. check your input carefully.
    !      the programs are not completely ansi compatible, so use caution.
    !
    !  (5) time series are currently limited to 65535 points.
    !
    !
    ! reference: torrence, c. and g. p. compo, 1998: a practical guide to
    !            wavelet analysis. <i>bull. amer. meteor. soc.</i>, 79, 61-78.
    !
    ! notice: please acknowledge the use of this software in any publications:
    !    ``wavelet software was provided by c. torrence and g. compo,
    !      and is available at url: http://paos.colorado.edu/research/wavelets/''.
    !
    !  copyright (c) 1998, christopher torrence
    ! this software may be used, copied, or redistributed as long as it is not
    ! sold and this copyright notice is reproduced on each copy made.  this
    ! routine is provided as is without any express or implied warranties
    ! whatsoever.
    !
    !
    ! modified: november 1999 by arjan van dijk to include implicit none and
    !           to convert all routines to double precision.
    !****************************************************************************

    !****************************************************************************
    ! wavelet: computes the wavelet transform of a time series,
    !          with appropriate parameters.
    !
    !
    ! inputs:
    !
    !  n [int] = the number of points in "y".
    !
    !  y [double precision] = the time series of length "n".
    !
    !  dt [double precision] = amount of time between each y value, i.e. the sampling time.
    !
    !  mother [int] = an integer giving the mother wavelet to use.
    !                   0='morlet'
    !                   1='paul'
    !                   2='dog' (derivative of gaussian)
    !               if (mother<0 or >2) then default is 'morlet'.
    !
    !  param [double precision] = mother wavelet parameter. if <0 then default is used.
    !            for 'morlet' this is k0 (wavenumber), default is 6.
    !            for 'paul' this is m (order), default is 4.
    !            for 'dog' this is m (m-th derivative), default is 2.
    !
    !
    !  s0 [double precision] = the smallest scale of the wavelet.  typically = 2*dt.
    !         note: for accurate reconstruction and variance computation
    !             set s0=dt for morlet; s0=dt/4 for paul
    !
    !  dj [double precision] = the spacing between discrete scales. typically = 0.25.
    !         a smaller # will give better scale resolution, but be slower.
    !
    !  jtot [int] = the # of scales.
    !             scales range from s0 up to s0*2^[(jtot-1)*dj],
    !             typically jtot=1+(log2(n dt/s0))/dj
    !
    !  npad [int] = the total number of points (including padding) to
    !             use for the wavelet transform. typically this is some
    !             power of 2. it must be greater or equal to "n".
    !             if npad>n, then zeroes are padded onto the end
    !             of the time series.
    !
    !
    ! outputs:
    !
    !  wave [dcmplx(n,jtot)] = 2d array of the real & imaginary parts
    !                 of the wavelet transform, versus time & scale.
    !                 cabs(wave) gives the wavelet amplitude,
    !                 atan2(aimag(wave),dble(wave)) gives wavelet phase.
    !                 the wavelet power spectrum is cabs(wave)**2.
    !
    !  scale [double precision(jtot)] = the wavelet scales that were used.
    !
    !  period [double precision(jtot)] = the "fourier" periods (in time units) corresponding
    !            to "scale".
    !
    !  coi [double precision(n)] = the e-folding factor used for the cone of influence.
    !
    !
    ! requires:   wave_function, cfftpack
    !
    !  copyright (c) 1998, christopher torrence and gilbert p. compo
    ! this software may be used, copied, or redistributed as long as it is not
    ! sold and this copyright notice is reproduced on each copy made.  this
    ! routine is provided as is without any express or implied warranties
    ! whatsoever.

    subroutine wavelet(n, y, dt, mother, param, s0, dj, jtot, npad, wave, scale, period, coi)
        implicit none

        integer n, mother, jtot, npad
        real(rk) :: y(n), dt, param, s0, dj, scale(jtot), period(jtot), coi(n)
        double complex wave(n, jtot)

        integer i, j, k, nk
        real(rk) :: ymean, freq1, pi, period1, coi1

!** initialize work arrays
        parameter(nk=65535)
        real(rk) :: wsave(4*nk + 15), kwave(nk)
        double complex yfft(nk), daughter(nk)

        pi = 4._rk*atan(1._rk)

        if (npad < n) then
            print *, '**wavelet: "npad" must be greater than or equal to "n"'
            return
        end if

        if ((mother < 0) .or. (mother > 2)) mother = 0

!** find the time-series mean & remove it
        ymean = 0._rk
        do i = 1, n
            ymean = ymean + y(i)
        end do
        ymean = ymean/n
        do i = 1, n
            yfft(i) = y(i) - ymean
        end do

!** if desired, pad with extra zeroes
        do i = n + 1, npad
            yfft(i) = 0._rk
        end do

!** find the fft of the time series [eqn(3)]
        call cffti(npad, wsave)
        call cfftf(npad, yfft, wsave)
        do k = 1, npad
            yfft(k) = yfft(k)/npad
        end do

!** construct the wavenumber array [eqn(5)]
        freq1 = 2._rk*pi/(dble(npad)*dt)
        kwave(1) = 0._rk
        do i = 2, npad/2 + 1
            kwave(i) = (dble(i) - 1._rk)*freq1
        end do
        do i = npad/2 + 2, npad
            kwave(i) = -kwave(npad - i + 2)
        end do

!**----- main wavelet loop
!      print '(a8,2a12,/,5x,27("-"))','j','scale','period'
        do j = 1, jtot
            scale(j) = s0*(2._rk**(dble(j - 1)*dj))
            call wave_function(npad, dt, mother, param, scale(j), kwave, period1, coi1, daughter)
            period(j) = period1
!**    multiply the daughter by the time-series fft
            do k = 1, npad
                daughter(k) = daughter(k)*yfft(k)
            end do
!**    inverse fft [eqn(4)]
            call cfftb(npad, daughter, wsave)
!**    store the wavelet transform, discard zero-padding at end
            do i = 1, n
                wave(i, j) = daughter(i)
            end do
!        print '(i8,2f12.3)',j,scale(j),period(j)
        end do
!**----- end loop

!** construct the cone of influence
        do i = 1, (n + 1)/2
            coi(i) = coi1*dt*(dble(i) - 1._rk)
            coi(n - i + 1) = coi(i)
        end do

        return
    end subroutine wavelet

!****************************************************************************
! wave_function: computes the daughter wavelets for a particular
!              wavelet function, with appropriate parameters.
!
!
! inputs:
!
!  nk [int] = the number of points in "kwave"
!
!  dt [double precision] = amount of time between each y value, i.e. the sampling time.
!
!  mother [int] = an integer giving the mother wavelet to use.
!                   0='morlet'
!                   1='paul'
!                   2='dog' (derivative of gaussian)
!
!  param [double precision] = mother wavelet parameter. if <0 then default is used.
!            for 'morlet' this is k0 (wavenumber), default is 6.
!            for 'paul' this is m (order), default is 4.
!            for 'dog' this is m (m-th derivative), default is 2.
!
!  scale1 [double precision] = the wavelet scale used to construct the daughter.
!
!  kwave [double precision(n)] = vector of wavenumbers, used to construct daughter.
!
!
! outputs:
!
!  period1 [double precision] = the "fourier" period (in time units) that corresponds
!            to "scale1".
!
!  coi1 [double precision] = the e-folding factor used for the cone of influence.
!
!  daughter [dcmplx(nk)] = real & imaginary parts of the wavelet function
!                    at "scale1" and "kwave".
!
!
! requires:   factorial, chisqr
!
!
! reference: tables 1 & 2 in
!            torrence, c. and g. p. compo, 1998: a practical guide to
!            wavelet analysis. <i>bull. amer. meteor. soc.</i>, 79, 61-78.
!
!  copyright (c) 1998, christopher torrence and gilbert p. compo
! this software may be used, copied, or redistributed as long as it is not
! sold and this copyright notice is reproduced on each copy made.  this
! routine is provided as is without any express or implied warranties
! whatsoever.

    subroutine wave_function(nk, dt, mother, param, scale1, kwave, period1, coi1, daughter)
        implicit none

        integer nk, mother
        real(rk) :: dt, kwave(nk), param, scale1, period1, coi1
        double complex daughter(nk), norm

        real(rk) :: expnt, sk, pi, fourier_factor
        integer k, m, factorial
        real(rk) :: gamma

        pi = 4._rk*atan(1._rk)

        if (mother == 0) then
!*******************************************   morlet wavelet
            if (param < 0) param = 6._rk
            norm = sqrt(2._rk*pi*scale1/dt)*(pi**(-0.25_rk))
            do k = 1, nk/2 + 1
                expnt = -0.5_rk*(scale1*kwave(k) - param)**2
                daughter(k) = dcmplx(norm*exp(expnt))
            end do
            do k = nk/2 + 2, nk
                daughter(k) = dcmplx(0._rk)
            end do
            fourier_factor = (4._rk*pi)/(param + sqrt(2._rk + param**2))
            period1 = scale1*fourier_factor
            coi1 = fourier_factor/sqrt(2._rk)
        else if (mother == 1) then
!*******************************************   paul wavelet
            if (param < 0) param = 4._rk
            m = int(param)
            norm = sqrt(2._rk*pi*scale1/dt)*(2._rk**m/sqrt(dble(m*factorial(2*m - 1))))
            do k = 1, nk/2 + 1
                expnt = -scale1*kwave(k)
                daughter(k) = dcmplx(norm*(scale1*kwave(k))**m*exp(expnt))
            end do
            do k = nk/2 + 2, nk
                daughter(k) = dcmplx(0._rk)
            end do
            fourier_factor = (4._rk*pi)/(2._rk*m + 1._rk)
            period1 = scale1*fourier_factor
            coi1 = fourier_factor*sqrt(2._rk)
        else if (mother == 2) then
!*******************************************   dog wavelet
            if (param < 0) param = 2._rk
            m = int(param)
            norm = sqrt(2._rk*pi*scale1/dt)*sqrt(1._rk/gamma(m + 0.5_rk))
            norm = -norm*(dcmplx(0._rk, 1._rk)**m)
            do k = 1, nk
                sk = scale1*kwave(k)
                daughter(k) = norm*(sk**m)*exp(-0.5_rk*sk**2)
            end do
            fourier_factor = 2._rk*pi*sqrt(2._rk/(2._rk*m + 1._rk))
            period1 = scale1*fourier_factor
            coi1 = fourier_factor/sqrt(2._rk)
        else
            stop
        end if
        return
    end subroutine wave_function

!****************************************************************************
! factorial: compute the factorial (n!) of an integer n
!  copyright (c) 1998, christopher torrence
    function factorial(n)
        implicit none
        integer factorial, n, i

        factorial = 1
        do i = 1, n
            factorial = factorial*i
        end do
    end function factorial

!****************************************************************************
! wave_signif: computes the significance levels for a wavelet transform.
!
!
! inputs:
!
!  isigtest [int] = 0, 1, or 2.
!
!          if 0, then just do a regular chi-square test,
!              i.e. eqn (18) from torrence & compo.
!          if 1, then do a "time-average" test, i.e. eqn (23).
!              in this case, dof(j) should be set to na, the number
!              of local wavelet spectra that were averaged together
!              at each scale. for the global wavelet spectrum,
!              this would be dof(j)=n-scale(j),
!              where n is the number of points in your time series.
!          if 2, then do a "scale-average" test, i.e. eqns (25)-(28).
!              in this case, "dof(1)" and "dof(2)" should be set to the
!              smallest (s1) and largest (s2) scales that were averaged
!              together, respectively.
!              e.g. if you scale-averaged scales between 2 and 8,
!                   then dof(1)=2.0 and dof(2)=8.0
!
!
!  n [int] = the number of points in "y".
!
!  y [double precision] = the time series of length "n".
!
!  dt [double precision] = amount of time between each y value, i.e. the sampling time.
!
!  mother [int] = an integer giving the mother wavelet to use.
!                   0='morlet'
!                   1='paul'
!                   2='dog' (derivative of gaussian)
!
!  param [double precision] = mother wavelet parameter.
!
!  s0 [double precision] = the smallest scale of the wavelet.
!
!  dj [double precision] = the spacing between discrete scales.
!
!  jtot [int] = the # of scales.
!
!  scale [double precision(jtot)] = the wavelet scales that were used.
!
!  period [double precision(jtot)] = the "fourier" periods corresponding to "scale".
!
!  lag1 [double precision] = lag 1 autocorrelation, used for signif levels.
!              default is 0.0, which corresponds to white-noise.
!
!  siglvl [double precision] = significance level to use. default is 0.05 (the "5%" level)
!
!  dof [double precision(jtot)] = degrees-of-freedom for signif test.
!     if sigtest=0, then the input dof is ignored.
!     if sigtest=1, then dof(j) = na, the number of times averaged together.
!     if sigtest=2, then dof(1)=s1, dof(2)=s2, the range of scales averaged.
!
!
! outputs:
!
!  dof [double precision(jtot)] = degrees-of-freedom that were actually used.
!     if sigtest=0, then dof(j) = 2 (or 1 for the 'dog')
!     if sigtest=1, then dof(j) = degrees-of-freedom versus scale.
!     if sigtest=2, then dof(1)=degrees-of-freedom, dof(2...jtot)=0.0
!
!  fft_theor [double precision(jtot)] = theoretical red-noise spectrum vs scale.
!     if sigtest=2, then fft_theor(1) = the average spectrum from s1-->s2
!                   fft_theor(2...jtot) = 0.0
!
!  signif [double precision(jtot)] = significance levels vs scale.
!     if sigtest=2, then signif(1) = the significance level
!                   signif(2...jtot) = 0.0
!
!  ymean [double precision] = the mean of the time series.
!
!  variance [double precision] = the variance of the time series.
!
!  cdelta [double precision] = the constant "cdelta" for the mother wavelet (table 2).
!
!  psi0[double precision] = the constant 'psi(0)' for the mother wavelet (table 2)
!
! requires:   chisqr
!
!  copyright (c) 1998, christopher torrence and gilbert p. compo
! this software may be used, copied, or redistributed as long as it is not
! sold and this copyright notice is reproduced on each copy made.  this
! routine is provided as is without any express or implied warranties
! whatsoever.

    subroutine wave_signif(isigtest, n, y, dt, mother, param, dj, jtot, scale, period, &
                           lag1, siglvl, dof, fft_theor, signif, ymean, variance, cdelta, psi0)
        implicit none

        integer isigtest, n, mother, jtot
        real(rk) :: y(n), dt, param, dj, scale(jtot), period(jtot)
        real(rk) :: lag1, siglvl, dof(jtot), fft_theor(jtot), signif(jtot)
        real(rk) :: ymean, variance, cdelta, psi0

        integer i, j, m, status, javg1, javg2, navg
        real(rk) :: pi, freq1, dofmin, gammafac, dj0, savg, smid
        real(rk) :: fft_theor1
        real(rk) :: chisqr, p, q, bound

        pi = acos(-1._rk)

        if (siglvl <= 0.) siglvl = 0.05_rk
        if (lag1 <= 0._rk) lag1 = 0._rk

        cdelta = -1._rk
        gammafac = -1._rk
        dj0 = -1._rk
        psi0 = -1._rk

        if (mother == 0) then
!*******************************************   morlet wavelet
            dofmin = 2._rk
            if (param == 6._rk) then
                cdelta = 0.776_rk
                gammafac = 2.32_rk
                dj0 = 0.60_rk
                psi0 = pi**(-0.25_rk)
            end if
        else if (mother == 1) then
!*******************************************   paul wavelet
            m = int(param)
            dofmin = 2._rk
            if (m == 4) then
                cdelta = 1.132_rk
                gammafac = 1.17_rk
                dj0 = 1.5_rk
                psi0 = 1.079_rk
            end if
        else if (mother == 2) then
!*******************************************   dog wavelet
            m = int(param)
            dofmin = 1._rk
            if (m == 2) then
                cdelta = 3.541_rk
                gammafac = 1.43_rk
                dj0 = 1.4_rk
                psi0 = 0.867_rk
            else if (m == 6) then
                cdelta = 1.966_rk
                gammafac = 1.37_rk
                dj0 = 0.97_rk
                psi0 = 0.884_rk
            end if
        else
            stop
        end if

!** find the time-series variance
        ymean = 0._rk
        do i = 1, n
            ymean = ymean + y(i)
        end do
        ymean = ymean/n
        variance = 0._rk
        do i = 1, n
            variance = variance + (y(i) - ymean)**2
        end do
        variance = variance/(dble(n))
!** construct theoretical red(white)-noise power spectrum [eqn(16)]
! - 1._rk)
        do j = 1, jtot
            freq1 = dt/period(j)
            fft_theor(j) = variance*(1._rk - lag1**2)/(1._rk - 2._rk*lag1*cos(freq1*2._rk*pi) + lag1**2)
        end do

        q = dble(siglvl)
        p = 1_rk - q

        if (isigtest == 0) then
!*******************************************   no smoothing, dof=dofmin
!   see eqn(18)
            do j = 1, jtot
                dof(j) = dofmin
                call cdfchi(2, p, q, chisqr, dble(dofmin), status, bound)
                signif(j) = fft_theor(j)*chisqr/dofmin
            end do
        else if (isigtest == 1) then
!***********************************   time-averaged, dof depend on scale
            if (gammafac <= 0._rk) then
                print *, '**wave_signif: "gammafac" undefined for this wavelet'
                return
            end if
            do j = 1, jtot
                if (dof(j) < 1.) dof(j) = 1._rk
!   see eqn(23)
                dof(j) = dofmin*sqrt(1._rk + (dof(j)*dt/gammafac/scale(j))**2)
                if (dof(j) < dofmin) dof(j) = dofmin
                call cdfchi(2, p, q, chisqr, dble(dof(j)), status, bound)
                signif(j) = fft_theor(j)*chisqr/dof(j)
            end do
        else if (isigtest == 2) then
!***********************************   scale-averaged, dof depend on scale
            if (cdelta <= 0.) then
                print *, '**wave_signif: "cdelta" and "dj0" '//'undefined for this wavelet'
                return
            end if
            javg1 = 0
            javg2 = 0
            do j = 1, jtot
                if ((scale(j) >= dof(1)) .and. (javg1 == 0)) javg1 = j
                if (scale(j) <= dof(2)) javg2 = j
            end do
            if ((javg1 == 0) .or. (javg2 == 0) .or. (javg1 > javg2)) then
                print *, '**wave_signif: scales in "dof(1)" & "dof(2)" '//'are out of range.'
                return
            end if
            navg = javg2 - javg1 + 1
!   see eqn(25)
            savg = 0._rk
            do j = javg1, javg2
                savg = savg + 1._rk/scale(j)
            end do
            savg = 1._rk/savg
!   see eqn(27)
            fft_theor1 = 0._rk
            do j = javg1, javg2
                fft_theor1 = fft_theor1 + fft_theor(j)/scale(j)
            end do
            fft_theor(1) = savg*fft_theor1
!   see eqn(28)
            smid = exp(0.5_rk*(log(scale(javg1)) + log(scale(javg2))))
            dof(1) = (dofmin*navg*savg/smid)*sqrt(1 + (navg*dj/dj0)**2)
!   see eqn(26)
            call cdfchi(2, p, q, chisqr, dble(dof(1)), status, bound)
            signif(1) = (dj*dt/cdelta/savg)*fft_theor(1)*chisqr/dof(1)
            do j = 2, jtot
                dof(j) = 0._rk
                fft_theor(j) = 0._rk
                signif(j) = 0._rk
            end do
        else
            stop
        end if

    end subroutine wave_signif

end module wavepack
