!****************************************************************************
! WAVETEST: Example Fortran program for WAVELET, using NINO3 SST dataset
!
! COMPILE:   f77 chisqr.f cfftpack.f wavelet.f wavetest.f
!
! See "http://paos.colorado.edu/research/wavelets/"
!
!  Copyright (C) 1998, Christopher Torrence and Gilbert P. Compo
! This software may be used, copied, or redistributed as long as it is not
! sold and this copyright notice is reproduced on each copy made.  This
! routine is provided as is without any express or implied warranties
! whatsoever.
!
! Modified: November 1999 by Arjan van Dijk to include IMPLICIT NONE and
!           to convert all routines to DOUBLE precision.
!****************************************************************************

program wavetest

    use wavepack, only: rk, wavelet, wave_signif
    implicit none

    integer n, subscale, jtot
    real(rk) :: dt, s0, dj

    ! these parameters depend on the particular time series
    parameter(n=504, dt=0.25_rk, s0=dt)
    parameter(subscale=4)
    parameter(dj=1._rk/subscale, jtot=11*subscale)
    ! Note: for accurate reconstruction and wavelet-derived variance
    !     do not pad with zeroes, set s0=dt (for Paul set s0=dt/4), and use
    !     a large "jtot" (even though the extra scales will be within
    !     the cone of influence).
    !     For plotting purposes, it is only necessary to use
    !     s0=2dt (for Morlet) and "jtot" from Eqn(10) Torrence&Compo(1998).

    integer mother, ibase2, npad
    real(rk) :: sst(n), recon_sst(n), param, pi
    real(rk) :: scale(jtot), period(jtot), coi(n)
    double complex wave(n, jtot)

    integer :: i, j, isigtest, javg1, javg2
    real(rk) :: lag1, siglvl, dof(jtot)
    real(rk) :: fft_theor(jtot), signif(jtot), ymean, variance
    real(rk) :: recon_mean, recon_vari
    real(rk) :: cdelta, psi0
    real(rk) :: global_ws(jtot), global_signif(jtot)
    real(rk) :: savg_dof(jtot), savg_signif(jtot), sstenso(n)

    pi = 4._rk*atan(1._rk)
    ibase2 = nint(log(dble(n))/log(2._rk)) + 1
    npad = int(2._rk**ibase2)
    !      npad = n  ! this is for no padding with zeroes

    !*************************************************** Wavelet transform

    !** let the WAVELET subroutine choose the defaults for these:
    mother = 0
    param = 6._rk

    !** read in the NINO3 SST data
    open (unit=11, file='data/sst_nino3.dat', status='old')
    read (11, *) sst
    close (11)
    print '(/,"sst(1)=",f6.2,"  sst(n) = ",f6.2,/)', sst(1), sst(n)

    !** get the wavelet transform
    call wavelet(n, sst, dt, mother, param, s0, dj, jtot, npad, wave, scale, period, coi)

    !*************************************************** Significance testing

    !** local significance test
    isigtest = 0
    lag1 = 0.72_rk
    siglvl = 0.05_rk
    call wave_signif(isigtest, n, sst, dt, mother, param, dj, jtot, scale, period, &
                     lag1, siglvl, dof, fft_theor, signif, ymean, variance, cdelta, psi0)

    !** global wavelet spectrum & significance test
    isigtest = 1
    lag1 = 0.72_rk
    siglvl = 0.05_rk
    do j = 1, jtot
        do i = 1, n
            global_ws(j) = global_ws(j) + abs(wave(i, j))**2
        end do
        global_ws(j) = global_ws(j)/n
        dof(j) = n - scale(j)
    end do

    call wave_signif(isigtest, n, sst, dt, mother, param, dj, jtot, scale, period, &
                     lag1, siglvl, dof, fft_theor, global_signif, ymean, variance, cdelta, psi0)

    !** scale-average time series & significance test
    isigtest = 2
    lag1 = 0.72_rk
    siglvl = 0.05_rk
    !    scale average between 2 and 7.9 years
    savg_dof(1) = 2.0_rk
    savg_dof(2) = 7.9_rk
    !    find the "j"-values that correspond to savg_dof(1) & savg_dof(2)
    javg1 = 0
    javg2 = 0
    do j = 1, jtot
        if ((scale(j) >= savg_dof(1)) .and. (javg1 == 0)) javg1 = j
        if (scale(j) <= savg_dof(2)) javg2 = j
    end do
    !   call wave_signif first, to get the value of "Cdelta"
    call wave_signif(isigtest, n, sst, dt, mother, param, dj, jtot, scale, period, &
                     lag1, siglvl, savg_dof, fft_theor, savg_signif, ymean, variance, cdelta, psi0)
    !   construct the scale-averaged time series [Eqn(24)]
    do i = 1, n
        sstenso(i) = 0._rk
        do j = javg1, javg2
            sstenso(i) = sstenso(i) + (abs(wave(i, j))**2)/scale(j)
        end do
        sstenso(i) = dj*dt*sstenso(i)/cdelta
    end do

    !************************************************************* print results
    print *, ' n=', n
    print *, ' dt=', dt
    print *, ' mother=', mother
    print *, ' param=', param
    print *, ' s0=', s0
    print *, ' dj=', dj
    print *, ' jtot=', jtot
    print *, ' npad=', npad
    print '(/,"let w = wave(n/2,j)",/)'
    print '(a4,7a10)', 'j', 'scale', 'period', 'abs(w)^2', 'phase(w)', '5%signif', 'global', 'gws5%sig'
    print '(i4,7f10.3)', (j, scale(j), period(j), abs(wave(n/2, j))**2, &
                          atan2(dimag(wave(n/2, j)), dble(wave(n/2, j)))*180._rk/pi, &
                          signif(j), global_ws(j), global_signif(j), j=1, jtot)
    print '(/,a,f10.3)', ' scale-average degrees of freedom = ', savg_dof(1)
    print '(a,f10.3,/)', ' scale-avg 5% significance level  = ', savg_signif(1)

    !************************************************************ Reconstruction

    !** construct the wavelet derived variance (Parseval's theorem)  [Eqn(14)]
    !   Cdelta & psi0 are returned from WAVE_SIGNIF
    recon_vari = 0._rk
    do i = 1, n
        do j = 1, jtot
            recon_vari = recon_vari + (abs(wave(i, j))**2)/scale(j)
        end do
    end do
    recon_vari = dj*dt*recon_vari/(cdelta*n)
    print '(a,f14.5)', ' reconstructed variance=', recon_vari
    print '(a,f14.5)', ' original variance   =', variance
    print '(a,f14.5,a,/)', ' ratio = ', recon_vari/variance, ' (this is low due to padding with zeroes)'

    !** reconstruct the time series [Eqn(11)]
    !   check mean and RMS difference of reconstructed time series
    recon_mean = 0._rk
    recon_vari = 0._rk
    do i = 1, n
        recon_sst(i) = 0._rk
        do j = 1, jtot
            recon_sst(i) = recon_sst(i) + (dble(wave(i, j)))/sqrt(scale(j))
        end do
        recon_sst(i) = dj*sqrt(dt)*recon_sst(i)/(cdelta*psi0)
        recon_vari = recon_vari + (sst(i) - ymean - recon_sst(i))**2
        recon_mean = recon_mean + recon_sst(i)
    end do
    recon_mean = recon_mean/n
    recon_vari = sqrt(recon_vari/n)

    print '(a,f14.6)', ' reconstructed mean=', recon_mean
    print '(a,f14.6)', ' original mean   =', ymean
    print '(a,f14.6,/)', ' root-mean-square difference of time series=', recon_vari

end program wavetest
