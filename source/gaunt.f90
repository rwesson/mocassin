    ! this subroutine computes the gaunt factors for any charge
    ! it generates thermally averaged free-free non-relativistic gaunt
    ! factor for a hydrogenic ion of charge z, with a maximum relative
    ! error of 0.007, (rms fitting error = 0.001) for tmps and freqs in
    ! intervals:
    !          10^-4 <= U <= 10^1.5
    !          10^-3 <= Gams <= 10^3~
    ! where U = h*nu/k*T and gams = Z^2 * Ryd / k*T. To obtain the stated
    ! accuracy the full number of significant figures must be retained.
    !
    ! this subroutine uses a two-dimensional chebyshev epansion computed 
    ! from expressions given bu Karzas and Latter (ApJSuppl., V.6, P.167, 
    ! 1961) augmented by various limiting forms of energy-specific gaunt-
    ! factors.
    ! D.G.Hummer, Jila, May 1987. ApJ 327, 477
    ! modified with correct limits, J Ferguson, July 94
    ! modified for MOCASSIN by B Ercolano (Ercolano et al 2003, MNRAS 340, 1136) 
    subroutine getGauntFF(z, log10Te, xlf, g, iflag)
        implicit none

        integer, intent(out) :: iflag                ! explanation given in each area        

        real, intent(in) :: log10Te                  ! log10(Te)
        real, intent(in) ::  z                       ! log10 of nuclear charge
        real, dimension(:), intent(in) :: xlf        ! array of log10(h*nu/Ryd)
        real, dimension(size(xlf)), intent(out) :: g ! array of values of g

        ! local variables 
        integer ::  i, ir, j

        real, dimension(11) :: b
        real, dimension(8) :: c
        real :: con
        real(kind=8), dimension(88) :: dd
        real(kind=8), dimension(8, 11) :: d
        real :: gamma2                               ! gamma^2 (gams = Z^2 * Ryd / k*T)
        real :: slope                                ! slope
        real :: txg                                  ! hummer variable related to gamma^2
        real :: txu                                  ! hummer variable related to U
        real :: u                                    ! U = h*nu/k*T
        real :: xlrkt                                ! log(ryd/kt)

        dd = (/8.986940175e+00, -4.009515855e+00,  8.808871266e-01,& 
&          2.640245111e-02, -4.580645915e-02, -3.568055702e-03,&
&      	  2.827798067e-03,  3.365860195e-04, -8.006936989e-01,&
&      	  9.466021705e-01,  9.043402532e-02, -9.608451450e-02,&
&         -1.885629865e-02,  1.050313890e-02,  2.800889961e-03,&
&         -1.078209202e-03, -3.781305103e-01,  1.102726332e-01,&
&         -1.543619180e-02,  8.310561114e-03,  2.179620525e-02,&
&          4.259726289e-03, -4.181588794e-03, -1.770208330e-03,&
&          1.877213132e-02, -1.004885705e-01, -5.483366378e-02,& 
&         -4.520154409e-03,  8.366530426e-03,  3.700273930e-03,&
&          6.889320423e-04,  9.460313195e-05,  7.300158392e-02,&
&          3.576785497e-03, -4.545307025e-03, -1.017965604e-02,&
&         -9.530211924e-03, -3.450186162e-03,  1.040482914e-03,&
&          1.407073544e-03, -1.744671550e-03,  2.864013856e-02,&
&          1.903394837e-02,  7.091074494e-03, -9.668371391e-04,&
&         -2.999107465e-03, -1.820642230e-03, -3.874082085e-04,&
&         -1.707268366e-02, -4.694254776e-03,  1.311691517e-03,&
&          5.316703136e-03,  5.178193095e-03,  2.451228935e-03,&
&         -2.277321615e-05, -8.182359057e-04,  2.567331664e-04,&
&         -9.155339970e-03, -6.997479192e-03, -3.571518641e-03,&
&         -2.096101038e-04,  1.553822487e-03,  1.509584686e-03,&
&          6.212627837e-04,  4.098322531e-03,  1.635218463e-03,&
&         -5.918883504e-04, -2.333091048e-03, -2.484138313e-03,&
&         -1.359996060e-03, -5.371426147e-05,  5.553549563e-04,& 
&          3.837562402e-05,  2.938325230e-03,  2.393747064e-03,&
&          1.328839809e-03,  9.135013312e-05, -7.137252303e-04,&
&         -7.656848158e-04, -3.504683798e-04, -8.491991820e-04,&
&         -3.615327726e-04,  3.148015257e-04,  8.909207650e-04,&
&          9.869737522e-04,  6.134671184e-04,  1.068883394e-04,&
&         -2.046080100e-04/)
      
        d = reshape( dd, (/8, 11/) )

        ! compute temperature dependent coeffcients for U expansion

        ! xlrxt is log(ryd/kt), code note valid for Te > 10^8
        ! xlrxt is log(ryd/kt), code note valid for Te < 10^2.5
        if ( log10Te < 2.5 ) then
            xlrkt = 5.1983649 - 2.5
        else if ( log10Te > 8. ) then
            xlrkt = 5.1983649 - 8.
        else
            xlrkt = 5.1983649 - log10Te
        end if

        ! set txg
        txg = 0.66666667*(2.0*z+xlrkt)
        gamma2 = 10**(txg*1.5)

        con = 0.72727273*xlrkt+0.90909091
        do j=1,8
            ir = 9
            b(11) = d(j,11)
            b(10) = txg*b(11)+d(j,10)
            do i=1,9
                b(ir) = txg*b(ir+1)-b(ir+2)+d(j,ir)
                ir = ir-1
            end do
            c(j) = 0.25*(b(1)-b(3))
        end do
        
        ! sum U expansion
        ! loop through energy at fixed temperature
        do i = 1, size(xlf)
            txu = 0.72727273*xlf(i)+con
            u = 10**((txu - .90909091)/.72727273)
            ! criteria set by hummer limits. it is a wedge from
            ! log(hnu),log(T)) =
            ! (-5,2.5) to (-5,4) to (-1,8) to (4,8) to (-1.5,2.5).
            ! these limits correspond to the gamma^2 and U limits 
            ! given above
            if(abs(txu)<=2.0) then
                ir = 6
                b(8) = c(8)
                b(7) = txu*b(8)+c(7)
                do j=1,6
                    b(ir) = txu*b(ir+1)-b(ir+2)+c(ir)
                    ir = ir-1
                end do
                g(i) = b(1)-b(3)
                if( (log10Te>=2.5) .and. (log10Te<=8.0)) iflag = 0
                ! On the bottom side of the hummer box,u<-4 and gamma2>.33
            else if( (log10(u)<-4.0) .and. (gamma2<0.3) ) then
               g(i) = 0.551329 * log(2.24593/u)
              if( (log10Te>=2.5) .and. (log10Te<=8.0) ) iflag = 2
            ! On the bottom side of the box,u<-4 and gamma2<.33
            else if( (log10(u)<-4.0) .and. (gamma2>=0.3) ) then
                g(i) = 0.551329 * alog( 0.944931/u/sqrt(gamma2) )
                if( (log10Te>=2.5) .and.( log10Te<=8.0) ) iflag = 3
                ! Now on the bottom side of the box
                else if (txu>2.0) then
                ! Top of the box high T first
                if (log10(gamma2)<-3.) then
                    g(i) =  sqrt(0.9549297/u)
                    if(log10Te>=2.5.and.log10Te<=8.0) iflag = 4
                    ! Must interpolate between two asymptotes
                else if(log10(gamma2)<0.0) then
                    slope = ( sqrt(12*gamma2/u) - sqrt(0.9549297/u) ) / 3.0
                    g(i) =  sqrt(12*gamma2/u) + slope * log10(gamma2)
                    if((log10Te>=2.5) .and. (log10Te<=8.0)) iflag = 7
                else
                ! The top side with wedge of 1.05 where region 6 fails
                    g(i) = sqrt(12. * gamma2 / u)
                    g(i) = min(1.05,g(i))
                    if( (g(i)==1.05) .and. (log10Te>=2.5) .and. (log10Te<=8.0)) &
                         & iflag = 5
                    if( (g(i)<1.05) .and. (log10Te>=2.5) .and. (log10Te<=8.0)) &
                         & iflag = 6
                end if
            end if
            u = 0.0
        end do
    end subroutine getGauntFF
