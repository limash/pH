module ph
!adopted from Munhoven2013

    use fabm_types, only: rk

    implicit none
    private

    public ph_solver

contains

    subroutine ph_solver(alktot, dictot, bortot, &
            po4tot, siltot, nh4tot, hstot, so4tot, flutot, &
            kc1, kc2, kb, kp1, kp2, kp3, ks, kn, khs, kw)

        real(rk), intent(in):: alktot !total alkalinity
        real(rk), intent(in):: dictot !total dissolved
                                      !inorganic carbon
        real(rk), intent(in):: bortot !total boron
        real(rk), intent(in):: po4tot !total po4
        real(rk), intent(in):: siltot !total silicate
        real(rk), intent(in):: nh4tot !total nh3
        real(rk), intent(in):: hstot  !total hs
        real(rk), intent(in):: so4tot !total so4
        real(rk), intent(in):: flutot !total fluor
        !dissociation constants
        real(rk), intent(in):: kc1 !1st of carbonic acid
        real(rk), intent(in):: kc2 !2nd of carbonic acid
        real(rk), intent(in):: kb  !boric acid
        real(rk), intent(in):: kp1, kp2, kp3 !phosphoric acid 
        real(rk), intent(in):: ks  !silicic acid
        real(rk), intent(in):: kn  !Ammonia
        real(rk), intent(in):: khs !Hydrogen sulfide
        real(rk), intent(in):: kw  !water

        !discriminant of the main equation R([H+])=0
        !r - value of main eq., dr - its derivative
        real(rk):: delta, r, dr
        ![H+] starting value and etc.
        real(rk):: initial_h, h, h_prev
        !bracketing bounds for alk
        real(rk):: alkinf, alksup
        !bracketing bounds for [H+]
        real(rk):: h_min, h_max
        !auxiliary
        real(rk):: iter, absmin

        !calculate initial [H+] value
        call initial_h_do(alktot, dictot, bortot, &
            kc1, kc2, kb, initial_h)

        !calculate bracketing bounds for alk
        alkinf = -po4tot-so4tot-flutot
        alksup = dictot+dictot+bortot+ &
                 po4tot+po4tot+siltot+ &
                 nh4tot+hstot
 
        !calculate discriminant for lower bound
        delta = (alktot-alkinf)**2+4._rk*kw
        !calculate lower bound
        if (alktot >= alkinf) then
            h_min = 2._rk*kw/(alktot-alkinf+sqrt(delta))
        else
            h_min = (-(alktot-alkinf)+sqrt(delta))/2._rk
        end if

        !calculate discriminant for upper bound
        delta = (alktot-alksup)**2+4._rk*kw
        !calculate upper bound
        if (alktot <= alksup) then
            h_max = (-(alktot-alksup)+sqrt(delta))/2._rk
        else
            h_max = 2._rk*kw/(alktot-alksup+sqrt(delta))
        end if

        !main algorithm
        h = max(min(h_max, initial_h), h_min)
        iter = 0
        absmin = huge(1._rk)
        do
            h_prev = h
            !return value (r) of main equation for given [H+]
            !and its derivative (dr)
            call r_calc(h, alktot, dictot, bortot, po4tot, &
                siltot, nh4tot, hstot, so4tot, flutot, &
                kc1, kc2, kb, kp1, kp2, kp3, ks, kn, khs, kw, &
                r, dr)
            !adapt bracketing interval
            if (r > 0._rk) then
                h_min = h_prev
            else if (r < 0._rk) then
                h_max = h_prev
            else
                !h is root
                exit
            end if

            iter = iter + 1

            !bisection method pH-Alk space
            if (abs(r) >= 0.5_rk*absmin) then
                h = sqrt(h_max * h_min)
                !required for convergence test
                h_fac = (h-h_prev)/h_prev

        end do

    end subroutine ph_solver

    subroutine initial_h_do(alktot, dictot, bortot, &
            kc1, kc2, kb, initial_h)
    !calculates initial value for [H+]

        real(rk), intent(in):: alktot
        real(rk), intent(in):: dictot
        real(rk), intent(in):: bortot
        real(rk), intent(in):: kc1
        real(rk), intent(in):: kc2
        real(rk), intent(in):: kb
        real(rk), intent(out):: initial_h

        real(rk):: dic_alk, bor_alk !auxiliary variables
        real(rk):: a0, a1, a2 !coefficients of the cubic polynomial
        real(rk):: d, d_sqrt !discriminant and square of it
        real(rk):: h_min !local minimum

        if (alktot <= 0._rk) then
            initial_h = 1.e-3_rk
        else if (alktot >= (2._rk*dictot+bortot)) then
            initial_h = 1.e-10_rk
        else
            dic_alk = dictot/alktot
            bor_alk = bortot/alktot

            !coefficients of the cubic polynomial equation
            ![H+]*3+a2[H+]*2+a1[H+]+a0=0
            !derived from alkalinity-pH equation
            !for given alkalinity, dictot, bortot
            a2 = kb*(1._rk-bor_alk)+kc1*(1._rk-dic_alk)
            a1 = kc1*kb*(1._rk-bor_alk-dic_alk) &
                +kc2*(1._rk-(dic_alk+dic_alk))
            a0 = kc2*kb*(1._rk-bor_alk-(dic_alk+dic_alk))

            !discriminant of the quadratic equation
            !(derivative of the cubic equation)
            !for the minimum close to the root
            d = a2*a2-3._rk*a1

            if (d > 0._rk) then
                d_sqrt = sqrt(d)
                if (a2 < 0) then
                    h_min = (-a2+d_sqrt)/3._rk
                else
                    h_min = -a1/(a2+d_sqrt)
                end if

                initial_h = h_min+sqrt(-(a0+h_min*(a1+h_min(a2+hmin))) &
                    /d_sqrt)
            else
                initial_h = 1.e-7_rk
            end if
        end if

    end subroutine initial_h_do

    subroutine r_calc(h, alktot, dictot, bortot, po4tot, &
            siltot, nh4tot, hstot, so4tot, flutot, &
            kc1, kc2, kb, kp1, kp2, kp3, ks, kn, khs, kw, &
            r, dr)
    !return value of main equation for given [H+]
    !and value of its derivative

        real(rk), intent(in):: h ![H+]
        real(rk), intent(in):: alktot !total alkalinity
        real(rk), intent(in):: dictot !total dissolved
                                      !inorganic carbon
        real(rk), intent(in):: bortot !total boron
        real(rk), intent(in):: po4tot !total po4
        real(rk), intent(in):: siltot !total silicate
        real(rk), intent(in):: nh4tot !total nh3
        real(rk), intent(in):: hstot  !total hs
        real(rk), intent(in):: so4tot !total so4
        real(rk), intent(in):: flutot !total fluor
        !dissociation constants
        real(rk), intent(in):: kc1 !1st of carbonic acid
        real(rk), intent(in):: kc2 !2nd of carbonic acid
        real(rk), intent(in):: kb  !boric acid
        real(rk), intent(in):: kp1, kp2, kp3 !phosphoric acid 
        real(rk), intent(in):: ks  !silicic acid
        real(rk), intent(in):: kn  !Ammonia
        real(rk), intent(in):: khs !Hydrogen sulfide
        real(rk), intent(in):: kw  !water
        !output variables
        real(rk), intent(out):: r, rd

        real(rk):: dic1, dic2, dic, dddic, ddic
        real(rk):: bor1, bor2, bor, ddbor, dbor
        real(rk):: po4_1, po4_2, po4, ddpo4, dpo4
        real(rk):: sil1, sil2, sil, ddsil, dsil
        real(rk):: nh4_1, nh4_2, nh4, ddnh4, dnh4
        real(rk):: h2s_1, h2s_2, h2s, ddh2s, dh2s
        real(rk):: wat

        !H2CO3 - HCO3 - CO3
        dic1 = 2._wp*kc2 + h*       kc1
        dic2 =       kc2 + h*(      kc1 + h)
        dic  = dictot * (dic1/dic2)
        !B(OH)3 - B(OH)4
        bor1 =       kb
        bor2 =       kb + h
        bor  = bortot * (bor1/bor2)
        !H3PO4 - H2PO4 - HPO4 - PO4
        po4_1 = 3._wp*kp3 + h*(2._wp*kp2 + h* kp1)
        po4_2 =       kp3 + h*(      kp2 + h*(kp1 + h))
        po4   = po4tot * (po4_1/po4_2 - 1._wp) ! Zero level of H3PO4 = 1
        !H4SiO4 - H3SiO4
        sil1 =       ks
        sil2 =       ks + h
        sil  = siltot * (sil1/sil2)
        !NH4 - NH3
        nh4_1 =       kn
        nh4_2 =       kn + h
        nh4   = nh4tot * (nh4_1/nh4_2)
        !H2S - HS
        h2s_1 =       khs
        h2s_2 =       khs + h
        h2s   = hstot * (h2s_1/h2s_2)
        !H2O - OH
        wat   = kw/h - h

        r = dic + bor + po4 + sil &
          + nh4 + h2s &
          + wat - alktot

        !H2CO3 - HCO3 - CO3
        dddic = kc1*kc2 + h*(4._wp*kc2+h*kc1)
        ddic  = -dictot*(dddic/dic2**2)
        !B(OH)3 - B(OH)4
        ddbor = kb
        dbor  = -bortot*(ddbor/bor2**2)
        !H3PO4 - H2PO4 - HPO4 - PO4
        ddpo4 = kp2*kp3 + h*(4._wp*kp1*kp3 &
              + h*(9._wp*kp3 + kp1*kp2     &
              + h*(4._wp*kp2               &
              + h*       kp1)))
        dpo4   = -po4tot * (ddpo4/po4_2**2)
        !H4SiO4 - H3SiO4
        ddsil = ks
        dsil  = -siltot * (ddsil/sil2**2)
        !NH4 - NH3
        ddnh4 = kn
        dnh4  = -nh4tot * (ddnh4/nh4_2**2)
        !H2S - HS
        ddh2s = khs
        dh2s  = -h2stot * (ddh2s/h2s_2**2)

        dr = ddic + dbor + dpo4 + dsil &
           + dnh4 + dh2s &
           - kw/h**2 - 1._rk

    end subroutine r_calc

end module ph
