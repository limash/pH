module ph
!adopted from Munhoven2013

    use fabm_types, only: rk

    implicit none
    private

    public ph_solver

contains

    subroutine ph_solver(alktot, dictot, bortot, &
            po4tot, siltot, nh3tot, hstot, so4tot, flutot, &
            kc1, kc2, kb, kw)

        real(rk), intent(in):: alktot !total alkalinity
        real(rk), intent(in):: dictot !total dissolved
                                      !inorganic carbon
        real(rk), intent(in):: bortot !total boron
        real(rk), intent(in):: po4tot !total po4
        real(rk), intent(in):: siltot !total silicate
        real(rk), intent(in):: nh3tot !total nh3
        real(rk), intent(in):: hstot  !total hs
        real(rk), intent(in):: so4tot !total so4
        real(rk), intent(in):: flutot !total fluor
        !dissociation constants
        real(rk), intent(in):: kc1 !1st of carbonic acid
        real(rk), intent(in):: kc2 !2nd of carbonic acid
        real(rk), intent(in):: kb  !boric acid
        real(rk), intent(in):: kw  !water
        !auxiliary
        real(rk), intent(in):: iter, absmin

        !discriminant of the main equation R([H+])=0
        real(rk):: delta
        ![H+] starting value and etc.
        real(rk):: initial_h, h, h_prev
        !bracketing bounds for alk
        real(rk):: alkinf, alksup
        !bracketing bounds for [H+]
        real(rk):: h_min, h_max

        !calculate initial [H+] value
        call initial_h_do(alktot, dictot, bortot, &
            kc1, kc2, kb, initial_h)

        !calculate bracketing bounds for alk
        alkinf = -po4tot-so4tot-flutot
        alksup = dictot+dictot+bortot+ &
                 po4tot+po4tot+siltot+ &
                 nh3tot+hstot
 
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

end module ph
