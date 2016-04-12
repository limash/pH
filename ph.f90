module ph
!adopted from Munhoven2013

    use fabm_types, only: rk

    implicit none
    private

    public solver

contains

    subroutine solver(alktot, dictot, bortot, &
            kc1, kc2, kb)

        real(rk), intent(in):: alktot !total alkalinity
        real(rk), intent(in):: dictot !total dissolved
                                      !inorganic carbon
        real(rk), intent(in):: bortot !total boron
        !dissociation constants
        real(rk), intent(in):: kc1 !1st of carbonic acid
        real(rk), intent(in):: kc2 !2nd of carbonic acid
        real(rk), intent(in):: kb  !boric acid

        ![H+] starting value
        real(rk):: initial_h

        !calculates initial [H+] value
        call initial_h_do(alktot, dictot, bortot, &
            kc1, kc2, kb, initial_h)

    end subroutine solver

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
            !devided by 4
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
