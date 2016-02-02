!======MODULE CODE FOR test======

subroutine calcenergyforces(ANumList, Pos, BoxL, AID, CutSq, Const, CalcForce, CalcDUParam, CalcDWParam, CalcTerms, PEnergy, &
  & Virial, Force, DUParam, DDUParam, DWParam, DDWParam, Terms, NList, NConst, NAtom, Dim, NDParam, NDDParam, NTerm)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: NConst
    integer, intent(in) :: NAtom
    integer, intent(in) :: Dim
    integer, intent(in) :: NDParam
    integer, intent(in) :: NDDParam
    integer, intent(in) :: NTerm
    integer, dimension(NList), intent(in) :: ANumList
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Pos
    real(8), dimension(0:Dim-1), intent(in) :: BoxL
    integer, dimension(0:NAtom-1), intent(in) :: AID
    real(8), dimension(0:NTerm-1), intent(in) :: CutSq
    real(8), dimension(0:NConst-1), intent(in) :: Const
    logical, intent(in) :: CalcForce
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcTerms
    real(8), intent(out) :: PEnergy
    real(8), intent(out) :: Virial
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Force
!f2py intent(in, out) :: Force
    real(8), dimension(0:NDParam-1), intent(inout) :: DUParam
!f2py intent(in, out) :: DUParam
    real(8), dimension(0:NDDParam-1), intent(inout) :: DDUParam
!f2py intent(in, out) :: DDUParam
    real(8), dimension(0:NDParam-1), intent(inout) :: DWParam
!f2py intent(in, out) :: DWParam
    real(8), dimension(0:NDDParam-1), intent(inout) :: DDWParam
!f2py intent(in, out) :: DDWParam
    real(8), dimension(0:NTerm-1), intent(inout) :: Terms
!f2py intent(in, out) :: Terms
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: AIDi
    integer :: AIDj
    integer :: ListIndi
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8), dimension(0:Dim-1) :: Forcei
    real(8) :: ThisU
    real(8) :: ThisW
    real(8) :: ThisWsum
    real(8) :: ThisdU
    real(8), dimension(0:Dim-1) :: ThisForce
    real(8) :: ThisLocalRho
    real(8) :: ThisLocalRhoF
    real(8) :: LocalRho
    real(8) :: LocalRhoW
    real(8) :: SPt
    integer :: SPInd
    real(8) :: SPx
    real(8) :: SPdm1
    real(8) :: SPd0
    real(8) :: SPdp1
    real(8) :: SPdp2
    real(8) :: RhoSc_3
    real(8) :: val1
    real(8) :: LocalRho_3
    real(8) :: LocalRhoW_3
    real(8), dimension(0:NATOM-1,0:DIM-1) :: LocalRhoF_3

    !compute initial quantities
    PEnergy = 0.d0
    Virial = 0.d0
    if (CalcForce) then
        if (NList == 0) then
            Force = 0.d0
        else
            do ListIndi = 1, NList
                i = ANumList(ListIndi)
                Force(i,:) = 0.d0
            enddo
        endif
    endif
    if (CalcDUParam) then
        DUParam = 0.d0
        DDUParam = 0.d0
    endif
    if (CalcDWParam) then
        DWParam = 0.d0
        DDWParam = 0.d0
    endif
    if (CalcTerms) Terms = 0.d0

    !initialization-commands for potential <LD_BW>
    RhoSc_3 = Const(415)

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    !loop over i
    i = -1
    ListIndi = 0
    do while ((NList == 0 .and. i < NAtom-1) .or. (ListIndi < NList))
        if (NList == 0) then
            i = i + 1
        else
            ListIndi = ListIndi + 1
            i = ANumList(ListIndi)
        endif

        Posi = Pos(i,:)
        AIDi = AID(i)

        if (AIDi==0) then
            LocalRho_3 = 0.d0
            LocalRhoW_3 = 0.d0
            if (CalcForce) LocalRhoF_3 = 0.d0
        end if

        !loop over j
        do j = 0, NAtom-1
            if (i == j) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            !zero the virial sum
            ThisWsum = 0.d0

            if (i < j) then
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    if (dijsq < CutSq(0)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_BB>
                        SPx = DIJ * Const(1)
                        SPind = max(min(int(SPx), 24), 0)
                        SPt = SPx - real(SPind)
                        THISU = Const(2 + SPInd) + SPt * (Const(27 + SPInd) + SPt * (Const(52 + SPInd) + SPt * &
                          & Const(77 + SPInd)))
                        THISW = (Const(27 + SPInd) + SPt * (2.*Const(52 + SPInd) + SPt * 3.*Const(77 + SPInd))) * SPx
                        PEnergy = PEnergy + ThisU
                        Virial = Virial + ThisW
                        ThisWsum = ThisWsum + ThisW
                        if (CalcTerms) then
                            Terms(0) = Terms(0) + ThisU
                        endif
                        if (CalcDUParam) then
                            SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                            SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                            SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                            SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                            if (SPind == 0) then
                                !note: y(-1) = 2 y(0) - y(1) for linearity condition
                                DUParam(2 + SPInd) = DUParam(2 + SPInd) + SPd0 + 2.0d0 * SPdm1
                                DUParam(2 + SPInd + 1) = DUParam(2 + SPInd + 1) + SPdp1 - SPdm1
                                DUParam(2 + SPInd + 2) = DUParam(2 + SPInd + 2) + SPdp2
                            elseif (SPind == 24) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DUParam(2 + SPInd - 1) = DUParam(2 + SPInd - 1) + SPdm1
                                DUParam(2 + SPInd) = DUParam(2 + SPInd) + SPd0 - 0.5d0 * SPdp1 + SPdp2
                            elseif (SPind == 23) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DUParam(2 + SPInd - 1) = DUParam(2 + SPInd - 1) + SPdm1
                                DUParam(2 + SPInd) = DUParam(2 + SPInd) + SPd0
                                DUParam(2 + SPInd + 1) = DUParam(2 + SPInd + 1) + SPdp1 - 0.5d0 * SPdp2
                            else
                                DUParam(2 + SPInd - 1) = DUParam(2 + SPInd - 1) + SPdm1
                                DUParam(2 + SPInd) = DUParam(2 + SPInd) + SPd0
                                DUParam(2 + SPInd + 1) = DUParam(2 + SPInd + 1) + SPdp1
                                DUParam(2 + SPInd + 2) = DUParam(2 + SPInd + 2) + SPdp2
                            endif
                        endif
                        if (CalcDWParam) then
                            SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
                            SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
                            SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
                            SPdp2 = (SPt * SPt * 0.5d0) * SPx
                            if (SPind == 0) then
                                !note: y(-1) = 2 y(0) - y(1) for linearity condition
                                DWParam(2 + SPInd) = DWParam(2 + SPInd) + SPd0 + 2.0d0 * SPdm1
                                DWParam(2 + SPInd + 1) = DWParam(2 + SPInd + 1) + SPdp1 - SPdm1
                                DWParam(2 + SPInd + 2) = DWParam(2 + SPInd + 2) + SPdp2
                            elseif (SPind == 24) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DWParam(2 + SPInd - 1) = DWParam(2 + SPInd - 1) + SPdm1
                                DWParam(2 + SPInd) = DWParam(2 + SPInd) + SPd0 - 0.5d0 * SPdp1 + SPdp2
                            elseif (SPind == 23) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DWParam(2 + SPInd - 1) = DWParam(2 + SPInd - 1) + SPdm1
                                DWParam(2 + SPInd) = DWParam(2 + SPInd) + SPd0
                                DWParam(2 + SPInd + 1) = DWParam(2 + SPInd + 1) + SPdp1 - 0.5d0 * SPdp2
                            else
                                DWParam(2 + SPInd - 1) = DWParam(2 + SPInd - 1) + SPdm1
                                DWParam(2 + SPInd) = DWParam(2 + SPInd) + SPd0
                                DWParam(2 + SPInd + 1) = DWParam(2 + SPInd + 1) + SPdp1
                                DWParam(2 + SPInd + 2) = DWParam(2 + SPInd + 2) + SPdp2
                            endif
                        endif
                    end if
                end if
            end if

            if (i < j) then
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    if (dijsq < CutSq(1)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_WW>
                        SPx = DIJ * Const(103)
                        SPind = max(min(int(SPx), 24), 0)
                        SPt = SPx - real(SPind)
                        THISU = Const(104 + SPInd) + SPt * (Const(129 + SPInd) + SPt * (Const(154 + SPInd) + SPt * &
                          & Const(179 + SPInd)))
                        THISW = (Const(129 + SPInd) + SPt * (2.*Const(154 + SPInd) + SPt * 3.*Const(179 + SPInd))) * SPx
                        PEnergy = PEnergy + ThisU
                        Virial = Virial + ThisW
                        ThisWsum = ThisWsum + ThisW
                        if (CalcTerms) then
                            Terms(1) = Terms(1) + ThisU
                        endif
                        if (CalcDUParam) then
                            SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                            SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                            SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                            SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                            if (SPind == 0) then
                                !note: y(-1) = 2 y(0) - y(1) for linearity condition
                                DUParam(27 + SPInd) = DUParam(27 + SPInd) + SPd0 + 2.0d0 * SPdm1
                                DUParam(27 + SPInd + 1) = DUParam(27 + SPInd + 1) + SPdp1 - SPdm1
                                DUParam(27 + SPInd + 2) = DUParam(27 + SPInd + 2) + SPdp2
                            elseif (SPind == 24) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DUParam(27 + SPInd - 1) = DUParam(27 + SPInd - 1) + SPdm1
                                DUParam(27 + SPInd) = DUParam(27 + SPInd) + SPd0 - 0.5d0 * SPdp1 + SPdp2
                            elseif (SPind == 23) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DUParam(27 + SPInd - 1) = DUParam(27 + SPInd - 1) + SPdm1
                                DUParam(27 + SPInd) = DUParam(27 + SPInd) + SPd0
                                DUParam(27 + SPInd + 1) = DUParam(27 + SPInd + 1) + SPdp1 - 0.5d0 * SPdp2
                            else
                                DUParam(27 + SPInd - 1) = DUParam(27 + SPInd - 1) + SPdm1
                                DUParam(27 + SPInd) = DUParam(27 + SPInd) + SPd0
                                DUParam(27 + SPInd + 1) = DUParam(27 + SPInd + 1) + SPdp1
                                DUParam(27 + SPInd + 2) = DUParam(27 + SPInd + 2) + SPdp2
                            endif
                        endif
                        if (CalcDWParam) then
                            SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
                            SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
                            SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
                            SPdp2 = (SPt * SPt * 0.5d0) * SPx
                            if (SPind == 0) then
                                !note: y(-1) = 2 y(0) - y(1) for linearity condition
                                DWParam(27 + SPInd) = DWParam(27 + SPInd) + SPd0 + 2.0d0 * SPdm1
                                DWParam(27 + SPInd + 1) = DWParam(27 + SPInd + 1) + SPdp1 - SPdm1
                                DWParam(27 + SPInd + 2) = DWParam(27 + SPInd + 2) + SPdp2
                            elseif (SPind == 24) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DWParam(27 + SPInd - 1) = DWParam(27 + SPInd - 1) + SPdm1
                                DWParam(27 + SPInd) = DWParam(27 + SPInd) + SPd0 - 0.5d0 * SPdp1 + SPdp2
                            elseif (SPind == 23) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DWParam(27 + SPInd - 1) = DWParam(27 + SPInd - 1) + SPdm1
                                DWParam(27 + SPInd) = DWParam(27 + SPInd) + SPd0
                                DWParam(27 + SPInd + 1) = DWParam(27 + SPInd + 1) + SPdp1 - 0.5d0 * SPdp2
                            else
                                DWParam(27 + SPInd - 1) = DWParam(27 + SPInd - 1) + SPdm1
                                DWParam(27 + SPInd) = DWParam(27 + SPInd) + SPd0
                                DWParam(27 + SPInd + 1) = DWParam(27 + SPInd + 1) + SPdp1
                                DWParam(27 + SPInd + 2) = DWParam(27 + SPInd + 2) + SPdp2
                            endif
                        endif
                    end if
                end if
            end if

            if (i < j) then
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    if (dijsq < CutSq(2)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_BW>
                        SPx = DIJ * Const(205)
                        SPind = max(min(int(SPx), 24), 0)
                        SPt = SPx - real(SPind)
                        THISU = Const(206 + SPInd) + SPt * (Const(231 + SPInd) + SPt * (Const(256 + SPInd) + SPt * &
                          & Const(281 + SPInd)))
                        THISW = (Const(231 + SPInd) + SPt * (2.*Const(256 + SPInd) + SPt * 3.*Const(281 + SPInd))) * SPx
                        PEnergy = PEnergy + ThisU
                        Virial = Virial + ThisW
                        ThisWsum = ThisWsum + ThisW
                        if (CalcTerms) then
                            Terms(2) = Terms(2) + ThisU
                        endif
                        if (CalcDUParam) then
                            SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                            SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                            SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                            SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                            if (SPind == 0) then
                                !note: y(-1) = 2 y(0) - y(1) for linearity condition
                                DUParam(52 + SPInd) = DUParam(52 + SPInd) + SPd0 + 2.0d0 * SPdm1
                                DUParam(52 + SPInd + 1) = DUParam(52 + SPInd + 1) + SPdp1 - SPdm1
                                DUParam(52 + SPInd + 2) = DUParam(52 + SPInd + 2) + SPdp2
                            elseif (SPind == 24) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DUParam(52 + SPInd - 1) = DUParam(52 + SPInd - 1) + SPdm1
                                DUParam(52 + SPInd) = DUParam(52 + SPInd) + SPd0 - 0.5d0 * SPdp1 + SPdp2
                            elseif (SPind == 23) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DUParam(52 + SPInd - 1) = DUParam(52 + SPInd - 1) + SPdm1
                                DUParam(52 + SPInd) = DUParam(52 + SPInd) + SPd0
                                DUParam(52 + SPInd + 1) = DUParam(52 + SPInd + 1) + SPdp1 - 0.5d0 * SPdp2
                            else
                                DUParam(52 + SPInd - 1) = DUParam(52 + SPInd - 1) + SPdm1
                                DUParam(52 + SPInd) = DUParam(52 + SPInd) + SPd0
                                DUParam(52 + SPInd + 1) = DUParam(52 + SPInd + 1) + SPdp1
                                DUParam(52 + SPInd + 2) = DUParam(52 + SPInd + 2) + SPdp2
                            endif
                        endif
                        if (CalcDWParam) then
                            SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
                            SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
                            SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
                            SPdp2 = (SPt * SPt * 0.5d0) * SPx
                            if (SPind == 0) then
                                !note: y(-1) = 2 y(0) - y(1) for linearity condition
                                DWParam(52 + SPInd) = DWParam(52 + SPInd) + SPd0 + 2.0d0 * SPdm1
                                DWParam(52 + SPInd + 1) = DWParam(52 + SPInd + 1) + SPdp1 - SPdm1
                                DWParam(52 + SPInd + 2) = DWParam(52 + SPInd + 2) + SPdp2
                            elseif (SPind == 24) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DWParam(52 + SPInd - 1) = DWParam(52 + SPInd - 1) + SPdm1
                                DWParam(52 + SPInd) = DWParam(52 + SPInd) + SPd0 - 0.5d0 * SPdp1 + SPdp2
                            elseif (SPind == 23) then
                                !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                                DWParam(52 + SPInd - 1) = DWParam(52 + SPInd - 1) + SPdm1
                                DWParam(52 + SPInd) = DWParam(52 + SPInd) + SPd0
                                DWParam(52 + SPInd + 1) = DWParam(52 + SPInd + 1) + SPdp1 - 0.5d0 * SPdp2
                            else
                                DWParam(52 + SPInd - 1) = DWParam(52 + SPInd - 1) + SPdm1
                                DWParam(52 + SPInd) = DWParam(52 + SPInd) + SPd0
                                DWParam(52 + SPInd + 1) = DWParam(52 + SPInd + 1) + SPdp1
                                DWParam(52 + SPInd + 2) = DWParam(52 + SPInd + 2) + SPdp2
                            endif
                        endif
                    end if
                end if
            end if

            !update the forces
            if (CalcForce .and. ThisWsum /= 0.d0) then
                Forcei = rij * ThisWsum / dijsq
                Force(i,:) = Force(i,:) + Forcei
                Force(j,:) = Force(j,:) - Forcei
            endif

            if ((AIDi==0) .and. (AIDj==1)) then
                if (dijsq < CutSq(3)) then
                    !local density loop for potential <LD_BW>
                    if (DIJSQ < Const(407)) then
                        THISLOCALRHO = RhoSc_3
                        THISLOCALRHOF = 0.d0
                    else
                        THISLOCALRHO = RhoSc_3 * (Const(408) + DIJSQ * (Const(409) + DIJSQ * (Const(410) + Const(411) * DIJSQ)))
                        THISLOCALRHOF = RhoSc_3 * (Const(412) + DIJSQ * (Const(413) + DIJSQ * Const(414)))
                    endif
                    LocalRho_3 = LocalRho_3 + ThisLocalRho
                    LocalRhoW_3 = LocalRhoW_3 + (ThisLocalRhoF * dijsq)
                    if (CalcForce) then
                        ThisForce = ThisLocalRhoF * rij
                        LocalRhoF_3(j,:) = -ThisForce
                        LocalRhoF_3(i,:) = LocalRhoF_3(i,:) + ThisForce
                    endif
                end if
            end if

        !end of loop j
        enddo

        if (AIDi==0) then
            !local density potential <LD_BW>
            LocalRho = LocalRho_3
            LocalRhoW = LocalRhoW_3
            SPx = (LOCALRHO - Const(404)) * Const(307)
            SPind = max(min(int(SPx), 23), 0)
            SPt = max(min(SPx - real(SPind), 1.d0), 0.d0)
            THISU = Const(308 + SPInd) + SPt * (Const(332 + SPInd) + SPt * (Const(356 + SPInd) + SPt * Const(380 + SPInd)))
            THISDU = (Const(332 + SPInd) + SPt * (2.*Const(356 + SPInd) + SPt * 3.*Const(380 + SPInd))) * Const(307)
            ThisW = ThisdU * LocalRhoW
            PEnergy = PEnergy + ThisU
            Virial = Virial + ThisW
            if (CalcTerms) then
                Terms(3) = Terms(3) + ThisU
            endif
            if (CalcDUParam) then
                SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                if (SPind == 0) then
                    !note y(-1) = y(1) from zero slope condition
                    DUParam(77 + SPInd) = DUParam(77 + SPInd) + SPd0
                    DUParam(77 + SPInd + 1) = DUParam(77 + SPInd + 1) + SPdp1 + SPdm1
                    DUParam(77 + SPInd + 2) = DUParam(77 + SPInd + 2) + SPdp2
                elseif (SPind == 23) then
                    !note y(N) = y(N-2) from zero slope condition
                    DUParam(77 + SPInd - 1) = DUParam(77 + SPInd - 1) + SPdm1
                    DUParam(77 + SPInd) = DUParam(77 + SPInd) + SPd0 + SPdp2
                    DUParam(77 + SPInd + 1) = DUParam(77 + SPInd + 1) + SPdp1
                else
                    DUParam(77 + SPInd - 1) = DUParam(77 + SPInd - 1) + SPdm1
                    DUParam(77 + SPInd) = DUParam(77 + SPInd) + SPd0
                    DUParam(77 + SPInd + 1) = DUParam(77 + SPInd + 1) + SPdp1
                    DUParam(77 + SPInd + 2) = DUParam(77 + SPInd + 2) + SPdp2
                endif
            endif
            if (CalcDWParam) then
                val1 = LOCALRHOW * Const(307)
                SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * val1
                SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * val1
                SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * val1
                SPdp2 = (SPt * SPt * 0.5d0) * val1
                if (SPind == 0) then
                    !note y(-1) = y(1) from zero slope condition
                    DWParam(77 + SPInd) = DWParam(77 + SPInd) + SPd0
                    DWParam(77 + SPInd + 1) = DWParam(77 + SPInd + 1) + SPdp1 + SPdm1
                    DWParam(77 + SPInd + 2) = DWParam(77 + SPInd + 2) + SPdp2
                elseif (SPind == 23) then
                    !note y(N) = y(N-2) from zero slope condition
                    DWParam(77 + SPInd - 1) = DWParam(77 + SPInd - 1) + SPdm1
                    DWParam(77 + SPInd) = DWParam(77 + SPInd) + SPd0 + SPdp2
                    DWParam(77 + SPInd + 1) = DWParam(77 + SPInd + 1) + SPdp1
                else
                    DWParam(77 + SPInd - 1) = DWParam(77 + SPInd - 1) + SPdm1
                    DWParam(77 + SPInd) = DWParam(77 + SPInd) + SPd0
                    DWParam(77 + SPInd + 1) = DWParam(77 + SPInd + 1) + SPdp1
                    DWParam(77 + SPInd + 2) = DWParam(77 + SPInd + 2) + SPdp2
                endif
            endif
            if (CalcForce) then
                Force = Force + ThisdU * LocalRhoF_3
            endif
        end if

    !end of loop i
    enddo

end subroutine


subroutine calcargstats(ANumList, Pos, BoxL, AID, CutSq, Const, Weight, ArgMin, ArgMax, ArgAvg, ArgStd, ArgSum, ArgSumSq, &
  & ArgCount, NList, Dim, NAtom, NTerm, NConst, NArg)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: Dim
    integer, intent(in) :: NAtom
    integer, intent(in) :: NTerm
    integer, intent(in) :: NConst
    integer, intent(in) :: NArg
    integer, dimension(NList), intent(in) :: ANumList
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Pos
    real(8), dimension(0:Dim-1), intent(in) :: BoxL
    integer, dimension(0:NAtom-1), intent(in) :: AID
    real(8), dimension(0:NTerm-1), intent(in) :: CutSq
    real(8), dimension(0:NConst-1), intent(in) :: Const
    real(8), intent(in) :: Weight
    real(8), dimension(0:NArg-1), intent(inout) :: ArgMin
!f2py intent(in, out) :: ArgMin
    real(8), dimension(0:NArg-1), intent(inout) :: ArgMax
!f2py intent(in, out) :: ArgMax
    real(8), dimension(0:NArg-1), intent(out) :: ArgAvg
    real(8), dimension(0:NArg-1), intent(out) :: ArgStd
    real(8), dimension(0:NArg-1), intent(inout) :: ArgSum
!f2py intent(in, out) :: ArgSum
    real(8), dimension(0:NArg-1), intent(inout) :: ArgSumSq
!f2py intent(in, out) :: ArgSumSq
    real(8), dimension(0:NArg-1), intent(inout) :: ArgCount
!f2py intent(in, out) :: ArgCount
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: AIDi
    integer :: AIDj
    integer :: ListIndi
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8) :: ThisLocalRho
    real(8) :: ThisLocalRhoF
    integer :: ArgType
    real(8) :: ArgVal
    real(8) :: RhoSc_3
    real(8) :: LocalRho_3

    !initialization for potential <SP_BB>

    !initialization for potential <SP_WW>

    !initialization for potential <SP_BW>

    !initialization for potential <LD_BW>
    RhoSc_3 = Const(415)

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    !loop over i
    i = -1
    ListIndi = 0
    do while ((NList == 0 .and. i < NAtom-1) .or. (ListIndi < NList))
        if (NList == 0) then
            i = i + 1
        else
            ListIndi = ListIndi + 1
            i = ANumList(ListIndi)
        endif

        Posi = Pos(i,:)
        AIDi = AID(i)

        if (AIDi==0) then
            LocalRho_3 = 0.d0
        end if

        !loop over j
        do j = 0, NAtom-1
            if (i == j) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (i < j) then
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    if (dijsq < CutSq(0)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_BB>
                        ARGVAL = DIJ
                        ARGTYPE = 0
                        if (Weight > 0.d0) then
                            ArgMin(ArgType) = min(ArgMin(ArgType), ArgVal)
                            ArgMax(ArgType) = max(ArgMax(ArgType), ArgVal)
                            ArgCount(ArgType) = ArgCount(ArgType) + Weight
                            ArgSum(ArgType) = ArgSum(ArgType) + ArgVal * Weight
                            ArgSumSq(ArgType) = ArgSumSq(ArgType) + ArgVal * ArgVal * Weight
                        endif
                    end if
                end if
            end if

            if (i < j) then
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    if (dijsq < CutSq(1)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_WW>
                        ARGVAL = DIJ
                        ARGTYPE = 0
                        if (Weight > 0.d0) then
                            ArgMin(1 + ArgType) = min(ArgMin(1 + ArgType), ArgVal)
                            ArgMax(1 + ArgType) = max(ArgMax(1 + ArgType), ArgVal)
                            ArgCount(1 + ArgType) = ArgCount(1 + ArgType) + Weight
                            ArgSum(1 + ArgType) = ArgSum(1 + ArgType) + ArgVal * Weight
                            ArgSumSq(1 + ArgType) = ArgSumSq(1 + ArgType) + ArgVal * ArgVal * Weight
                        endif
                    end if
                end if
            end if

            if (i < j) then
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    if (dijsq < CutSq(2)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_BW>
                        ARGVAL = DIJ
                        ARGTYPE = 0
                        if (Weight > 0.d0) then
                            ArgMin(2 + ArgType) = min(ArgMin(2 + ArgType), ArgVal)
                            ArgMax(2 + ArgType) = max(ArgMax(2 + ArgType), ArgVal)
                            ArgCount(2 + ArgType) = ArgCount(2 + ArgType) + Weight
                            ArgSum(2 + ArgType) = ArgSum(2 + ArgType) + ArgVal * Weight
                            ArgSumSq(2 + ArgType) = ArgSumSq(2 + ArgType) + ArgVal * ArgVal * Weight
                        endif
                    end if
                end if
            end if

            if ((AIDi==0) .and. (AIDj==1)) then
                if (dijsq < CutSq(3)) then
                    !local density loop for potential <LD_BW>
                    if (DIJSQ < Const(407)) then
                        THISLOCALRHO = RhoSc_3
                        THISLOCALRHOF = 0.d0
                    else
                        THISLOCALRHO = RhoSc_3 * (Const(408) + DIJSQ * (Const(409) + DIJSQ * (Const(410) + Const(411) * DIJSQ)))
                        THISLOCALRHOF = RhoSc_3 * (Const(412) + DIJSQ * (Const(413) + DIJSQ * Const(414)))
                    endif
                    LocalRho_3 = LocalRho_3 + ThisLocalRho
                end if
            end if

        !end of loop j
        enddo

        !local density potential <LD_BW>
        ARGVAL = LocalRho_3
        ARGTYPE = 0
        if (Weight > 0.d0) then
            ArgMin(3 + ArgType) = min(ArgMin(3 + ArgType), ArgVal)
            ArgMax(3 + ArgType) = max(ArgMax(3 + ArgType), ArgVal)
            ArgCount(3 + ArgType) = ArgCount(3 + ArgType) + Weight
            ArgSum(3 + ArgType) = ArgSum(3 + ArgType) + ArgVal * Weight
            ArgSumSq(3 + ArgType) = ArgSumSq(3 + ArgType) + ArgVal * ArgVal * Weight
        endif

    !end of loop i
    enddo

    ArgAvg = ArgSum / max(ArgCount, 1.d0)
    ArgStd = ArgSumSq / max(ArgCount, 1.d0)
    ArgStd = sqrt(max(ArgStd - ArgAvg*ArgAvg, 0.d0))

end subroutine


subroutine calcarghist(ANumList, Pos, BoxL, AID, CutSq, Const, ArgHistMin, ArgHistiBinw, Weight, ArgHist, NList, Dim, NAtom, &
  & NTerm, NConst, NArg, NHist)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: Dim
    integer, intent(in) :: NAtom
    integer, intent(in) :: NTerm
    integer, intent(in) :: NConst
    integer, intent(in) :: NArg
    integer, intent(in) :: NHist
    integer, dimension(NList), intent(in) :: ANumList
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Pos
    real(8), dimension(0:Dim-1), intent(in) :: BoxL
    integer, dimension(0:NAtom-1), intent(in) :: AID
    real(8), dimension(0:NTerm-1), intent(in) :: CutSq
    real(8), dimension(0:NConst-1), intent(in) :: Const
    real(8), dimension(0:NArg-1), intent(in) :: ArgHistMin
    real(8), dimension(0:NArg-1), intent(in) :: ArgHistiBinw
    real(8), intent(in) :: Weight
    real(8), dimension(0:NHist-1), intent(inout) :: ArgHist
!f2py intent(in, out) :: ArgHist
    real(8), dimension(0:Dim-1) :: iBoxL
    logical :: DoMinImage
    integer :: i
    integer :: j
    integer :: AIDi
    integer :: AIDj
    integer :: ListIndi
    real(8), dimension(0:Dim-1) :: Posi
    real(8), dimension(0:Dim-1) :: Posj
    real(8), dimension(0:Dim-1) :: rij
    real(8) :: dijsq
    real(8) :: dij
    real(8) :: ThisLocalRho
    real(8) :: ThisLocalRhoF
    integer :: ArgType
    real(8) :: ArgVal
    integer :: BinInd
    real(8) :: RhoSc_3
    real(8) :: LocalRho_3

    !initialization for potential <SP_BB>

    !initialization for potential <SP_WW>

    !initialization for potential <SP_BW>

    !initialization for potential <LD_BW>
    RhoSc_3 = Const(415)

    DoMinImage = any(BoxL > 0.d0)
    iBoxL = 1.d0 / max(1.d-300, BoxL)

    !loop over i
    i = -1
    ListIndi = 0
    do while ((NList == 0 .and. i < NAtom-1) .or. (ListIndi < NList))
        if (NList == 0) then
            i = i + 1
        else
            ListIndi = ListIndi + 1
            i = ANumList(ListIndi)
        endif

        Posi = Pos(i,:)
        AIDi = AID(i)

        if (AIDi==0) then
            LocalRho_3 = 0.d0
        end if

        !loop over j
        do j = 0, NAtom-1
            if (i == j) cycle

            Posj = Pos(j,:)
            AIDj = AID(j)

            rij = Posj - Posi
            if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)
            dijsq = dot_product(rij, rij)
            dij = -1.d0

            if (i < j) then
                if (((AIDi==0) .and. (AIDj==0)) .or. ((AIDj==0) .and. (AIDi==0))) then
                    if (dijsq < CutSq(0)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_BB>
                        ARGVAL = DIJ
                        ARGTYPE = 0
                        if (10000 > 0) then
                            BinInd = int((ArgVal - ArgHistMin(ArgType)) * ArgHistiBinw(ArgType))
                            if (BinInd >= 0 .and. BinInd < 10000) ArgHist(10000*ArgType + BinInd) = &
                              & ArgHist(10000*ArgType + BinInd) + Weight
                        endif
                    end if
                end if
            end if

            if (i < j) then
                if (((AIDi==1) .and. (AIDj==1)) .or. ((AIDj==1) .and. (AIDi==1))) then
                    if (dijsq < CutSq(1)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_WW>
                        ARGVAL = DIJ
                        ARGTYPE = 0
                        if (10000 > 0) then
                            BinInd = int((ArgVal - ArgHistMin(1 + ArgType)) * ArgHistiBinw(1 + ArgType))
                            if (BinInd >= 0 .and. BinInd < 10000) ArgHist(10000 + 10000*ArgType + BinInd) = &
                              & ArgHist(10000 + 10000*ArgType + BinInd) + Weight
                        endif
                    end if
                end if
            end if

            if (i < j) then
                if (((AIDi==1) .and. (AIDj==0)) .or. ((AIDj==1) .and. (AIDi==0))) then
                    if (dijsq < CutSq(2)) then
                        if (dij < 0.d0) dij = sqrt(dijsq)
                        !pair interactions for potential <SP_BW>
                        ARGVAL = DIJ
                        ARGTYPE = 0
                        if (10000 > 0) then
                            BinInd = int((ArgVal - ArgHistMin(2 + ArgType)) * ArgHistiBinw(2 + ArgType))
                            if (BinInd >= 0 .and. BinInd < 10000) ArgHist(20000 + 10000*ArgType + BinInd) = &
                              & ArgHist(20000 + 10000*ArgType + BinInd) + Weight
                        endif
                    end if
                end if
            end if

            if ((AIDi==0) .and. (AIDj==1)) then
                if (dijsq < CutSq(3)) then
                    !local density loop for potential <LD_BW>
                    if (DIJSQ < Const(407)) then
                        THISLOCALRHO = RhoSc_3
                        THISLOCALRHOF = 0.d0
                    else
                        THISLOCALRHO = RhoSc_3 * (Const(408) + DIJSQ * (Const(409) + DIJSQ * (Const(410) + Const(411) * DIJSQ)))
                        THISLOCALRHOF = RhoSc_3 * (Const(412) + DIJSQ * (Const(413) + DIJSQ * Const(414)))
                    endif
                    LocalRho_3 = LocalRho_3 + ThisLocalRho
                end if
            end if

        !end of loop j
        enddo

        !local density potential <LD_BW>
        ARGVAL = LocalRho_3
        ARGTYPE = 0
        if (10000 > 0) then
            BinInd = int((ArgVal - ArgHistMin(3 + ArgType)) * ArgHistiBinw(3 + ArgType))
            if (BinInd >= 0 .and. BinInd < 10000) ArgHist(30000 + 10000*ArgType + BinInd) = &
              & ArgHist(30000 + 10000*ArgType + BinInd) + Weight
        endif

    !end of loop i
    enddo

end subroutine


subroutine calcargeval(CutSq, Const, CalcDUParam, CalcDWParam, CalcTerms, ArgHistMin, ArgHistBinw, ArgHist, PEnergy, Virial, &
  & DUParam, DDUParam, DWParam, DDWParam, Terms, NConst, NArg, NHist, NDParam, NDDParam, NTerm)
    implicit none
    integer, intent(in) :: NConst
    integer, intent(in) :: NArg
    integer, intent(in) :: NHist
    integer, intent(in) :: NDParam
    integer, intent(in) :: NDDParam
    integer, intent(in) :: NTerm
    real(8), dimension(0:NTerm-1), intent(in) :: CutSq
    real(8), dimension(0:NConst-1), intent(in) :: Const
    logical, intent(in) :: CalcDUParam
    logical, intent(in) :: CalcDWParam
    logical, intent(in) :: CalcTerms
    real(8), dimension(0:NArg-1), intent(in) :: ArgHistMin
    real(8), dimension(0:NArg-1), intent(in) :: ArgHistBinw
    real(8), dimension(0:NHist-1), intent(in) :: ArgHist
    real(8), intent(out) :: PEnergy
    real(8), intent(out) :: Virial
    real(8), dimension(0:NDParam-1), intent(inout) :: DUParam
!f2py intent(in, out) :: DUParam
    real(8), dimension(0:NDDParam-1), intent(inout) :: DDUParam
!f2py intent(in, out) :: DDUParam
    real(8), dimension(0:NDParam-1), intent(inout) :: DWParam
!f2py intent(in, out) :: DWParam
    real(8), dimension(0:NDDParam-1), intent(inout) :: DDWParam
!f2py intent(in, out) :: DDWParam
    real(8), dimension(0:NTerm-1), intent(inout) :: Terms
!f2py intent(in, out) :: Terms
    integer :: i
    real(8) :: dijsq
    real(8) :: dij
    real(8) :: ThisU
    real(8) :: ThisW
    real(8) :: LocalRho
    real(8) :: LocalRhoW
    integer :: ArgType
    real(8) :: ArgVal
    real(8) :: ThisHist
    real(8) :: SPt
    integer :: SPInd
    real(8) :: SPx
    real(8) :: SPdm1
    real(8) :: SPd0
    real(8) :: SPdp1
    real(8) :: SPdp2
    real(8) :: val1

    !compute initial quantities
    PEnergy = 0.d0
    Virial = 0.d0
    if (CalcDUParam) then
        DUParam = 0.d0
        DDUParam = 0.d0
    endif
    if (CalcDWParam) then
        DWParam = 0.d0
        DDWParam = 0.d0
    endif
    if (CalcTerms) Terms = 0.d0

    !initialization for potential <SP_BB>

    !initialization for potential <SP_WW>

    !initialization for potential <SP_BW>

    !initialization for potential <LD_BW>

    !potential <SP_BB>
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = ArgHist(10000*ArgType + i)
            ArgVal = ArgHistMin(ArgType) + ArgHistBinw(ArgType) * (0.5d0 + i)
            if (ArgVal*ArgVal > CutSq(0)) cycle
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            SPx = DIJ * Const(1)
            SPind = max(min(int(SPx), 24), 0)
            SPt = SPx - real(SPind)
            THISU = Const(2 + SPInd) + SPt * (Const(27 + SPInd) + SPt * (Const(52 + SPInd) + SPt * Const(77 + SPInd)))
            THISW = (Const(27 + SPInd) + SPt * (2.*Const(52 + SPInd) + SPt * 3.*Const(77 + SPInd))) * SPx
            PEnergy = PEnergy + ThisU * ThisHist
            if (CalcTerms) then
                Terms(0) = Terms(0) + ThisU * ThisHist
            endif
            if (CalcDUParam) then
                SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                if (SPind == 0) then
                    !note: y(-1) = 2 y(0) - y(1) for linearity condition
                    DUParam(2 + SPInd) = DUParam(2 + SPInd) + (ThisHist) * (SPd0 + 2.0d0 * SPdm1)
                    DUParam(2 + SPInd + 1) = DUParam(2 + SPInd + 1) + (ThisHist) * (SPdp1 - SPdm1)
                    DUParam(2 + SPInd + 2) = DUParam(2 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 24) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DUParam(2 + SPInd - 1) = DUParam(2 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(2 + SPInd) = DUParam(2 + SPInd) + (ThisHist) * (SPd0 - 0.5d0 * SPdp1 + SPdp2)
                elseif (SPind == 23) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DUParam(2 + SPInd - 1) = DUParam(2 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(2 + SPInd) = DUParam(2 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(2 + SPInd + 1) = DUParam(2 + SPInd + 1) + (ThisHist) * (SPdp1 - 0.5d0 * SPdp2)
                else
                    DUParam(2 + SPInd - 1) = DUParam(2 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(2 + SPInd) = DUParam(2 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(2 + SPInd + 1) = DUParam(2 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DUParam(2 + SPInd + 2) = DUParam(2 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
            if (CalcDWParam) then
                SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
                SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
                SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
                SPdp2 = (SPt * SPt * 0.5d0) * SPx
                if (SPind == 0) then
                    !note: y(-1) = 2 y(0) - y(1) for linearity condition
                    DWParam(2 + SPInd) = DWParam(2 + SPInd) + (ThisHist) * (SPd0 + 2.0d0 * SPdm1)
                    DWParam(2 + SPInd + 1) = DWParam(2 + SPInd + 1) + (ThisHist) * (SPdp1 - SPdm1)
                    DWParam(2 + SPInd + 2) = DWParam(2 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 24) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DWParam(2 + SPInd - 1) = DWParam(2 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(2 + SPInd) = DWParam(2 + SPInd) + (ThisHist) * (SPd0 - 0.5d0 * SPdp1 + SPdp2)
                elseif (SPind == 23) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DWParam(2 + SPInd - 1) = DWParam(2 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(2 + SPInd) = DWParam(2 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(2 + SPInd + 1) = DWParam(2 + SPInd + 1) + (ThisHist) * (SPdp1 - 0.5d0 * SPdp2)
                else
                    DWParam(2 + SPInd - 1) = DWParam(2 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(2 + SPInd) = DWParam(2 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(2 + SPInd + 1) = DWParam(2 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DWParam(2 + SPInd + 2) = DWParam(2 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
        enddo
    enddo

    !potential <SP_WW>
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = ArgHist(10000 + 10000*ArgType + i)
            ArgVal = ArgHistMin(1 + ArgType) + ArgHistBinw(1 + ArgType) * (0.5d0 + i)
            if (ArgVal*ArgVal > CutSq(1)) cycle
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            SPx = DIJ * Const(103)
            SPind = max(min(int(SPx), 24), 0)
            SPt = SPx - real(SPind)
            THISU = Const(104 + SPInd) + SPt * (Const(129 + SPInd) + SPt * (Const(154 + SPInd) + SPt * Const(179 + SPInd)))
            THISW = (Const(129 + SPInd) + SPt * (2.*Const(154 + SPInd) + SPt * 3.*Const(179 + SPInd))) * SPx
            PEnergy = PEnergy + ThisU * ThisHist
            if (CalcTerms) then
                Terms(1) = Terms(1) + ThisU * ThisHist
            endif
            if (CalcDUParam) then
                SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                if (SPind == 0) then
                    !note: y(-1) = 2 y(0) - y(1) for linearity condition
                    DUParam(27 + SPInd) = DUParam(27 + SPInd) + (ThisHist) * (SPd0 + 2.0d0 * SPdm1)
                    DUParam(27 + SPInd + 1) = DUParam(27 + SPInd + 1) + (ThisHist) * (SPdp1 - SPdm1)
                    DUParam(27 + SPInd + 2) = DUParam(27 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 24) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DUParam(27 + SPInd - 1) = DUParam(27 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(27 + SPInd) = DUParam(27 + SPInd) + (ThisHist) * (SPd0 - 0.5d0 * SPdp1 + SPdp2)
                elseif (SPind == 23) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DUParam(27 + SPInd - 1) = DUParam(27 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(27 + SPInd) = DUParam(27 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(27 + SPInd + 1) = DUParam(27 + SPInd + 1) + (ThisHist) * (SPdp1 - 0.5d0 * SPdp2)
                else
                    DUParam(27 + SPInd - 1) = DUParam(27 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(27 + SPInd) = DUParam(27 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(27 + SPInd + 1) = DUParam(27 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DUParam(27 + SPInd + 2) = DUParam(27 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
            if (CalcDWParam) then
                SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
                SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
                SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
                SPdp2 = (SPt * SPt * 0.5d0) * SPx
                if (SPind == 0) then
                    !note: y(-1) = 2 y(0) - y(1) for linearity condition
                    DWParam(27 + SPInd) = DWParam(27 + SPInd) + (ThisHist) * (SPd0 + 2.0d0 * SPdm1)
                    DWParam(27 + SPInd + 1) = DWParam(27 + SPInd + 1) + (ThisHist) * (SPdp1 - SPdm1)
                    DWParam(27 + SPInd + 2) = DWParam(27 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 24) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DWParam(27 + SPInd - 1) = DWParam(27 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(27 + SPInd) = DWParam(27 + SPInd) + (ThisHist) * (SPd0 - 0.5d0 * SPdp1 + SPdp2)
                elseif (SPind == 23) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DWParam(27 + SPInd - 1) = DWParam(27 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(27 + SPInd) = DWParam(27 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(27 + SPInd + 1) = DWParam(27 + SPInd + 1) + (ThisHist) * (SPdp1 - 0.5d0 * SPdp2)
                else
                    DWParam(27 + SPInd - 1) = DWParam(27 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(27 + SPInd) = DWParam(27 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(27 + SPInd + 1) = DWParam(27 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DWParam(27 + SPInd + 2) = DWParam(27 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
        enddo
    enddo

    !potential <SP_BW>
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = ArgHist(20000 + 10000*ArgType + i)
            ArgVal = ArgHistMin(2 + ArgType) + ArgHistBinw(2 + ArgType) * (0.5d0 + i)
            if (ArgVal*ArgVal > CutSq(2)) cycle
            DIJ = ARGVAL
            DIJSQ = DIJ * DIJ
            SPx = DIJ * Const(205)
            SPind = max(min(int(SPx), 24), 0)
            SPt = SPx - real(SPind)
            THISU = Const(206 + SPInd) + SPt * (Const(231 + SPInd) + SPt * (Const(256 + SPInd) + SPt * Const(281 + SPInd)))
            THISW = (Const(231 + SPInd) + SPt * (2.*Const(256 + SPInd) + SPt * 3.*Const(281 + SPInd))) * SPx
            PEnergy = PEnergy + ThisU * ThisHist
            if (CalcTerms) then
                Terms(2) = Terms(2) + ThisU * ThisHist
            endif
            if (CalcDUParam) then
                SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                if (SPind == 0) then
                    !note: y(-1) = 2 y(0) - y(1) for linearity condition
                    DUParam(52 + SPInd) = DUParam(52 + SPInd) + (ThisHist) * (SPd0 + 2.0d0 * SPdm1)
                    DUParam(52 + SPInd + 1) = DUParam(52 + SPInd + 1) + (ThisHist) * (SPdp1 - SPdm1)
                    DUParam(52 + SPInd + 2) = DUParam(52 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 24) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DUParam(52 + SPInd - 1) = DUParam(52 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(52 + SPInd) = DUParam(52 + SPInd) + (ThisHist) * (SPd0 - 0.5d0 * SPdp1 + SPdp2)
                elseif (SPind == 23) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DUParam(52 + SPInd - 1) = DUParam(52 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(52 + SPInd) = DUParam(52 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(52 + SPInd + 1) = DUParam(52 + SPInd + 1) + (ThisHist) * (SPdp1 - 0.5d0 * SPdp2)
                else
                    DUParam(52 + SPInd - 1) = DUParam(52 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(52 + SPInd) = DUParam(52 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(52 + SPInd + 1) = DUParam(52 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DUParam(52 + SPInd + 2) = DUParam(52 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
            if (CalcDWParam) then
                SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
                SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
                SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
                SPdp2 = (SPt * SPt * 0.5d0) * SPx
                if (SPind == 0) then
                    !note: y(-1) = 2 y(0) - y(1) for linearity condition
                    DWParam(52 + SPInd) = DWParam(52 + SPInd) + (ThisHist) * (SPd0 + 2.0d0 * SPdm1)
                    DWParam(52 + SPInd + 1) = DWParam(52 + SPInd + 1) + (ThisHist) * (SPdp1 - SPdm1)
                    DWParam(52 + SPInd + 2) = DWParam(52 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 24) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DWParam(52 + SPInd - 1) = DWParam(52 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(52 + SPInd) = DWParam(52 + SPInd) + (ThisHist) * (SPd0 - 0.5d0 * SPdp1 + SPdp2)
                elseif (SPind == 23) then
                    !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
                    DWParam(52 + SPInd - 1) = DWParam(52 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(52 + SPInd) = DWParam(52 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(52 + SPInd + 1) = DWParam(52 + SPInd + 1) + (ThisHist) * (SPdp1 - 0.5d0 * SPdp2)
                else
                    DWParam(52 + SPInd - 1) = DWParam(52 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(52 + SPInd) = DWParam(52 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(52 + SPInd + 1) = DWParam(52 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DWParam(52 + SPInd + 2) = DWParam(52 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
        enddo
    enddo

    !potential <LD_BW>
    do ArgType = 0, 1 - 1
        do i = 0, 10000 - 1
            ThisHist = ArgHist(30000 + 10000*ArgType + i)
            ArgVal = ArgHistMin(3 + ArgType) + ArgHistBinw(3 + ArgType) * (0.5d0 + i)
            LOCALRHO = ARGVAL
            SPx = (LOCALRHO - Const(404)) * Const(307)
            SPind = max(min(int(SPx), 23), 0)
            SPt = max(min(SPx - real(SPind), 1.d0), 0.d0)
            THISU = Const(308 + SPInd) + SPt * (Const(332 + SPInd) + SPt * (Const(356 + SPInd) + SPt * Const(380 + SPInd)))
            PEnergy = PEnergy + ThisU * ThisHist
            if (CalcTerms) then
                Terms(3) = Terms(3) + ThisU * ThisHist
            endif
            if (CalcDUParam) then
                SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
                SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
                SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
                SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
                if (SPind == 0) then
                    !note y(-1) = y(1) from zero slope condition
                    DUParam(77 + SPInd) = DUParam(77 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(77 + SPInd + 1) = DUParam(77 + SPInd + 1) + (ThisHist) * (SPdp1 + SPdm1)
                    DUParam(77 + SPInd + 2) = DUParam(77 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 23) then
                    !note y(N) = y(N-2) from zero slope condition
                    DUParam(77 + SPInd - 1) = DUParam(77 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(77 + SPInd) = DUParam(77 + SPInd) + (ThisHist) * (SPd0 + SPdp2)
                    DUParam(77 + SPInd + 1) = DUParam(77 + SPInd + 1) + (ThisHist) * (SPdp1)
                else
                    DUParam(77 + SPInd - 1) = DUParam(77 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DUParam(77 + SPInd) = DUParam(77 + SPInd) + (ThisHist) * (SPd0)
                    DUParam(77 + SPInd + 1) = DUParam(77 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DUParam(77 + SPInd + 2) = DUParam(77 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
            if (CalcDWParam) then
                val1 = LOCALRHOW * Const(307)
                SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * val1
                SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * val1
                SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * val1
                SPdp2 = (SPt * SPt * 0.5d0) * val1
                if (SPind == 0) then
                    !note y(-1) = y(1) from zero slope condition
                    DWParam(77 + SPInd) = DWParam(77 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(77 + SPInd + 1) = DWParam(77 + SPInd + 1) + (ThisHist) * (SPdp1 + SPdm1)
                    DWParam(77 + SPInd + 2) = DWParam(77 + SPInd + 2) + (ThisHist) * (SPdp2)
                elseif (SPind == 23) then
                    !note y(N) = y(N-2) from zero slope condition
                    DWParam(77 + SPInd - 1) = DWParam(77 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(77 + SPInd) = DWParam(77 + SPInd) + (ThisHist) * (SPd0 + SPdp2)
                    DWParam(77 + SPInd + 1) = DWParam(77 + SPInd + 1) + (ThisHist) * (SPdp1)
                else
                    DWParam(77 + SPInd - 1) = DWParam(77 + SPInd - 1) + (ThisHist) * (SPdm1)
                    DWParam(77 + SPInd) = DWParam(77 + SPInd) + (ThisHist) * (SPd0)
                    DWParam(77 + SPInd + 1) = DWParam(77 + SPInd + 1) + (ThisHist) * (SPdp1)
                    DWParam(77 + SPInd + 2) = DWParam(77 + SPInd + 2) + (ThisHist) * (SPdp2)
                endif
            endif
        enddo
    enddo

end subroutine


subroutine calcmeasures(Pos, BoxL, SID, MeasureAll, NMeasure, StepNum, CycleNum, Weight, KEnergy, PEnergy, TEnergy, NDOF, &
  & ConservesMomentum, kB, TempSet, Virial, NMolID, DUParam, DDUParam, DWParam, DDWParam, MasterF, MasterI, SDCPosOrigin, &
  & SDCSums, SDCCounts, NAtom, NMID, NDParam, NDDParam, NMasterF, NMasterI, SDCNOrigin, SDCNAtom, Dim, SDCNSID, SDCNBin)
    implicit none
    integer, intent(in) :: NAtom
    integer, intent(in) :: NMID
    integer, intent(in) :: NDParam
    integer, intent(in) :: NDDParam
    integer, intent(in) :: NMasterF
    integer, intent(in) :: NMasterI
    integer, intent(in) :: SDCNOrigin
    integer, intent(in) :: SDCNAtom
    integer, intent(in) :: Dim
    integer, intent(in) :: SDCNSID
    integer, intent(in) :: SDCNBin
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Pos
    real(8), dimension(0:Dim-1), intent(in) :: BoxL
    integer, dimension(0:NAtom-1), intent(in) :: SID
    logical, intent(in) :: MeasureAll
    integer, intent(in) :: NMeasure
    integer, intent(in) :: StepNum
    integer, intent(in) :: CycleNum
    real(8), intent(in) :: Weight
    real(8), intent(in) :: KEnergy
    real(8), intent(in) :: PEnergy
    real(8), intent(in) :: TEnergy
    real(8), intent(in) :: NDOF
    logical, intent(in) :: ConservesMomentum
    real(8), intent(in) :: kB
    real(8), intent(in) :: TempSet
    real(8), intent(in) :: Virial
    integer, dimension(0:NMID-1), intent(in) :: NMolID
    real(8), dimension(0:NDParam-1), intent(in) :: DUParam
    real(8), dimension(0:NDDParam-1), intent(in) :: DDUParam
    real(8), dimension(0:NDParam-1), intent(in) :: DWParam
    real(8), dimension(0:NDDParam-1), intent(in) :: DDWParam
    real(8), dimension(0:NMasterF-1), intent(inout) :: MasterF
!f2py intent(in, out) :: MasterF
    integer, dimension(0:NMasterI-1), intent(inout) :: MasterI
!f2py intent(in, out) :: MasterI
    real(8), dimension(0:SDCNOrigin-1,0:SDCNAtom-1,0:Dim-1), intent(inout) :: SDCPosOrigin
!f2py intent(in, out) :: SDCPosOrigin
    real(8), dimension(0:SDCNSID-1,0:SDCNBin-1), intent(inout) :: SDCSums
!f2py intent(in, out) :: SDCSums
    real(8), dimension(0:SDCNSID-1,0:SDCNBin-1), intent(inout) :: SDCCounts
!f2py intent(in, out) :: SDCCounts
    integer :: i
    integer :: j
    integer :: SIDi
    logical, dimension(0:NMeasure-1) :: UseMeasure
    integer :: v1
    integer :: v2
    real(8) :: Val0
    integer :: SDCOr
    integer :: SDCInd
    real(8) :: SDCrsq

    !###### determine which measures to use ######
    !measure KEnergy
    if (MasterI(2)==0) then
        UseMeasure(0) = .false.
    else
        if (MeasureAll) then
            UseMeasure(0) = .true.
            MasterF(0) = 0.d0
        elseif (MasterI(0) > 0 .and. mod(StepNum, MasterI(0))==0) then
            UseMeasure(0) = .true.
            MasterF(0) = 0.d0
        elseif (MasterI(1) > 0 .and. mod(CycleNum, MasterI(1))==0) then
            UseMeasure(0) = .true.
            MasterF(0) = 0.d0
        else
            UseMeasure(0) = .false.
        endif
    endif

    !measure PEnergy
    if (MasterI(5)==0) then
        UseMeasure(1) = .false.
    else
        if (MeasureAll) then
            UseMeasure(1) = .true.
            MasterF(4) = 0.d0
        elseif (MasterI(3) > 0 .and. mod(StepNum, MasterI(3))==0) then
            UseMeasure(1) = .true.
            MasterF(4) = 0.d0
        elseif (MasterI(4) > 0 .and. mod(CycleNum, MasterI(4))==0) then
            UseMeasure(1) = .true.
            MasterF(4) = 0.d0
        else
            UseMeasure(1) = .false.
        endif
    endif

    !measure TEnergy
    if (MasterI(8)==0) then
        UseMeasure(2) = .false.
    else
        if (MeasureAll) then
            UseMeasure(2) = .true.
            MasterF(8) = 0.d0
        elseif (MasterI(6) > 0 .and. mod(StepNum, MasterI(6))==0) then
            UseMeasure(2) = .true.
            MasterF(8) = 0.d0
        elseif (MasterI(7) > 0 .and. mod(CycleNum, MasterI(7))==0) then
            UseMeasure(2) = .true.
            MasterF(8) = 0.d0
        else
            UseMeasure(2) = .false.
        endif
    endif

    !measure KTemp
    if (MasterI(11)==0) then
        UseMeasure(3) = .false.
    else
        if (MeasureAll) then
            UseMeasure(3) = .true.
            MasterF(12) = 0.d0
        elseif (MasterI(9) > 0 .and. mod(StepNum, MasterI(9))==0) then
            UseMeasure(3) = .true.
            MasterF(12) = 0.d0
        elseif (MasterI(10) > 0 .and. mod(CycleNum, MasterI(10))==0) then
            UseMeasure(3) = .true.
            MasterF(12) = 0.d0
        else
            UseMeasure(3) = .false.
        endif
    endif

    !measure Pressure
    if (MasterI(14)==0) then
        UseMeasure(4) = .false.
    else
        if (MeasureAll) then
            UseMeasure(4) = .true.
            MasterF(16) = 0.d0
        elseif (MasterI(12) > 0 .and. mod(StepNum, MasterI(12))==0) then
            UseMeasure(4) = .true.
            MasterF(16) = 0.d0
        elseif (MasterI(13) > 0 .and. mod(CycleNum, MasterI(13))==0) then
            UseMeasure(4) = .true.
            MasterF(16) = 0.d0
        else
            UseMeasure(4) = .false.
        endif
    endif

    !measure Virial
    if (MasterI(17)==0) then
        UseMeasure(5) = .false.
    else
        if (MeasureAll) then
            UseMeasure(5) = .true.
            MasterF(20) = 0.d0
        elseif (MasterI(15) > 0 .and. mod(StepNum, MasterI(15))==0) then
            UseMeasure(5) = .true.
            MasterF(20) = 0.d0
        elseif (MasterI(16) > 0 .and. mod(CycleNum, MasterI(16))==0) then
            UseMeasure(5) = .true.
            MasterF(20) = 0.d0
        else
            UseMeasure(5) = .false.
        endif
    endif

    !measure Rho:B
    if (MasterI(20)==0) then
        UseMeasure(6) = .false.
    else
        if (MeasureAll) then
            UseMeasure(6) = .true.
            MasterF(24) = 0.d0
        elseif (MasterI(18) > 0 .and. mod(StepNum, MasterI(18))==0) then
            UseMeasure(6) = .true.
            MasterF(24) = 0.d0
        elseif (MasterI(19) > 0 .and. mod(CycleNum, MasterI(19))==0) then
            UseMeasure(6) = .true.
            MasterF(24) = 0.d0
        else
            UseMeasure(6) = .false.
        endif
    endif

    !measure Rho:W
    if (MasterI(23)==0) then
        UseMeasure(7) = .false.
    else
        if (MeasureAll) then
            UseMeasure(7) = .true.
            MasterF(28) = 0.d0
        elseif (MasterI(21) > 0 .and. mod(StepNum, MasterI(21))==0) then
            UseMeasure(7) = .true.
            MasterF(28) = 0.d0
        elseif (MasterI(22) > 0 .and. mod(CycleNum, MasterI(22))==0) then
            UseMeasure(7) = .true.
            MasterF(28) = 0.d0
        else
            UseMeasure(7) = .false.
        endif
    endif

    !measure DUParam
    if (MasterI(26)==0) then
        UseMeasure(8) = .false.
    else
        if (MeasureAll) then
            UseMeasure(8) = .true.
            MasterF(32:133) = 0.d0
        elseif (MasterI(24) > 0 .and. mod(StepNum, MasterI(24))==0) then
            UseMeasure(8) = .true.
            MasterF(32:133) = 0.d0
        elseif (MasterI(25) > 0 .and. mod(CycleNum, MasterI(25))==0) then
            UseMeasure(8) = .true.
            MasterF(32:133) = 0.d0
        else
            UseMeasure(8) = .false.
        endif
    endif

    !measure DDUParam
    if (MasterI(131)==0) then
        UseMeasure(9) = .false.
    else
        if (MeasureAll) then
            UseMeasure(9) = .true.
            MasterF(10641:13344) = 0.d0
        elseif (MasterI(129) > 0 .and. mod(StepNum, MasterI(129))==0) then
            UseMeasure(9) = .true.
            MasterF(10641:13344) = 0.d0
        elseif (MasterI(130) > 0 .and. mod(CycleNum, MasterI(130))==0) then
            UseMeasure(9) = .true.
            MasterF(10641:13344) = 0.d0
        else
            UseMeasure(9) = .false.
        endif
    endif

    !measure DWParam
    if (MasterI(134)==0) then
        UseMeasure(10) = .false.
    else
        if (MeasureAll) then
            UseMeasure(10) = .true.
            MasterF(18754:18855) = 0.d0
        elseif (MasterI(132) > 0 .and. mod(StepNum, MasterI(132))==0) then
            UseMeasure(10) = .true.
            MasterF(18754:18855) = 0.d0
        elseif (MasterI(133) > 0 .and. mod(CycleNum, MasterI(133))==0) then
            UseMeasure(10) = .true.
            MasterF(18754:18855) = 0.d0
        else
            UseMeasure(10) = .false.
        endif
    endif

    !measure DDWParam
    if (MasterI(239)==0) then
        UseMeasure(11) = .false.
    else
        if (MeasureAll) then
            UseMeasure(11) = .true.
            MasterF(19061:21764) = 0.d0
        elseif (MasterI(237) > 0 .and. mod(StepNum, MasterI(237))==0) then
            UseMeasure(11) = .true.
            MasterF(19061:21764) = 0.d0
        elseif (MasterI(238) > 0 .and. mod(CycleNum, MasterI(238))==0) then
            UseMeasure(11) = .true.
            MasterF(19061:21764) = 0.d0
        else
            UseMeasure(11) = .false.
        endif
    endif

    !measure DUParamDWParam
    if (MasterI(242)==0) then
        UseMeasure(12) = .false.
    else
        if (MeasureAll) then
            UseMeasure(12) = .true.
            MasterF(27174:37577) = 0.d0
        elseif (MasterI(240) > 0 .and. mod(StepNum, MasterI(240))==0) then
            UseMeasure(12) = .true.
            MasterF(27174:37577) = 0.d0
        elseif (MasterI(241) > 0 .and. mod(CycleNum, MasterI(241))==0) then
            UseMeasure(12) = .true.
            MasterF(27174:37577) = 0.d0
        else
            UseMeasure(12) = .false.
        endif
    endif

    !measure SDC
    if (MasterI(245)==0) then
        UseMeasure(13) = .false.
    else
        if (MeasureAll) then
            UseMeasure(13) = .true.
        elseif (MasterI(243) > 0 .and. mod(StepNum, MasterI(243))==0) then
            UseMeasure(13) = .true.
        elseif (MasterI(244) > 0 .and. mod(CycleNum, MasterI(244))==0) then
            UseMeasure(13) = .true.
        else
            UseMeasure(13) = .false.
        endif
    endif

    !###### Non-atom measures ######
    if (UseMeasure(0)) then
        !measure KEnergy
        MasterF(0) = KEnergy
        Val0 = MasterF(0)
        MasterF(1) = MasterF(1) + Val0*Weight
        MasterF(2) = MasterF(2) + Val0*Val0*Weight
        MasterF(3) = MasterF(3) + Weight
    end if

    if (UseMeasure(1)) then
        !measure PEnergy
        MasterF(4) = PEnergy
        Val0 = MasterF(4)
        MasterF(5) = MasterF(5) + Val0*Weight
        MasterF(6) = MasterF(6) + Val0*Val0*Weight
        MasterF(7) = MasterF(7) + Weight
    end if

    if (UseMeasure(2)) then
        !measure TEnergy
        MasterF(8) = TEnergy
        Val0 = MasterF(8)
        MasterF(9) = MasterF(9) + Val0*Weight
        MasterF(10) = MasterF(10) + Val0*Val0*Weight
        MasterF(11) = MasterF(11) + Weight
    end if

    if (UseMeasure(3)) then
        !measure KTemp
        MasterF(12) = 2.d0 * KEnergy / (kB * merge(NDOF-Dim, NDOF, ConservesMomentum))
        Val0 = MasterF(12)
        MasterF(13) = MasterF(13) + Val0*Weight
        MasterF(14) = MasterF(14) + Val0*Val0*Weight
        MasterF(15) = MasterF(15) + Weight
    end if

    if (UseMeasure(4)) then
        !measure Pressure
        MasterF(16) = (kB*NDOF*TempSet - Virial) / (Dim * product(BoxL))
        Val0 = MasterF(16)
        MasterF(17) = MasterF(17) + Val0*Weight
        MasterF(18) = MasterF(18) + Val0*Val0*Weight
        MasterF(19) = MasterF(19) + Weight
    end if

    if (UseMeasure(5)) then
        !measure Virial
        MasterF(20) = Virial
        Val0 = MasterF(20)
        MasterF(21) = MasterF(21) + Val0*Weight
        MasterF(22) = MasterF(22) + Val0*Val0*Weight
        MasterF(23) = MasterF(23) + Weight
    end if

    if (UseMeasure(6)) then
        !measure Rho:B
        MasterF(24) = float(NMolID(0)) / product(BoxL)
        Val0 = MasterF(24)
        MasterF(25) = MasterF(25) + Val0*Weight
        MasterF(26) = MasterF(26) + Val0*Val0*Weight
        MasterF(27) = MasterF(27) + Weight
    end if

    if (UseMeasure(7)) then
        !measure Rho:W
        MasterF(28) = float(NMolID(1)) / product(BoxL)
        Val0 = MasterF(28)
        MasterF(29) = MasterF(29) + Val0*Weight
        MasterF(30) = MasterF(30) + Val0*Val0*Weight
        MasterF(31) = MasterF(31) + Weight
    end if

    if (UseMeasure(8)) then
        !measure DUParam
        MasterF(32:133) = DUParam
        MasterI(27:128) = MasterI(27:128) + merge(0, 1, DUParam==0.d0)
        MasterF(134:235) = MasterF(134:235) + (MasterF(32:133))*Weight
        do v1 = 0, 102 - 1
            do v2 = 0, 102 - 1
                MasterF(236 + 102*v1 + v2) = MasterF(236 + 102*v1 + v2) + MasterF(32 + v1)*MasterF(32 + v2)*Weight
            enddo
        enddo
        MasterF(10640) = MasterF(10640) + Weight
    end if

    if (UseMeasure(9)) then
        !measure DDUParam
        MasterF(10641:13344) = DDUParam
        MasterF(13345:16048) = MasterF(13345:16048) + (MasterF(10641:13344))*Weight
        MasterF(16049:18752) = MasterF(16049:18752) + (MasterF(10641:13344))*(MasterF(10641:13344))*Weight
        MasterF(18753) = MasterF(18753) + Weight
    end if

    if (UseMeasure(10)) then
        !measure DWParam
        MasterF(18754:18855) = DWParam
        MasterI(135:236) = MasterI(135:236) + merge(0, 1, DWParam==0.d0)
        MasterF(18856:18957) = MasterF(18856:18957) + (MasterF(18754:18855))*Weight
        MasterF(18958:19059) = MasterF(18958:19059) + (MasterF(18754:18855))*(MasterF(18754:18855))*Weight
        MasterF(19060) = MasterF(19060) + Weight
    end if

    if (UseMeasure(11)) then
        !measure DDWParam
        MasterF(19061:21764) = DDWParam
        MasterF(21765:24468) = MasterF(21765:24468) + (MasterF(19061:21764))*Weight
        MasterF(24469:27172) = MasterF(24469:27172) + (MasterF(19061:21764))*(MasterF(19061:21764))*Weight
        MasterF(27173) = MasterF(27173) + Weight
    end if

    if (UseMeasure(12)) then
        !measure DUParamDWParam
        do i = 0, NDParam - 1
            do j = 0, NDParam - 1
                MasterF(27174 + i*NDParam + j) = DUParam(i) * DWParam(j)
            enddo
        enddo
        MasterF(37578:47981) = MasterF(37578:47981) + (MasterF(27174:37577))*Weight
        MasterF(47982:58385) = MasterF(47982:58385) + (MasterF(27174:37577))*(MasterF(27174:37577))*Weight
        MasterF(58386) = MasterF(58386) + Weight
    end if

    if (UseMeasure(13)) then
        !measure SDC
        if (NAtom /= SDCNAtom) then
            print *, "Error: cannot use SDC calculations with variable atom number."
            stop
        endif
        do SDCOr = 0, MasterI(248) - 1
            if (StepNum == MasterI(249) * SDCOr + MasterI(243)) then
                SDCPosOrigin(SDCOr, :, :) = Pos
            endif
            SDCInd = int((StepNum - MasterI(246) - SDCOr * MasterI(249)) / MasterI(243) - 1)
            if (SDCInd < 0 .or. SDCInd >= SDCNBin) cycle
            do i = 0, NAtom - 1
                SIDi = SID(i)
                SDCrsq = sum((Pos(i,:) - SDCPosOrigin(SDCOr, i, :))**2)
                SDCSums(SIDi, SDCInd) = SDCSums(SIDi, SDCInd) + SDCrsq
                SDCCounts(SIDi, SDCInd) = SDCCounts(SIDi, SDCInd) + 1.d0
            enddo
        enddo
    end if

end subroutine


subroutine vvquench1(ANumList, dt, Force, iMass, PEnergy, Pos, Vel, KEnergy, TEnergy, NList, NAtom, Dim)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: NAtom
    integer, intent(in) :: Dim
    integer, dimension(NList), intent(in) :: ANumList
    real(8), intent(in) :: dt
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Force
    real(8), dimension(0:NAtom-1), intent(in) :: iMass
    real(8), intent(in) :: PEnergy
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Pos
!f2py intent(in, out) :: Pos
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Vel
!f2py intent(in, out) :: Vel
    real(8), intent(inout) :: KEnergy
!f2py intent(in, out) :: KEnergy
    real(8), intent(out) :: TEnergy
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: ListInd

    !velocity verlet quench 1
    dtsq2 = dt*dt*0.5
    dt2 = 0.5*dt
    idt = 1./dt
    KEnergy = 0.d0
    TEnergy = KEnergy + PEnergy

    if (NList == 0) then
        do i = 0, NAtom-1
            Accel = Force(i,:) * iMass(i)
            Vel(i,:) = 0.d0
            Pos(i,:) = Pos(i,:) + dtsq2*Accel
        enddo
    else
        do ListInd = 1, NList
            i = ANumList(ListInd)
            Accel = Force(i,:) * iMass(i)
            Vel(i,:) = 0.d0
            Pos(i,:) = Pos(i,:) + dtsq2*Accel
        enddo
    endif

end subroutine


subroutine vvupdatekenergy(ANumList, Mass, PEnergy, Vel, KEnergy, TEnergy, NList, NAtom, Dim)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: NAtom
    integer, intent(in) :: Dim
    integer, dimension(NList), intent(in) :: ANumList
    real(8), dimension(0:NAtom-1), intent(in) :: Mass
    real(8), intent(in) :: PEnergy
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Vel
!f2py intent(in, out) :: Vel
    real(8), intent(inout) :: KEnergy
!f2py intent(in, out) :: KEnergy
    real(8), intent(out) :: TEnergy
    integer :: i
    integer :: ListInd

    !update kinetic energy
    KEnergy = 0.d0
    if (NList == 0) then
        do i = 0, NAtom-1
            KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
        enddo
    else
        do ListInd = 1, NList
            i = ANumList(ListInd)
            KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
        enddo
    endif
    KEnergy = KEnergy * 0.5d0
    TEnergy = KEnergy + PEnergy

end subroutine


subroutine vvintegrate1(ANumList, dt, Thermostat, Force, Mass, iMass, sqrtMass, PEnergy, NDOF, TempSet, kB, &
  & AndersenCollisionFreq, AndersenStepFreq, RemoveCOMStepFreq, Pos, Vel, KEnergy, TEnergy, NH_Glogs, NH_Vlogs, NH_Xlogs, &
  & NH_QMass, AndersenStep, RemoveCOMStep, LangevinGamma, NList, NAtom, Dim, NH_N)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: NAtom
    integer, intent(in) :: Dim
    integer, intent(in) :: NH_N
    integer, dimension(NList), intent(in) :: ANumList
    real(8), intent(in) :: dt
    integer, intent(in) :: Thermostat
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Force
    real(8), dimension(0:NAtom-1), intent(in) :: Mass
    real(8), dimension(0:NAtom-1), intent(in) :: iMass
    real(8), dimension(0:NAtom-1), intent(in) :: sqrtMass
    real(8), intent(in) :: PEnergy
    integer, intent(in) :: NDOF
    real(8), intent(in) :: TempSet
    real(8), intent(in) :: kB
    real(8), intent(in) :: AndersenCollisionFreq
    integer, intent(in) :: AndersenStepFreq
    integer, intent(in) :: RemoveCOMStepFreq
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Pos
!f2py intent(in, out) :: Pos
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Vel
!f2py intent(in, out) :: Vel
    real(8), intent(inout) :: KEnergy
!f2py intent(in, out) :: KEnergy
    real(8), intent(out) :: TEnergy
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_Glogs
!f2py intent(in, out) :: NH_Glogs
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_Vlogs
!f2py intent(in, out) :: NH_Vlogs
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_Xlogs
!f2py intent(in, out) :: NH_Xlogs
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_QMass
!f2py intent(in, out) :: NH_QMass
    integer, intent(inout) :: AndersenStep
!f2py intent(in, out) :: AndersenStep
    integer, intent(inout) :: RemoveCOMStep
!f2py intent(in, out) :: RemoveCOMStep
    real(8), intent(inout) :: LangevinGamma
!f2py intent(in, out) :: LangevinGamma
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: ListInd
    real(8) :: NH_wdti0
    real(8) :: NH_wdti1
    real(8) :: NH_wdti2
    real(8) :: NH_wdti3
    real(8) :: NH_kT
    real(8) :: NH_NkT
    real(8) :: NH_scale
    real(8) :: AA
    real(8) :: NH_akin
    integer :: inos
    real(8) :: rn
    real(8), dimension(0:Dim-1) :: ranvec
    real(8) :: dtfreq
    real(8) :: sqrtkt
    real(8) :: langevin1
    real(8) :: langevin2
    real(8), parameter :: randomvelclip = -1.d0
    external :: ran2
    external :: ran2normarray

    if (Thermostat == 1) then
        !do an andersen massive collision update
        if (mod(AndersenStep, AndersenStepFreq)==0) then
            AndersenStep = 0
            sqrtkt = sqrt(kB * TempSet)
            if (NList == 0) then
                do i = 0, NAtom-1
                    call ran2normarray(Dim, ranvec)
                    !clip the limits on the random variate for stability
                    if (randomvelclip > 0.d0) then
                        ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                    endif
                    ranvec = ranvec * sqrtkt / sqrtMass(i)
                    Vel(i,:) = ranvec
                enddo
            else
                do ListInd = 1, NList
                    i = ANumList(ListInd)
                    call ran2normarray(Dim, ranvec)
                    !clip the limits on the random variate for stability
                    if (randomvelclip > 0.d0) then
                        ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                    endif
                    ranvec = ranvec * sqrtkt / sqrtMass(i)
                    Vel(i,:) = ranvec
                enddo
            endif
        endif
        AndersenStep = AndersenStep + 1
    endif

    if (Thermostat == 2) then
        !do an andersen particle collision update
        dtfreq = dt * AndersenCollisionFreq
        sqrtkt = sqrt(kB * TempSet)
        if (NList == 0) then
            do i = 0, NAtom-1
                call ran2(rn)
                if (rn < dtfreq) then
                    call ran2normarray(Dim, ranvec)
                    !clip the limits on the random variate for stability
                    if (randomvelclip > 0.d0) then
                        ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                    endif
                    ranvec = ranvec * sqrtkt / sqrtMass(i)
                    Vel(i,:) = ranvec
                endif
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                call ran2(rn)
                if (rn < dtfreq) then
                    call ran2normarray(Dim, ranvec)
                    !clip the limits on the random variate for stability
                    if (randomvelclip > 0.d0) then
                        ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                    endif
                    ranvec = ranvec * sqrtkt / sqrtMass(i)
                    Vel(i,:) = ranvec
                endif
            enddo
        endif
    endif

    !remove the center of mass
    if (RemoveCOMStepFreq > 0) then
        if (mod(RemoveCOMStep, RemoveCOMStepFreq)==0) then
            RemoveCOMStep = 0
            if (NList == 0) then
                ranvec = 0.d0
                do i = 0, NAtom-1
                    ranvec = ranvec + Mass(i) * Vel(i,:)
                enddo
                ranvec = ranvec / NAtom
                do i = 0, NAtom-1
                    Vel(i,:) = Vel(i,:) - ranvec * iMass(i)
                enddo
            else
                ranvec = 0.d0
                do ListInd = 1, NList
                    i = ANumList(ListInd)
                    ranvec = ranvec + Mass(i) * Vel(i,:)
                enddo
                ranvec = ranvec / NList
                do ListInd = 1, NList
                    i = ANumList(ListInd)
                    Vel(i,:) = Vel(i,:) - ranvec * iMass(i)
                enddo
            endif
        endif
        RemoveCOMStep = RemoveCOMStep + 1
    endif

    !NOSE-HOOVER ROUTINES
    if (Thermostat == 3) then
        !update kinetic energy
        KEnergy = 0.d0
        if (NList == 0) then
            do i = 0, NAtom-1
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
        endif
        KEnergy = KEnergy * 0.5d0
        TEnergy = KEnergy + PEnergy
        !set frequently used variables
        NH_wdti0 = dt
        NH_wdti1 = dt * 0.5d0
        NH_wdti2 = dt * 0.25d0
        NH_wdti3 = dt * 0.125
        NH_kT = TempSet * kB
        NH_NkT = NH_kT * dble(NDOF - Dim)
        NH_scale = 1.D0
        !get kinetic energy
        NH_akin = 2.d0 * KEnergy
        !update the forces
        NH_Glogs(0) = (NH_akin - NH_NkT) / NH_Qmass(0)
        !update the thermostat velocities
        NH_Vlogs(NH_N-1) = NH_Vlogs(NH_N-1) + NH_Glogs(NH_N-1) * NH_wdti2
        do inos = 1, NH_N - 1
            AA = exp( -NH_wdti3 * NH_Vlogs(NH_N-inos) )
            NH_Vlogs(NH_N-inos-1) = NH_Vlogs(NH_N-inos-1)*AA*AA + NH_wdti2*NH_Glogs(NH_N-inos-1)*AA
        enddo
        !update the particle velocities
        AA = exp( -NH_wdti1 * NH_Vlogs(0) )
        NH_scale = NH_scale * AA
        !update the forces
        NH_Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / NH_Qmass(0)
        !update the thermostat positions
        do inos = 0, NH_N - 1
            NH_Xlogs(inos) = NH_Xlogs(inos) + NH_Vlogs(inos) * NH_wdti1
        enddo
        !update the thermostat velocities
        do inos = 1, NH_N-1
            AA = exp( -NH_wdti3 * NH_Vlogs(inos) )
            NH_Vlogs(inos-1) = NH_Vlogs(inos-1)*AA*AA + NH_wdti2*NH_Glogs(inos-1)*AA
            NH_Glogs(inos) = (NH_Qmass(inos-1)*NH_Vlogs(inos-1)*NH_Vlogs(inos-1)-NH_kT) / NH_Qmass(inos)
        enddo
        NH_Vlogs(NH_N-1) = NH_Vlogs(NH_N-1) + NH_wdti2*NH_Glogs(NH_N-1)
        !update the particle velocities
        if (NH_scale > 0.) then
            Vel = Vel * NH_scale
        endif
    endif

    !velocity verlet integration 1
    dtsq2 = dt*dt*0.5
    dt2 = 0.5*dt
    idt = 1./dt

    if (Thermostat == 4) then
        !langevin thermostat
        !based on Bussi and Parrinello, Physical Review E 75, 056707, 2007
        langevin1 = exp(-0.5d0 * LangevinGamma * dt)
        langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
        if (NList == 0) then
            do i = 0, NAtom-1
                Accel = Force(i,:) * iMass(i)
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
                Pos(i,:) = Pos(i,:) + dt*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                Accel = Force(i,:) * iMass(i)
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
                Pos(i,:) = Pos(i,:) + dt*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif
    else
        !normal constant energy dynamics
        if (NList == 0) then
            do i = 0, NAtom-1
                Accel = Force(i,:) * iMass(i)
                Pos(i,:) = Pos(i,:) + dt*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                Accel = Force(i,:) * iMass(i)
                Pos(i,:) = Pos(i,:) + dt*Vel(i,:) + dtsq2*Accel
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif
    endif

end subroutine


subroutine vvintegrate2(ANumList, dt, Thermostat, Force, Mass, iMass, sqrtMass, PEnergy, NDOF, TempSet, kB, Vel, KEnergy, &
  & TEnergy, NH_Glogs, NH_Vlogs, NH_Xlogs, NH_QMass, LangevinGamma, TEnergySum, TEnergySqSum, NList, NAtom, Dim, NH_N)
    implicit none
    integer, intent(in) :: NList
    integer, intent(in) :: NAtom
    integer, intent(in) :: Dim
    integer, intent(in) :: NH_N
    integer, dimension(NList), intent(in) :: ANumList
    real(8), intent(in) :: dt
    integer, intent(in) :: Thermostat
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(in) :: Force
    real(8), dimension(0:NAtom-1), intent(in) :: Mass
    real(8), dimension(0:NAtom-1), intent(in) :: iMass
    real(8), dimension(0:NAtom-1), intent(in) :: sqrtMass
    real(8), intent(in) :: PEnergy
    integer, intent(in) :: NDOF
    real(8), intent(in) :: TempSet
    real(8), intent(in) :: kB
    real(8), dimension(0:NAtom-1,0:Dim-1), intent(inout) :: Vel
!f2py intent(in, out) :: Vel
    real(8), intent(inout) :: KEnergy
!f2py intent(in, out) :: KEnergy
    real(8), intent(out) :: TEnergy
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_Glogs
!f2py intent(in, out) :: NH_Glogs
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_Vlogs
!f2py intent(in, out) :: NH_Vlogs
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_Xlogs
!f2py intent(in, out) :: NH_Xlogs
    real(8), dimension(0:NH_N-1), intent(inout) :: NH_QMass
!f2py intent(in, out) :: NH_QMass
    real(8), intent(inout) :: LangevinGamma
!f2py intent(in, out) :: LangevinGamma
    real(8), intent(inout) :: TEnergySum
!f2py intent(in, out) :: TEnergySum
    real(8), intent(inout) :: TEnergySqSum
!f2py intent(in, out) :: TEnergySqSum
    real(8) :: dtsq2
    real(8) :: dt2
    real(8) :: idt
    real(8), dimension(0:Dim-1) :: Accel
    integer :: i
    integer :: ListInd
    real(8) :: NH_wdti0
    real(8) :: NH_wdti1
    real(8) :: NH_wdti2
    real(8) :: NH_wdti3
    real(8) :: NH_kT
    real(8) :: NH_NkT
    real(8) :: NH_scale
    real(8) :: AA
    real(8) :: NH_akin
    integer :: inos
    real(8), dimension(0:Dim-1) :: ranvec
    real(8) :: langevin1
    real(8) :: langevin2
    real(8), parameter :: randomvelclip = -1.d0
    external :: ran2normarray

    !velocity verlet integration 2
    dtsq2 = dt*dt*0.5
    dt2 = 0.5*dt
    idt = 1./dt

    if (Thermostat == 4) then
        !langevin thermostat
        langevin1 = exp(-0.5d0 * LangevinGamma * dt)
        langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
        if (NList == 0) then
            do i = 0, NAtom-1
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
                Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
            enddo
        endif
    else
        !normal constant energy dynamics
        if (NList == 0) then
            do i = 0, NAtom-1
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                Accel = Force(i,:) * iMass(i)
                Vel(i,:) = Vel(i,:) + dt2*Accel
            enddo
        endif
    endif

    !NOSE-HOOVER ROUTINES
    if (Thermostat == 3) then
        !update kinetic energy
        KEnergy = 0.d0
        if (NList == 0) then
            do i = 0, NAtom-1
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
        else
            do ListInd = 1, NList
                i = ANumList(ListInd)
                KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
            enddo
        endif
        KEnergy = KEnergy * 0.5d0
        TEnergy = KEnergy + PEnergy
        !set frequently used variables
        NH_wdti0 = dt
        NH_wdti1 = dt * 0.5d0
        NH_wdti2 = dt * 0.25d0
        NH_wdti3 = dt * 0.125
        NH_kT = TempSet * kB
        NH_NkT = NH_kT * dble(NDOF - Dim)
        NH_scale = 1.D0
        !get kinetic energy
        NH_akin = 2.d0 * KEnergy
        !update the forces
        NH_Glogs(0) = (NH_akin - NH_NkT) / NH_Qmass(0)
        !update the thermostat velocities
        NH_Vlogs(NH_N-1) = NH_Vlogs(NH_N-1) + NH_Glogs(NH_N-1) * NH_wdti2
        do inos = 1, NH_N - 1
            AA = exp( -NH_wdti3 * NH_Vlogs(NH_N-inos) )
            NH_Vlogs(NH_N-inos-1) = NH_Vlogs(NH_N-inos-1)*AA*AA + NH_wdti2*NH_Glogs(NH_N-inos-1)*AA
        enddo
        !update the particle velocities
        AA = exp( -NH_wdti1 * NH_Vlogs(0) )
        NH_scale = NH_scale * AA
        !update the forces
        NH_Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / NH_Qmass(0)
        !update the thermostat positions
        do inos = 0, NH_N - 1
            NH_Xlogs(inos) = NH_Xlogs(inos) + NH_Vlogs(inos) * NH_wdti1
        enddo
        !update the thermostat velocities
        do inos = 1, NH_N-1
            AA = exp( -NH_wdti3 * NH_Vlogs(inos) )
            NH_Vlogs(inos-1) = NH_Vlogs(inos-1)*AA*AA + NH_wdti2*NH_Glogs(inos-1)*AA
            NH_Glogs(inos) = (NH_Qmass(inos-1)*NH_Vlogs(inos-1)*NH_Vlogs(inos-1)-NH_kT) / NH_Qmass(inos)
        enddo
        NH_Vlogs(NH_N-1) = NH_Vlogs(NH_N-1) + NH_wdti2*NH_Glogs(NH_N-1)
        !update the particle velocities
        if (NH_scale > 0.) then
            Vel = Vel * NH_scale
        endif
    endif

    !update kinetic energy
    KEnergy = 0.d0
    if (NList == 0) then
        do i = 0, NAtom-1
            KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
        enddo
    else
        do ListInd = 1, NList
            i = ANumList(ListInd)
            KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
        enddo
    endif
    KEnergy = KEnergy * 0.5d0
    TEnergy = KEnergy + PEnergy

    !update running sums for total energy
    TEnergySum = TEnergySum + TEnergy
    TEnergySqSum = TEnergySqSum + TEnergy*TEnergy

end subroutine





subroutine ran2(r)
    implicit none
    real(8), intent(out) :: r
    integer, parameter :: NTAB = 32
    integer, parameter :: IM1=2147483563, IM2=2147483399, IMM1=2147483562
    integer, parameter :: IA1=40014, IA2=40692, IQ1=53668, IQ2=52774
    integer, parameter :: IR1=12211, IR2=3791, NDIV=1+IMM1/NTAB  
    real(8), parameter :: AM=1.d0/IM1, EPS=1.2d-7, RNMX=1.d0-EPS
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    integer :: j, k

    if (state(1).le.0) then
        state(1)=max(-state(1),1)
        state(2)=state(1)
        do j=NTAB+8,1,-1
             k=state(1)/IQ1
             state(1)=IA1*(state(1)-k*IQ1)-k*IR1
             if (state(1).lt.0) state(1)=state(1)+IM1
             if (j.le.NTAB) state(3+j)=state(1)
        enddo
        state(3)=state(4)
    endif
    k=state(1)/IQ1
    state(1)=IA1*(state(1)-k*IQ1)-k*IR1
    if (state(1).lt.0) state(1)=state(1)+IM1
    k=state(2)/IQ2
    state(2)=IA2*(state(2)-k*IQ2)-k*IR2
    if (state(2).lt.0) state(2)=state(2)+IM2
    j=1+state(3)/NDIV
    state(3)=state(3+j)-state(2)
    state(3+j)=state(1)
    if(state(3).lt.1)state(3)=state(3)+IMM1
    r=min(AM*state(3),RNMX)
end subroutine

subroutine ran2seed(seedval)
    implicit none
    integer, intent(in) :: seedval
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    state(1) = abs(seedval)
end subroutine

subroutine ran2array(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    integer :: i
    external :: ran2
    real(8) :: r
    do i = 1, n
        call ran2(r)
        rarray(i) = r
    enddo
end subroutine

subroutine ran2norm(r)
    implicit none
    real(8), intent(out) :: r
    external :: ran2
    real(8) :: r1, r2, rsq
    real(8), save :: rsaved = 0.1
    logical, save :: hassaved = .false.
    if (hassaved) then
        r = rsaved
        hassaved = .false.
    else
        rsq = 2.d0
        do while (rsq == 0.d0 .or. rsq >= 1.d0)
            call ran2(r1)
            call ran2(r2)
            r1 = 2.d0 * r1 - 1.d0
            r2 = 2.d0 * r2 - 1.d0
            rsq = r1*r1 + r2*r2
        enddo
        rsq = sqrt(-2.d0 * log(rsq) / rsq)
        r = r1 * rsq
        rsaved = r2 * rsq
        hassaved = .true.
    endif
end subroutine

subroutine ran2normarray(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    external :: ran2norm
    integer :: i
    real(8) :: r
    do i = 1, n
        call ran2norm(r)
        rarray(i) = r
    enddo
end subroutine


subroutine erfc(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the complementary error function erfc(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
          & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
          & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    if (x < 0.) value = 2. - value
end subroutine erfc

subroutine erf(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the error function erf(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
          & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
          & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    if (x < 0.) value = 2. - value
    value = 1. - value
end subroutine erf


integer function GetPairIndFromij(i, j)
    integer, intent(in) :: i, j
    if (i > j) then
        GetPairIndFromij = i * (i + 1) / 2 + j
    else
        GetPairIndFromij = j * (j + 1) / 2 + i
    endif      
end function

subroutine GetijFromPairInd(ind, i, j)
    integer, intent(in) :: ind
    integer, intent(out) :: i, j
    real(8) :: r
    r = dble(ind) + 0.001d0
    r = sqrt(1.d0 + 8.d0 * r) * 0.5d0 - 0.5d0
    i = int(r)
    j = ind - i*(i+1)/2
end subroutine



function modulehash()
    character(len=40) :: modulehash
    modulehash = '6b8261efc105f04485320e72350dccb34290ad5f'
end function
