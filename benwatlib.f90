SUBROUTINE SKBI(Pos, BoxL, NAtom, AtomTypes, AtomType_B, AtomType_W, L0, RandIter, G_BB, G_WW, G_BW)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, AtomType_B, AtomType_W, RandIter
    INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
    REAL(8), INTENT(IN) :: L0
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(OUT) :: G_BB, G_WW, G_BW
    
    INTEGER :: i, iter
    REAL(8) :: muB, muW, muBB, muWW, muBW, SVol
    REAL(8), DIMENSION(0:RandIter-1) :: NB, NW 
    REAL(8), DIMENSION(0:2) :: rij, invBoxL, Origin
    
    NB = 0.d0
    NW = 0.d0
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    SVol = L0 * L0 * L0
    
    DO iter = 0, RandIter - 1
        ! find origin of cubical subvolume
        Origin = (/ 0.5*L0 + (BoxL(0)-L0) * RAND(0), &
                    0.5*L0 + (BoxL(1)-L0) * RAND(0), &
                    0.5*L0 + (BoxL(2)-L0) * RAND(0)  /)
        
        ! calculate particle numbers within cube
        DO i = 0, NAtom-1
            rij = Pos(i,:) - Origin(:)
            rij = rij - BoxL * DNINT(invBoxL * rij) ! reimage
            IF ( (ABS(rij(0)) <= 0.5*L0) .AND. (ABS(rij(1)) <= 0.5*L0) .AND. (ABS(rij(2)) <= 0.5*L0) ) THEN
                IF (AtomTypes(i) == AtomType_B) NB(iter) = NB(iter) + 1
                IF (AtomTypes(i) == AtomType_W) NW(iter) = NW(iter) + 1
            END IF
        END DO
    END DO
    
    ! calculate KBIs from particle fluctuations
    muB = SUM(NB) / RandIter
    muW = SUM(NW) / RandIter
    muBB = SUM(NB*NB) / RandIter
    muWW = SUM(NW*NW) / RandIter
    muBW = SUM(NB*NW) / RandIter
    G_BB = SVol * ( (muBB - muB*muB)/(muB*muB) - 1.d0/muB )
    G_WW = SVol * ( (muWW - muW*muW)/(muW*muW) - 1.d0/muW )
    G_BW = SVol * ( (muBW - muB*muW)/(muB*muW) )
END SUBROUTINE
            

SUBROUTINE SLICEDENSITY(Pos, BoxL, NAtom, AtomTypes, AtomType_B, AtomType_W, rhoB, rhoW, NSlice)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, NSlice, AtomType_B, AtomType_W
    INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(OUT), DIMENSION(0:NSlice-1) :: rhoB, rhoW
    
    INTEGER :: i, iter
    REAL(8) :: zi, z0, dz
    REAL(8), DIMENSION(0:2) :: invBoxL
    
    dz = (BoxL(2) - 0.d0) / NSlice
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    
    DO iter = 0, NSlice-1
        z0 = 0.d0 + dz*iter
        DO i = 0, NAtom-1
            zi = Pos(i,2) - z0
            zi = zi - BoxL(2) * DNINT(invBoxL(2) * zi) ! reimage along Z axis
            IF (ABS(zi-z0) <= dz) THEN
                IF (AtomTypes(i) == AtomType_B) rhoB(iter) = rhoB(iter) + 1
                IF (AtomTypes(i) == AtomType_W) rhoW(iter) = rhoW(iter) + 1
            END IF
        END DO
    END DO            
END SUBROUTINE
    
