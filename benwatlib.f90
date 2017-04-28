SUBROUTINE SKBI(Pos, BoxL, NAtom, AtomTypes, AtomType_a, AtomType_b, L0, KBI, RandIter)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, AtomType_a, AtomType_b, RandIter
    INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
    REAL(8), INTENT(IN) :: L0
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(OUT), DIMENSION(0:2) :: KBI
    
    INTEGER :: iter, i
    REAL(8) :: BoxVol, SubVol
    REAL(8), DIMENSION(0:2) :: rij, invBoxL, Origin, rsq, Lsq
    REAL(8), DIMENSION(0:RandIter-1) :: Na, Nb
    
    KBI = 0.d0
    Lsq = (/L0*L0, L0*L0, L0*L0/)
    Nb = 0.d0
    Na = 0.d0
    
    Origin = 0.d0
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    BoxVol = BoxL(0) * BoxL(1) * BoxL(2)
    SubVol = L0 * L0 * L0

    ! calculate particle fluctuations in random (cubical) subvolumes 
    DO iter = 0, RandIter-1
        ! find origin
        Origin = (/ 0.5*L0 + (BoxL(0)-L0) * RAND(0), &
                    0.5*L0 + (BoxL(1)-L0) * RAND(0), &
                    0.5*L0 + (BoxL(2)-L0) * RAND(0)  /)
        
        DO i = 0, NAtom-1
            rij = Pos(i,:) - Origin(:)
            rij = rij - BoxL * DNINT(invBoxL * rij) ! reimage
            IF ( (ABS(rij(0)) <= L0) .AND. (ABS(rij(1)) <= L0) .AND. (ABS(rij(2)) <= L0) ) THEN
                IF (AtomTypes(i) == AtomType_a) Na(iter) = Na(iter) + 1
                IF (AtomTypes(i) == AtomType_b) Nb(iter) = Nb(iter) + 1
            END IF
        END DO
    END DO
    
    ! compute KBI from definition
    KBI(0) = SubVol * ( SUM(Na*Na) / (SUM(Na)*SUM(Na)) - 1 - RandIter/SUM(Na) )
    KBI(1) = SubVol * ( SUM(Nb*Nb) / (SUM(Nb)*SUM(Nb)) - 1 - RandIter/SUM(Nb) )
    KBI(2) = SubVol * ( SUM(Na*Nb) / (SUM(Na)*SUM(Nb)) - 1 )
    
END SUBROUTINE
            

            
    
    
