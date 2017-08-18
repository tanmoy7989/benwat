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


SUBROUTINE ZDENSITY(Pos, BoxL, NAtom, rho, NBins)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, NBins
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(OUT), DIMENSION(0:NBins-1) :: rho
    
    INTEGER :: i, iter
    REAL(8) :: di, z0, dz
    REAL(8), DIMENSION(0:2) :: invBoxL
    
    dz = (BoxL(2) - 0.d0) / NBins
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    DO iter = 0, NBins-1
        z0 = 0.d0 + dz * (iter + 0.5)
        DO i = 0, NAtom-1
            di = Pos(i,2) - z0
            di = di - BoxL(2) * DNINT(invBoxL(2) * di) ! reimage along Z axis
            IF ( ABS(di) <= 0.5 * dz) rho(iter) = rho(iter) + 1
        END DO
    END DO            
END SUBROUTINE


SUBROUTINE BINONGRID(Pos, BoxL, xcenters, ycenters, zcenters, NAtom, NBins, &
idx, idy, idz)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, NBins
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:NBins-1) :: xcenters, ycenters, zcenters
    INTEGER, INTENT(OUT), DIMENSION(0:NAtom-1) :: idx, idy ,idz
    
    INTEGER :: i, bin
    REAL(8) :: dist_x, dist_y, dist_z, x0, y0, z0, dx, dy, dz
    REAL(8), DIMENSION(0:2) :: Posi, invBoxL
    
    dx = xcenters(1) - xcenters(0)
    dy = ycenters(1) - ycenters(0)
    dz = zcenters(1) - zcenters(0)
    idx = 0.d0
    idy = 0.d0
    idz = 0.d0
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    
    DO i = 0, NAtom-1
        Posi = Pos(i,:)
        DO bin = 0, NBins - 1
            x0 = 0.d0 + (bin + 0.5) * dz
            y0 = 0.d0 + (bin + 0.5) * dy
            z0 = 0.d0 + (bin + 0.5) * dz
            dist_x = Posi(0) - x0
            dist_y = Posi(1) - y0
            dist_z = Posi(1) - z0
            dist_x = dist_x - BoxL(0) * DNINT(invBoxL(0) * dist_x)
            dist_y = dist_y - BoxL(1) * DNINT(invBoxL(1) * dist_y)
            dist_z = dist_z - BoxL(2) * DNINT(invBoxL(2) * dist_z)
            
            IF ( ABS(dist_x) <= 0.5 * dx) idx(i) = bin
            IF ( ABS(dist_y) <= 0.5 * dy) idy(i) = bin
            IF ( ABS(dist_z) <= 0.5 * dz) idz(i) = bin
        END DO 
    ENDDO
END SUBROUTINE
        

SUBROUTINE GRIDINSERT(Pos, BoxL, NAtom, InsPos, NIns, &
AtomTypes, AtomType_B, AtomType_W, Rcav, Ncav)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom, NIns, AtomType_B, AtomType_W
    REAL(8), INTENT(IN) :: Rcav
    INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(IN), DIMENSION(0:NIns-1, 0:2) :: InsPos
    INTEGER, INTENT(OUT) :: Ncav
    
    INTEGER :: i, j
    LOGICAL :: isCav, hasB, hasW
    REAL(8) :: rB, rW, rsq, r1sq, r2sq, Rcavsq
    REAL(8), DIMENSION(0:2) :: invBoxL, Posi, Posj, rij
    
    !VDW Radii (Bondii et. al.)
    rB = 1.77 ; rW = 1.52 ! HS radii of O2 (i.e. SPC/E water)
    r1sq = (rB + Rcav) * (rB + Rcav)
    r2sq = (rW + Rcav) * (rW + Rcav)
    Rcavsq = Rcav * Rcav
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    Ncav = 0.d0
    
    DO i = 0, NIns - 1
        isCav = .TRUE.
        Posi = InsPos(i,:)
        DO j = 0, NAtom - 1
            ! obtain cavity location from supplied grid
            Posj = Pos(j,:)
            rij = Posi - Posj
            rij = rij - BoxL * ANINT(rij * invBoxL) ! minimage
            rsq = SUM(rij * rij)
            !hasB = ( (AtomTypes(j) == AtomType_B) .AND. (rsq < r1sq) )
            !hasW = ( (AtomTypes(j) == AtomType_W) .AND. (rsq < r2sq) )
            IF (rsq < Rcavsq) THEN
                isCav = .FALSE.
                EXIT
            END IF
        END DO
        
        IF (isCav) Ncav = Ncav + 1
    END DO
END SUBROUTINE    
    
    
    
    
    
    
    
    
    
   
