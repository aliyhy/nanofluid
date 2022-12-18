/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "icoNanoFluidPolynomial.H"
#include "scalarIOField.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoNanoFluidPolynomial<Specie, PolySize>::icoNanoFluidPolynomial
(
    const Specie& sp,
    const Polynomial<PolySize>& rhoCoeffs,
    const scalar rhoP
)
:
    Specie(sp),
    rhoCoeffs_(rhoCoeffs),
    rhoP_(rhoP)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoNanoFluidPolynomial<Specie, PolySize>::icoNanoFluidPolynomial
(
    const word& name,
    const icoNanoFluidPolynomial<Specie, PolySize>& ip
)
:
    Specie(name, ip),
    rhoCoeffs_(ip.rhoCoeffs_),
    rhoP_(ip.rhoP_)
{}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoNanoFluidPolynomial<Specie, PolySize>>
Foam::icoNanoFluidPolynomial<Specie, PolySize>::clone() const
{
    return autoPtr<icoNanoFluidPolynomial<Specie, PolySize>>::New(*this);
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoNanoFluidPolynomial<Specie, PolySize>>
Foam::icoNanoFluidPolynomial<Specie, PolySize>::New(const dictionary& dict)
{
    return autoPtr<icoNanoFluidPolynomial<Specie, PolySize>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::rho
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return rhoP_*phiNP(mesh) + (1.-phiNP(mesh))*rhoBF(mesh, p, T);
}

template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::rhoBF
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return rhoCoeffs_.value(T);
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::H
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::Cp
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::phiNP
(
    const fvMesh& mesh
) const
{
   scalarIOField phiNP_ =
     mesh.objectRegistry::thisDb().lookupObjectRef<scalarIOField>("phiNP");
   return phiNP_[0];
//   return mesh.lookupObject<volScalarField>("phiNP")[mesh.nCells()];
}



template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::E
(
    const fvMesh& mesh, scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::Cv
(
    const fvMesh& mesh, scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::S
(
    const fvMesh& mesh, scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::psi
(
    const fvMesh& mesh, scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::Z
(
    const fvMesh& mesh, scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoNanoFluidPolynomial<Specie, PolySize>::CpMCv
(
    const fvMesh& mesh, scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline void Foam::icoNanoFluidPolynomial<Specie, PolySize>::operator+=
(
    const icoNanoFluidPolynomial<Specie, PolySize>& ip
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ip);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = ip.Y()/this->Y();

        rhoCoeffs_ = Y1*rhoCoeffs_ + Y2*ip.rhoCoeffs_;
    }
}


template<class Specie, int PolySize>
inline void Foam::icoNanoFluidPolynomial<Specie, PolySize>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::icoNanoFluidPolynomial<Specie, PolySize> Foam::operator+
(
    const icoNanoFluidPolynomial<Specie, PolySize>& ip1,
    const icoNanoFluidPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
      + static_cast<const Specie&>(ip2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return icoNanoFluidPolynomial<Specie, PolySize>
        (
            sp,
            ip1.rhoCoeffs_
        );
    }
    else
    {
        const scalar Y1 = ip1.Y()/sp.Y();
        const scalar Y2 = ip2.Y()/sp.Y();

        return icoNanoFluidPolynomial<Specie, PolySize>
        (
            sp,
            Y1*ip1.rhoCoeffs_ + Y2*ip2.rhoCoeffs_
        );
    }
}


template<class Specie, int PolySize>
Foam::icoNanoFluidPolynomial<Specie, PolySize> Foam::operator*
(
    const scalar s,
    const icoNanoFluidPolynomial<Specie, PolySize>& ip
)
{
    return icoNanoFluidPolynomial<Specie, PolySize>
    (
        s*static_cast<const Specie&>(ip),
        ip.rhoCoeffs_,
        ip.rhoP_
    );
}


template<class Specie, int PolySize>
Foam::icoNanoFluidPolynomial<Specie, PolySize> Foam::operator==
(
    const icoNanoFluidPolynomial<Specie, PolySize>& ip1,
    const icoNanoFluidPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
     == static_cast<const Specie&>(ip2)
    );

    const scalar Y1 = ip1.Y()/sp.Y();
    const scalar Y2 = ip2.Y()/sp.Y();

    return icoNanoFluidPolynomial<Specie, PolySize>
    (
        sp,
        Y2*ip2.rhoCoeffs_ - Y1*ip1.rhoCoeffs_,
        Y2*ip2.rhoP_ - Y1*ip1.rhoP_ 
    );
}


// ************************************************************************* //
