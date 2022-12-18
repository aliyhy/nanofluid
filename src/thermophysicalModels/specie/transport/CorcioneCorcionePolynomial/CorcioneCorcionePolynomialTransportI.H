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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::CorcioneCorcionePolynomialTransport
(
//    
    const Thermo& t,
    const Polynomial<PolySize>& muCoeffs,
    const Polynomial<PolySize>& kappaCoeffs,
    const Polynomial<PolySize>& CpCoeffs,
    const Polynomial<PolySize>& rhoCoeffs,
    const scalar dP,
    const scalar molWeight,
    const scalar rho0,
    const scalar Tfr,
    const Polynomial<PolySize>& kappaPCoeffs     
)
:
  //  
    Thermo(t),
    muCoeffs_(muCoeffs),
    kappaCoeffs_(kappaCoeffs),
    CpCoeffs_(CpCoeffs),
    rhoCoeffs_(rhoCoeffs),
    dP_(dP),
    molWeight_(molWeight),
    rho0_(rho0),
    Tfr_(Tfr),
    kappaPCoeffs_(kappaPCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::CorcioneCorcionePolynomialTransport
(
//    
    const word& name,
    const CorcioneCorcionePolynomialTransport& pt
)
:
    Thermo(name, pt),
//    
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_),
    CpCoeffs_(pt.CpCoeffs_),
    rhoCoeffs_(pt.rhoCoeffs_),
    dP_(pt.dP_),
    molWeight_(pt.molWeight_),
    rho0_(pt.rho0_),
    Tfr_(pt.Tfr_),
    kappaPCoeffs_(pt.kappaPCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>>
Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<CorcioneCorcionePolynomialTransport<Thermo, PolySize>>::New(*this);
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>>
Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::New(/**/ const dictionary& dict)
{
    return autoPtr<CorcioneCorcionePolynomialTransport<Thermo, PolySize>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::mu
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
  scalar df_ =  0.1*pow(6.0*molWeight_/(6.0221408*pow(10.,23.)*3.1416*rho0_),1./3.);
  return muBF(mesh, p, T)/(1.0-34.87*pow(dP_/df_,-0.3)*pow(this->phiNP(mesh),1.03));
}

template<class Thermo, int PolySize>
inline Foam::scalar Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::kappa
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    scalar Re_ =
     2.* rhoCoeffs_.value(T)*1.380649
     /pow(10.,23.)*T/(3.1416*pow(muBF(mesh, p, T),2.)*dP_);
    scalar Pr_ =  this->CpBF(mesh, p, T)*muBF(mesh, p, T)/kappaBF(mesh, p, T);

    return 
      kappaBF(mesh, p, T)*(1. + 4.4* pow (Re_,0.4)*pow(Pr_,0.66)
      *pow(T/Tfr_,10.)*pow(kappaPCoeffs_.value(T)
      /kappaBF(mesh, p, T),0.03)*pow(this->phiNP(mesh),0.66));
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::alphah
(
    const fvMesh& mesh, const scalar p,const scalar T
) const
{
    return kappa(mesh, p, T)/this->Cp(mesh, p, T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::muBF
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    return muCoeffs_.value(T);
}

template<class Thermo, int PolySize>
inline Foam::scalar Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::kappaBF
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    return kappaCoeffs_.value(T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::alphahBF
(
    const fvMesh& mesh, const scalar p,const scalar T
) const
{
    return kappaBF(mesh,p, T)/this->CpBF(mesh, p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::operator+=
(
    const CorcioneCorcionePolynomialTransport<Thermo, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = pt.Y()/this->Y();

        muCoeffs_ = Y1*muCoeffs_ + Y2*pt.muCoeffs_;
        kappaCoeffs_ = Y1*kappaCoeffs_ + Y2*pt.kappaCoeffs_;
    }
}


template<class Thermo, int PolySize>
inline void Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize> Foam::operator+
(
    const CorcioneCorcionePolynomialTransport<Thermo, PolySize>& pt1,
    const CorcioneCorcionePolynomialTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return CorcioneCorcionePolynomialTransport<Thermo>
        (
            t,
            0,
            pt1.muCoeffs_,
            pt1.kappaCoeffs_
        );
    }
    else
    {
        scalar Y1 = pt1.Y()/t.Y();
        scalar Y2 = pt2.Y()/t.Y();

        return CorcioneCorcionePolynomialTransport<Thermo, PolySize>
        (
            t,
            Y1*pt1.muCoeffs_ + Y2*pt2.muCoeffs_,
            Y1*pt1.kappaCoeffs_ + Y2*pt2.kappaCoeffs_
        );
    }
}


template<class Thermo, int PolySize>
inline Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const CorcioneCorcionePolynomialTransport<Thermo, PolySize>& pt
)
{
    return CorcioneCorcionePolynomialTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.muCoeffs_,
        pt.kappaCoeffs_,
        pt.CpCoeffs_,
        pt.rhoCoeffs_,
        pt.dP_,
        pt.molWeight_,
        pt.rho0_,
        pt.Tfr_,
        pt.kappaPCoeffs_
    );
}


// ************************************************************************* //
