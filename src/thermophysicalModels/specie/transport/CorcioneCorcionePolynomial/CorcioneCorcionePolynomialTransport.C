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
#include "CorcioneCocciaPolynomialTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::CorcioneCorcionePolynomialTransport
(
//    
    const dictionary& dict
)
:
  //  mesh_(mesh_),
    Thermo(dict),
    muCoeffs_(dict.subDict("transport").lookup(coeffsName("mu"))),
    kappaCoeffs_(dict.subDict("transport").lookup(coeffsName("kappa"))),
    CpCoeffs_(dict.subDict("thermodynamics").lookup(coeffsName("Cp"))),
    rhoCoeffs_(dict.subDict("equationOfState").lookup(coeffsName("rho"))),
    dP_(dict.subDict("nanoParticle").get<scalar>("dP")),
    molWeight_(dict.subDict("specie").get<scalar>("molWeight")),
    rho0_(dict.subDict("transport").get<scalar>("rho0")),
    Tfr_(dict.subDict("transport").get<scalar>("Tfr")),
    kappaPCoeffs_(dict.subDict("nanoParticle").lookup(coeffsName("kappaP")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::CorcioneCorcionePolynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    os.beginBlock(this->name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.writeEntry(coeffsName("mu"), muCoeffs_);
        os.writeEntry(coeffsName("kappa"), kappaCoeffs_);
        os.writeEntry(coeffsName("Cp"), CpCoeffs_);
        os.writeEntry(coeffsName("rho"), rhoCoeffs_);
        os.writeEntry("dP", dP_);
        os.writeEntry("molWeight", molWeight_);
        os.writeEntry("rho0", rho0_);
        os.writeEntry("Tfr", Tfr_);
        os.writeEntry(coeffsName("kappaP"), kappaPCoeffs_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CorcioneCorcionePolynomialTransport<Thermo, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
