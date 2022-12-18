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

#include "hNanoFluidPolynomialThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::hNanoFluidPolynomialThermo<EquationOfState, PolySize>::hNanoFluidPolynomialThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(dict.subDict("thermodynamics").get<scalar>("Hf")),
    Sf_(dict.subDict("thermodynamics").get<scalar>("Sf")),
    CpCoeffs_(dict.subDict("thermodynamics").lookup(coeffsName("Cp"))),
    hCoeffs_(),
    sCoeffs_(),
    rhoCoeffs_(dict.subDict("equationOfState").lookup(coeffsName("rho"))),
    rhoP_(dict.subDict("nanoParticle").get<scalar>("rhoP")),
    CpP_(dict.subDict("nanoParticle").get<scalar>("CpP"))
{
    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(Tstd);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
void Foam::hNanoFluidPolynomialThermo<EquationOfState, PolySize>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("thermodynamics");
        os.writeEntry("Hf", Hf_);
        os.writeEntry("Sf", Sf_);
        os.writeEntry(coeffsName("Cp"), CpCoeffs_);
        os.writeEntry(coeffsName("rho"), rhoCoeffs_);
        os.writeEntry("rhoP", rhoP_);
        os.writeEntry("CpP", CpP_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
