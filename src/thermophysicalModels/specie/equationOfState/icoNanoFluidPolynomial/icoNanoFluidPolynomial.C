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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
icoNanoFluidPolynomial<Specie, PolySize>::icoNanoFluidPolynomial(const dictionary& dict)
:
    Specie(dict),
    rhoCoeffs_(dict.subDict("equationOfState").lookup(coeffsName("rho"))),
    rhoP_(dict.subDict("nanoParticle").get<scalar>("rhoP"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
void icoNanoFluidPolynomial<Specie, PolySize>::write(Ostream& os) const
{
    Specie::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("equationOfState");
        os.writeEntry(coeffsName("rho"), rhoCoeffs_);
        os.writeEntry("rhoP", rhoP_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Ostream& operator<<(Ostream& os, const icoNanoFluidPolynomial<Specie, PolySize>& ip)
{
    ip.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
