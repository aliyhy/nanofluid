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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar B,
    const scalar C
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    B_(B),
    C_(C)
{}


inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.get<scalar>("A")),
    beta_(dict.get<scalar>("beta")),
    Ta_(dict.get<scalar>("Ta")),
    B_(dict.get<scalar>("B")),
    C_(dict.get<scalar>("C"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::LandauTellerReactionRate::operator()
(
    const fvMesh& mesh, const scalar p,
    const scalar T,
    const scalarField&
) const
{
    scalar lta = A_;

    if (mag(beta_) > VSMALL)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0.0;

    if (mag(Ta_) > VSMALL)
    {
        expArg -= Ta_/T;
    }

    if (mag(B_) > VSMALL)
    {
        expArg += B_/cbrt(T);
    }

    if (mag(C_) > VSMALL)
    {
        expArg += C_/pow(T, 2.0/3.0);
    }

    if (mag(expArg) > VSMALL)
    {
        lta *= exp(expArg);
    }

    return lta;
}


inline void Foam::LandauTellerReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("beta", beta_);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("B", B_);
    os.writeEntry("C", C_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LandauTellerReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
