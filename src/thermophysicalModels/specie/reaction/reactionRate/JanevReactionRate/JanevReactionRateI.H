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

inline Foam::JanevReactionRate::JanevReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nb_> b
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    b_(b)
{}


inline Foam::JanevReactionRate::JanevReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.get<scalar>("A")),
    beta_(dict.get<scalar>("beta")),
    Ta_(dict.get<scalar>("Ta")),
    b_(dict.lookup("b"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::JanevReactionRate::operator()
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

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    lta *= exp(expArg);

    return lta;
}


inline void Foam::JanevReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << nl;
    os.writeKeyword("beta") << beta_ << nl;
    os.writeKeyword("Ta") << Ta_ << nl;
    os.writeKeyword("b") << b_ << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const JanevReactionRate& jrr
)
{
    jrr.write(os);
    return os;
}


// ************************************************************************* //
