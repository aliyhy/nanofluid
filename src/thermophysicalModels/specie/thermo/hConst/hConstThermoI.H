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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const EquationOfState& st,
    const scalar cp,
    const scalar hf
)
:
    EquationOfState(st),
    Cp_(cp),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const word& name,
    const hConstThermo& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConstThermo<EquationOfState>>
Foam::hConstThermo<EquationOfState>::clone() const
{
    return autoPtr<hConstThermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConstThermo<EquationOfState>>
Foam::hConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hConstThermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Cp
(
    const fvMesh& mesh, const scalar p,
    const scalar T
) const
{
    return Cp_ + EquationOfState::Cp(mesh, p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Ha
(
    const fvMesh& mesh, const scalar p, const scalar T
) const
{
    return Hs(mesh, p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Hs
(
    const fvMesh& mesh, const scalar p, const scalar T
) const
{
    return Cp_*T + EquationOfState::H(mesh, p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::S
(
    const fvMesh& mesh, const scalar p, const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::S(mesh, p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::dGdT
(
    const fvMesh& mesh, const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::dCpdT
(
    const fvMesh& mesh, const scalar p, const scalar T
) const
{
    return 0;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hConstThermo<EquationOfState>::operator+=
(
    const hConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator+
(
    const hConstThermo<EquationOfState>& ct1,
    const hConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return hConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cp_,
            ct1.Hf_
        );
    }
    else
    {
        return hConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hConstThermo<EquationOfState>& ct
)
{
    return hConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator==
(
    const hConstThermo<EquationOfState>& ct1,
    const hConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
