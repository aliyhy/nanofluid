/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermophysicalProperties>
inline const ThermophysicalProperties&
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>
::properties() const
{
    return *propertiesPtr_;
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::W() const
{
    return propertiesPtr_->W();
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::limit
(
    const scalar T
) const
{
    return propertiesPtr_->limit(T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::rho
(
    const fvMesh& mesh, 
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->rho(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::psi
(
 const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->psi(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::CpMCv
(
 const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->CpMCv(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Cp
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Cp(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Ha
(
 const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Ha(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Hs
(
 const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Hs(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Hc() const
{
    return propertiesPtr_->Hc();
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Cv
(
 const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Cp(mesh, p, T) - CpMCv(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Es
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Hs(mesh, p, T) - p/propertiesPtr_->rho(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Ea
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Ha(mesh, p, T) - p/propertiesPtr_->rho(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::S
(
 const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->S(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::mu
(
    const fvMesh& mesh, 
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->mu(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::kappa
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->kappa(mesh, p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::alphah
(
    const fvMesh& mesh,
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->alphah(mesh, p, T);
}


// ************************************************************************* //
