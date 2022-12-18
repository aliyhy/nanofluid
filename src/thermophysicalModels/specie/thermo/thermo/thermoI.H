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

#include "thermo.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type>::thermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::T
(
    scalar f,
    const fvMesh& mesh,
    scalar p,
    scalar T0,
    scalar (thermo<Thermo, Type>::*F)(const fvMesh&, const scalar, const scalar) const,
    scalar (thermo<Thermo, Type>::*dFdT)(const fvMesh&, const scalar, const scalar)
        const,
    scalar (thermo<Thermo, Type>::*limit)(const scalar) const
) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(mesh ,p, Test) - f)/(this->*dFdT)(mesh, p, Test));

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type>::thermo
(
    const word& name,
    const thermo& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::species::thermo<Thermo, Type>::heName()
{
    return Type<thermo<Thermo, Type>>::energyName();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Cpv(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::Cpv(*this, mesh, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::gamma(const fvMesh& mesh, const scalar p, const scalar T) const
{
    #ifdef __clang__
    volatile const scalar Cp = this->Cp(mesh, p, T);
    #else
    const scalar Cp = this->Cp(mesh, p, T);
    #endif

    return Cp/(Cp - this->CpMCv(mesh, p, T));
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::CpByCpv
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    return Type<thermo<Thermo, Type>>::CpByCpv(*this, mesh, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::HE(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::HE(*this, mesh, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::G(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Ha(mesh, p, T) - T*this->S(mesh, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::A(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Ea(mesh, p, T) - T*this->S(mesh, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::cp(const fvMesh& mesh,  const scalar p, const scalar T) const
{
    return this->Cp(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::ha(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Ha(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::hs(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Hs(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::hc() const
{
    return this->Hc()*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::s(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->S(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::he(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->HE(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::cv(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Cv(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::es(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Es(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::ea(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->Ea(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::g(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->G(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::a(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->A(mesh, p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::K(const fvMesh& mesh, const scalar p, const scalar T) const
{
    scalar arg = -this->Y()*this->G(mesh, Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Kp(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return K(mesh, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Kc(const fvMesh& mesh, const scalar p, const scalar T) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(mesh, p, T);
    }
    else
    {
        return Kp(mesh, p, T)*pow(Pstd/(RR*T), nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::Kx
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(mesh, p, T);
    }
    else
    {
        return Kp(mesh, p, T)*pow(Pstd/p, nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::Kn
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(mesh, p, T);
    }
    else
    {
        return Kp(mesh, p, T)*pow(n*Pstd/p, nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::THE
(
    const scalar he,
    const fvMesh& mesh,
    const scalar p,
    const scalar T0
) const
{
    return Type<thermo<Thermo, Type>>::THE(*this, he, mesh, p, T0);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::THs
(
    const scalar hs,
    const fvMesh& mesh,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hs,
      mesh,
        p,
        T0,
        &thermo<Thermo, Type>::Hs,
        &thermo<Thermo, Type>::Cp,
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::THa
(
    const scalar ha,
    const fvMesh& mesh,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ha,
      mesh,
        p,
        T0,
        &thermo<Thermo, Type>::Ha,
        &thermo<Thermo, Type>::Cp,
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::TEs
(
    const scalar es,
    const fvMesh& mesh,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        es,
      mesh,
        p,
        T0,
        &thermo<Thermo, Type>::Es,
        &thermo<Thermo, Type>::Cv, //Cv
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::TEa
(
    const scalar ea,
    const fvMesh& mesh,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ea,
      mesh,
        p,
        T0,
        &thermo<Thermo, Type>::Ea,
        &thermo<Thermo, Type>::Cv,//Cv
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::dKcdTbyKc
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return -this->dGdT(mesh, Pstd, T)*this->Y()/RR;
    }
    else
    {
        return -(nm/T + this->dGdT(mesh, Pstd, T)*this->Y()/RR);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::dcpdT(const fvMesh& mesh, const scalar p, const scalar T) const
{
    return this->dCpdT(mesh, p, T)*this->W();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline void Foam::species::thermo<Thermo, Type>::operator+=
(
    const thermo<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::species::thermo<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator+
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return thermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator*
(
    const scalar s,
    const thermo<Thermo, Type>& st
)
{
    return thermo<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator==
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return thermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) == static_cast<const Thermo&>(st2)
    );
}


// ************************************************************************* //
