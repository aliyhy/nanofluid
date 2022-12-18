/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Class
    Foam::nanoFluidPolynomialTransport

Group
    grpSpecieTransport

Description
    Transport package using polynomial functions for \c mu and \c kappa.

Usage

    \table
        Property        | Description
        muCoeffs<8>     | Dynamic viscosity polynomial coefficients
        kappaCoeffs<8>  | Thermal conductivity polynomial coefficients
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        muCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
        kappaCoeffs<8>  ( 2000 -0.15 0.023 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expressions are evaluated as so:

        \f[
            \mu    = 1000 - 0.05 T + 0.003 T^2
        \f]

        \f[
            \kappa = 2000 - 0.15 T + 0.023 T^2
        \f]

Note
    - Dynamic viscosity polynomial coefficients evaluate to an expression in
      [Pa.s], but internally uses [Pa.s/kmol].
    - Thermal conductivity polynomial coefficients evaluate to an expression in
      [W/m/K], but internally uses [W/m/K/kmol].

SourceFiles
    nanoFluidPolynomialTransportI.H
    nanoFluidPolynomialTransport.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef nanoFluidPolynomialTransport_H
#define nanoFluidPolynomialTransport_H

#include "Polynomial.H"
#include "volFields.H"
#include "scalarField.H"
#include "fvcVolumeIntegrate.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo, int PolySize> class nanoFluidPolynomialTransport;

template<class Thermo, int PolySize>
inline nanoFluidPolynomialTransport<Thermo, PolySize> operator+
(
    const nanoFluidPolynomialTransport<Thermo, PolySize>&,
    const nanoFluidPolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline nanoFluidPolynomialTransport<Thermo, PolySize> operator*
(
    const scalar,
    const nanoFluidPolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const nanoFluidPolynomialTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class nanoFluidPolynomialTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class nanoFluidPolynomialTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity polynomial coefficients
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity polynomial coefficients
        Polynomial<PolySize> kappaCoeffs_;

        //- Specific heat at constant pressure polynomial coeffs
        Polynomial<PolySize> CpCoeffs_;

        //- density polynomial coeffs
        Polynomial<PolySize> rhoCoeffs_;


        //- Nano particle concentrasion

        //- Nano particle diameter
        scalar dP_;

        //- Fluid relative molecular weight
        scalar molWeight_;

        //- Fluid density at T = 293 K
        scalar rho0_;

        //- Nano particle Thermal conductivelty
        Polynomial<PolySize> kappaPCoeffs_;

        //- Nano particle shape factor
        scalar nP_;

 

    // Private Member Functions

        //- Coeffs name. Eg, "muCoeffs<10>"
        inline static word coeffsName(const char* name)
        {
            return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        }

        //- Construct from components
        inline nanoFluidPolynomialTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly,
            const Polynomial<PolySize>& CpPoly,
            const Polynomial<PolySize>& rhoPoly,
            const scalar dP,
            const scalar molWeight,
            const scalar rho0,
            const Polynomial<PolySize>& kappaPPoly,
            const scalar nP
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct as named copy
        inline nanoFluidPolynomialTransport(const word&, const nanoFluidPolynomialTransport&);

        //- Construct from dictionary
        explicit nanoFluidPolynomialTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<nanoFluidPolynomialTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<nanoFluidPolynomialTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CorcioneHamilton<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity of nanofluid[kg/ms]
        inline scalar mu(const fvMesh& mesh, const scalar p,const scalar T) const;

        //- Thermal conductivity of nanofluid[W/mK]
        inline scalar kappa(const fvMesh& mesh, const scalar p,const scalar T) const;

        //- Thermal diffusivity of enthalpy of nanofluid[kg/ms]
        inline scalar alphah(const fvMesh& mesh, const scalar p,const scalar T) const;

        //- Dynamic viscosity of basefluid[kg/ms]
        inline scalar muBF(const fvMesh& mesh, const scalar p,const scalar T) const;

        //- Thermal conductivity of base fluid[W/mK]
        inline scalar kappaBF(const fvMesh& mesh, const scalar p,const scalar T) const;

        //- Thermal diffusivity of enthalpy of Base fluid[kg/ms]
        inline scalar alphahBF(const fvMesh& mesh, const scalar p,const scalar T) const;

 
        // Species diffusivity
        //inline scalar D(const fvMesh& mesh, const scalar p,const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const nanoFluidPolynomialTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend nanoFluidPolynomialTransport operator+ <Thermo, PolySize>
        (
            const nanoFluidPolynomialTransport&,
            const nanoFluidPolynomialTransport&
        );

        friend nanoFluidPolynomialTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const nanoFluidPolynomialTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const nanoFluidPolynomialTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nanoFluidPolynomialTransportI.H"

#ifdef NoRepository
    #include "nanoFluidPolynomialTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
