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
    Foam::icoNanoFluidPolynomial

Group
    grpSpecieEquationOfState

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density.

Usage
    \table
        Property     | Description
        rhoCoeffs<8> | Density polynomial coefficients
    \endtable

    Example of the specification of the equation of state:
    \verbatim
    equationOfState
    {
        rhoCoeffs<8>    ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            \rho = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    Input in [kg/m3], but internally uses [kg/m3/kmol].

SourceFiles
    icoNanoFluidPolynomialI.H
    icoNanoFluidPolynomial.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef icoNanoFluidPolynomial_H
#define icoNanoFluidPolynomial_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "Polynomial.H"
#include "scalarField.H"
#include "fvcVolumeIntegrate.H"
#include "fvMesh.H"
#include "scalarField.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Specie, int PolySize>
class icoNanoFluidPolynomial;

template<class Specie, int PolySize>
icoNanoFluidPolynomial<Specie, PolySize> operator+
(
    const icoNanoFluidPolynomial<Specie, PolySize>&,
    const icoNanoFluidPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoNanoFluidPolynomial<Specie, PolySize> operator*
(
    const scalar,
    const icoNanoFluidPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoNanoFluidPolynomial<Specie, PolySize> operator==
(
    const icoNanoFluidPolynomial<Specie, PolySize>&,
    const icoNanoFluidPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const icoNanoFluidPolynomial<Specie, PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class icoNanoFluidPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Specie, int PolySize=8>
class icoNanoFluidPolynomial
:
    public Specie
{
    // Private Data

        //- Density polynomial coefficients
        Polynomial<PolySize> rhoCoeffs_;

        //- Nano particle concentrasion
       
        //- initial liquid volume to the total volume

        //- Nano particle density
        scalar rhoP_;




    // Private Member Functions

        //- Coeffs name. Eg, "rhoCoeffs<10>"
        inline static word coeffsName(const char* name)
        {
            return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        }


public:

    // Generated Methods: copy construct, copy assignment

    // Constructors

        //- Construct from components
        inline icoNanoFluidPolynomial
        (
            const Specie& sp,
            const Polynomial<PolySize>& rhoPoly,
            const scalar rhoP
        );

        //- Construct from dictionary
        explicit icoNanoFluidPolynomial(const dictionary& dict);

        //- Construct as named copy
        inline icoNanoFluidPolynomial(const word& name, const icoNanoFluidPolynomial&);

        //- Construct and return a clone
        inline autoPtr<icoNanoFluidPolynomial> clone() const;

        // Selector from dictionary
        inline static autoPtr<icoNanoFluidPolynomial> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "icoNanoFluidPolynomial<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density of nanofluid[kg/m^3]
            inline scalar rho(const fvMesh& mesh, scalar p, scalar T)  const;

            //- Return density of base fluid[kg/m^3]
            inline scalar rhoBF(const fvMesh& mesh, scalar p, scalar T) const;
 
            //- Return enthalpy departure [J/kg]
            inline scalar H(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]

           inline scalar Cp(const fvMesh& mesh, scalar p, scalar T) const;

           inline scalar phiNP(const fvMesh& mesh) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(const fvMesh& mesh, scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(const fvMesh& mesh, scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(const fvMesh& mesh, scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(const fvMesh& mesh, scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const icoNanoFluidPolynomial&);
        inline void operator*=(const scalar);


    // Friend Operators

        friend icoNanoFluidPolynomial operator+ <Specie, PolySize>
        (
            const icoNanoFluidPolynomial&,
            const icoNanoFluidPolynomial&
        );

        friend icoNanoFluidPolynomial operator* <Specie, PolySize>
        (
            const scalar s,
            const icoNanoFluidPolynomial&
        );

        friend icoNanoFluidPolynomial operator== <Specie, PolySize>
        (
            const icoNanoFluidPolynomial&,
            const icoNanoFluidPolynomial&
        );


    // IOstream Operators

        friend Ostream& operator<< <Specie, PolySize>
        (
            Ostream&,
            const icoNanoFluidPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIcoNanoFluidPolynomial(PolySize)                                            \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    icoNanoFluidPolynomial<Specie, PolySize>,                                           \
    "icoNanoFluidPolynomial<"#PolySize">",                                              \
    0                                                                          \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoNanoFluidPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "icoNanoFluidPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
