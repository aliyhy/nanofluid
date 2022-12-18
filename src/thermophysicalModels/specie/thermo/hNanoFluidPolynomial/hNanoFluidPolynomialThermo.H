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
    Foam::hNanoFluidPolynomialThermo

Group
    grpSpecieThermo

Description
    Thermodynamics package templated on the equation of state, using polynomial
    functions for \c cp, \c h and \c s.

    Polynomials for \c h and \c s derived from \c cp.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        CpCoeffs<8>  | Specific heat at constant pressure polynomial coeffs
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        CpCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            Cp = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    - Heat of formation is inputted in [J/kg], but internally uses [J/kmol]
    - Standard entropy is inputted in [J/kg/K], but internally uses [J/kmol/K]
    - Specific heat at constant pressure polynomial coefficients evaluate to an
      expression in [J/(kg.K)].

SourceFiles
    hNanoFluidPolynomialThermoI.H
    hNanoFluidPolynomialThermo.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef hNanoFluidPolynomialThermo_H
#define hNanoFluidPolynomialThermo_H

#include "scalar.H"
#include "Polynomial.H"

#include "fvMesh.H"
#include "volFields.H"
#include "Polynomial.H"
#include "scalarField.H"
#include "fvcVolumeIntegrate.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class EquationOfState, int PolySize>
class hNanoFluidPolynomialThermo;

template<class EquationOfState, int PolySize>
inline hNanoFluidPolynomialThermo<EquationOfState, PolySize> operator+
(
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>&,
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hNanoFluidPolynomialThermo<EquationOfState, PolySize> operator*
(
    const scalar,
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hNanoFluidPolynomialThermo<EquationOfState, PolySize> operator==
(
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>&,
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const hNanoFluidPolynomialThermo<EquationOfState, PolySize>&
);


/*---------------------------------------------------------------------------*\
                      Class hNanoFluidPolynomialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState, int PolySize=8>
class hNanoFluidPolynomialThermo
:
    public EquationOfState
{
    // Private Data

        //- Heat of formation
        scalar Hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant pressure polynomial coeffs
        Polynomial<PolySize> CpCoeffs_;

        //- Enthalpy polynomial coeffs - derived from cp [J/kg]
        //  NOTE: relative to Tstd
        typename Polynomial<PolySize>::intPolyType hCoeffs_;

        //- Entropy - derived from Cp [J/(kg.K)] - relative to Tstd
        Polynomial<PolySize> sCoeffs_;

        //- Density polynomial coefficients
        Polynomial<PolySize> rhoCoeffs_;


        //- Nano particle density
        scalar rhoP_;

        //- Nano particle specific heta capacity
        scalar CpP_;


    // Private Member Functions

        //- Coeffs name. Eg, "CpCoeffs<10>"
        inline static word coeffsName(const char* name)
        {
            return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        }

        //- Construct from components
        inline hNanoFluidPolynomialThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
            const Polynomial<PolySize>& CpCoeffs,
            const typename Polynomial<PolySize>::intPolyType& hCoeffs,
            const Polynomial<PolySize>& sCoeffs,
            const Polynomial<PolySize>& rhoPoly,
            const scalar rhoP,
            const scalar CpP
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct from dictionary
        hNanoFluidPolynomialThermo(const dictionary& dict);

        //- Construct as a named copy
        inline hNanoFluidPolynomialThermo(const word&, const hNanoFluidPolynomialThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hNanoFluidPolynomial<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


        // Fundamental properties

            //- Heat capacity of nanofluid at constant pressure [J/(kg K)]
            inline scalar Cp(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Heat capacity of basefluid at constant pressure [J/(kg K)]
            inline scalar CpBF(const fvMesh& mesh, const scalar p, const scalar T) const;

            inline scalar phiNP(const fvMesh& mesh) const; 

 
            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const fvMesh& mesh, const scalar p, const scalar T) const;


            #include "HtoEthermo.H"


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const fvMesh& mesh, const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const fvMesh& mesh, const scalar p, const scalar T) const;



        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const hNanoFluidPolynomialThermo&);
        inline void operator*=(const scalar);


    // Friend Operators

        friend hNanoFluidPolynomialThermo operator+ <EquationOfState, PolySize>
        (
            const hNanoFluidPolynomialThermo&,
            const hNanoFluidPolynomialThermo&
        );

        friend hNanoFluidPolynomialThermo operator* <EquationOfState, PolySize>
        (
            const scalar,
            const hNanoFluidPolynomialThermo&
        );

        friend hNanoFluidPolynomialThermo operator== <EquationOfState, PolySize>
        (
            const hNanoFluidPolynomialThermo&,
            const hNanoFluidPolynomialThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState, PolySize>
        (
            Ostream&,
            const hNanoFluidPolynomialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hNanoFluidPolynomialThermoI.H"

#ifdef NoRepository
    #include "hNanoFluidPolynomialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
