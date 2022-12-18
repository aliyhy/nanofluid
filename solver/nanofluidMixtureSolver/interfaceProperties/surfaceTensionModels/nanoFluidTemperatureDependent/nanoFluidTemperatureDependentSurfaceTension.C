/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "nanoFluidTemperatureDependentSurfaceTension.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(nanoFluidTemperatureDependent, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        nanoFluidTemperatureDependent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::nanoFluidTemperatureDependent::nanoFluidTemperatureDependent
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    TName_(dict.getOrDefault<word>("T", "T")),
    sigma_(Function1<scalar>::New("sigma", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::nanoFluidTemperatureDependent::~nanoFluidTemperatureDependent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::nanoFluidTemperatureDependent::sigma() const
{
  auto tsigma = tmp<volScalarField>::New
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimSigma
    );
    auto& sigma = tsigma.ref();

    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    scalarIOField phiNP_ = 
      mesh_.objectRegistry::thisDb().lookupObjectRef<scalarIOField>("phiNP");
    scalar a_ = 7.673 / pow(10,7);
    scalar b_ = -7.773/pow(10,3);
    scalar c_ = phiNP_[0];
    
    sigma.field() =
      sigma_->value(T.field())
      *(1. - b_*Foam::log(c_/a_ + 1.0)/Foam::log(2.71828182845));

    volScalarField::Boundary& sigmaf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaf, patchi)
    {
        sigmaf[patchi] = sigma_->value(TBf[patchi]);
    }
    return tsigma;
}



bool Foam::surfaceTensionModels::nanoFluidTemperatureDependent::readDict
(
    const dictionary& dict
)
{
    const dictionary& sigmaDict = surfaceTensionModel::sigmaDict(dict);

    TName_ = sigmaDict.getOrDefault<word>("T", "T");
    sigma_ = Function1<scalar>::New("sigma", sigmaDict);

    return true;
}


bool Foam::surfaceTensionModels::nanoFluidTemperatureDependent::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_() << token::END_STATEMENT << nl;
        return os.good();
    }

    return false;
}


// ************************************************************************* //
