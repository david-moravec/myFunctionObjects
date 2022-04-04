/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "dissipation.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dissipation, 0);
    addToRunTimeSelectionTable(functionObject, dissipation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dissipation::dissipation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, false),
    phaseName_(word::null),
    volScalarField dissipation_
    (
        IOobject
        (
            typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        ),
        mesh_,
        dimensionedScalar(dimset(0 2 -2 0 0 0 0), 0), 
    )

{
    read(dict);
    resetLocalObjectName(IOobject::groupName(type(), phaseName_));


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::dissipation::~dissipation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dissipation::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    //phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


bool Foam::functionObjects::dissipation::execute()
{
    if(!isA<fvMesh>(mesh_))
    {
        WarningIn
        (
            "dissipation needs fvMesh, deactivating"
        );
        return true;
    }

    std::cout << "ekk\n";
    volVectorField U = mesh_.lookupObject<volVectorField>("U");

    volScalarField nut = mesh_.lookupObject<volScalarField>("nut");
    volScalarField nu = mesh_.lookupObject<volScalarField>("nu");

    volSymmTensorField shearRate = dev(twoSymm(fvc::grad(U)));

    dissipation_ = nu * shearRate;
    
    store("R", nut * shearRate);
    store("Tau", nu * shearRate);

    return true;
}


bool Foam::functionObjects::dissipation::write()
{
    return true;
}


// ************************************************************************* //
