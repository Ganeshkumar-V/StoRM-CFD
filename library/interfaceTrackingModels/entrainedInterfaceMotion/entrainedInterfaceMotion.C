/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

#include "entrainedInterfaceMotion.H"
#include "phasePair.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "processorFvPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTrackingModels
{
    defineTypeNameAndDebug(entrainedInterfaceMotion, 0);
    addToRunTimeSelectionTable(interfaceTrackingModel, entrainedInterfaceMotion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::entrainedInterfaceMotion::entrainedInterfaceMotion
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceTrackingModel(dict, pair),
    n(dict.get<scalar>("n")),
    f(dict.get<scalar>("f")),
    a("", pow(dimLength*dimTime*dimTime/dimMass, n)*dimVelocity, dict.get<scalar>("a")),
    flame_
    (
        pair_.phase1().mesh(),
        "flame",
        n, f, a
    ),
    bed_
    (
        pair_.phase1().mesh(),
        "bed",
        n, f, a
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::entrainedInterfaceMotion::~entrainedInterfaceMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::interfaceTrackingModels::entrainedInterfaceMotion::correct()
{}

void Foam::interfaceTrackingModels::entrainedInterfaceMotion::regress
(
    volScalarField& alpha,
    const volScalarField& alphaOld
)
{
    // Pressure for burning Rate
    const volScalarField& p(pair_.phase1().thermo().p());

    // 1. Regress flame surface and get dmdt
    tmp<volScalarField> dmdtf = flame_->regressInterface(p, bed.iNeighbours());

    // 2. Regress bed surface based on flame surface dmdt
    tmp<volScalarField> dmdtb = bed_->regressInterface(p, dmdtf(), flame.iNeighbours());

}

void Foam::interfaceTrackingModels::entrainedInterfaceMotion::store()
{
    // store old time step alpha
    alphaOld_ = alpha_;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::entrainedInterfaceMotion::rb() const
{
    return bed_.rb();
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::entrainedInterfaceMotion::As() const
{
    return bed_.As();
}
// ************************************************************************* //
