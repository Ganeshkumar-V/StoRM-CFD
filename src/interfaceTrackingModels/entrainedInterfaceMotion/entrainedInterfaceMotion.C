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
        dict
    ),
    bed_
    (
        pair_.phase1().mesh(),
        "bed",
        dict
    ),
    dmdt_
    (
      new volScalarField
      (
        IOobject("dmdt", pair_.phase1().mesh()),
        pair_.phase1().mesh(),
        dimensionedScalar("", dimVelocity/dimLength, 0.0)
      )
    ),
    MR(dict.get<scalar>("MR"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::entrainedInterfaceMotion::~entrainedInterfaceMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::interfaceTrackingModels::entrainedInterfaceMotion::correct()
{}

void Foam::interfaceTrackingModels::entrainedInterfaceMotion::regress
(
    const scalar fp,
    volScalarField& alpha
)
{
    // Pressure for burning Rate
    const volScalarField& p(pair_.phase1().thermo().p());

    // 1. Regress flame surface and get dmdt
    tmp<volScalarField> dmdtf = flame_.regressInterface(p, bed_.iNeighbours());

    // 2. Regress bed surface based on flame surface dmdt
    dmdt_ = bed_.regressInterface(p, dmdtf(), flame_.iNeighbours(), fp, MR);

    // 3. Update propellant alpha
    alpha = bed_.alpha();

    // Info << "Flame Volume: "
    //   << sum(flame_.alpha().internalField()*flame_.alpha().mesh().V())*10000 << endl;
    // Info << "bed Volume: "
    //   << sum(bed_.alpha().internalField()*bed_.alpha().mesh().V())*10000 << endl;
    //
    // const scalar dt = alpha.mesh().time().deltaTValue();
    // Info << "Flame: sum(dmdt*V*dt): "
    //   << sum(dmdtf().internalField()*flame_.alpha().mesh().V())*dt*10000 << endl;
    // Info << "Bed: sum(dmdt*V*dt): "
    //   << sum(dmdt_().internalField()*bed_.alpha().mesh().V())*dt*10000 << endl;

}

void Foam::interfaceTrackingModels::entrainedInterfaceMotion::store()
{
    // store old time step alpha
    bed_.store();
    flame_.store();
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

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::entrainedInterfaceMotion::dmdt() const
{
    return dmdt_;
}
// ************************************************************************* //
