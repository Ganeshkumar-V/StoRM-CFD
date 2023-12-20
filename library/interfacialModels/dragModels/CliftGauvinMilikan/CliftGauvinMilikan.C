/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "CliftGauvinMilikan.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(CliftGauvinMilikan, 0);
    addToRunTimeSelectionTable(particleDragModel, CliftGauvinMilikan, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::CliftGauvinMilikan::CliftGauvinMilikan
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    particleDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    R_
    (
      "R", dimVelocity*dimVelocity/dimTemperature,
      dict.get<scalar>("R")
    ),
    gamma_(dict.get<scalar>("gamma"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModels::CliftGauvinMilikan::~CliftGauvinMilikan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::CliftGauvinMilikan::CdRe() const
{
    const tmp<volScalarField> tRe(pair_.Re());
    const volScalarField& Re(tRe());

    const tmp<volScalarField> tT(pair_.continuous().thermo().T());
    const volScalarField& T(tT());

    volScalarField M(max(pair_.magUr()/sqrt(gamma_*R_*T), SMALL));

    const tmp<volScalarField> tKn(sqrt(constant::mathematical::pi)*sqrt(gamma_/2)*M/max(Re, SMALL));
    const volScalarField& Kn(tKn());

    const tmp<volScalarField> tFact(max(1.0 + Kn*(2.49 + 0.84*exp(-1.74/Kn)), SMALL));
    const volScalarField& Fact(tFact());
    Info << "Correction Factor found" << endl;

    return 24.0*(1.0 + 0.15*pow(Re, 0.687) + 0.0175*Re/(1.0 + 42500/max(pow(Re, 1.16), SMALL)))
          /Fact;
}


// ************************************************************************* //
