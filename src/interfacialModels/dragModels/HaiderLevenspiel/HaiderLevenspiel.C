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

#include "HaiderLevenspiel.H"
#include "mathematicalConstants.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(HaiderLevenspiel, 0);
    addToRunTimeSelectionTable(particleDragModel, HaiderLevenspiel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::HaiderLevenspiel::HaiderLevenspiel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    particleDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    phi_("phi", dimless, dict.get<scalar>("phi")),
    A_("A", dimless, 0),
    B_("B", dimless, 0),
    C_("C", dimless, 0),
    D_("D", dimless, 0)
{
    A_ = exp(2.3288 - 6.4581*phi_ + 2.4486*pow(phi_, 2));
    B_ = 0.0964 + 0.5565*phi_;
    C_ = exp(4.905 - 13.8944*phi_ + 18.4222*pow(phi_, 2) - 10.2599*pow(phi_, 3));
    D_ = exp(1.4681 + 12.2584*phi_ - 20.7322*pow(phi_, 2) + 15.8855*pow(phi_, 3));

    // Print Coefficients
    Info << "Haider-Levenspiel Coefficients: " << endl;
    Info << "         A = " << A_.value() << endl;
    Info << "         B = " << B_.value() << endl;
    Info << "         C = " << C_.value() << endl;
    Info << "         D = " << D_.value() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModels::HaiderLevenspiel::~HaiderLevenspiel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::particleDragModels::HaiderLevenspiel::CdRe() const
{
    const tmp<volScalarField> tRe(max(pair_.Re(), SMALL));
    const volScalarField& Re(tRe());

    return   24*(1.0 + A_*pow(Re, B_)) + C_*sqr(Re)/(D_ + Re);
}

// ************************************************************************* //
