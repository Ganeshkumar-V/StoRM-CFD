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

#include "Brenner.H"
#include "mathematicalConstants.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sharpInterfaceHeatTransferModels
{
    defineTypeNameAndDebug(Brenner, 0);
    addToRunTimeSelectionTable(sharpInterfaceHeatTransferModel, Brenner, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::Brenner::Brenner
(
    const dictionary& dict,
    const phasePair& pair
)
:
    sharpInterfaceHeatTransferModel(dict, pair),
    cutoff(dict.get<scalar>("cutoff")),
    d_("d", dimLength, dict.get<scalar>("d")),
    t_("t", dimLength, dict.get<scalar>("t")),
    dV_("dV", dimLength, dict.get<scalar>("dV")),
    Xcorr_("Xcorr", dimless, 0.0)
{
    Xcorr_ = (1.0/3.0)*(sqr(dV_)/sqr(d_))*(2.0 + d_/t_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::Brenner::~Brenner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModels::Brenner::K(const scalar residualAlpha) const
{
  // volScalarField Re(pair_.Re());
  // forAll(Re, i)
  // {
  //   if (Re[i] < 0)
  //   {
  //     Info << "Cell: " << i << " Re[i]: " << Re[i] << endl;
  //   }
  // }

    const tmp<volScalarField> tPe(max(pair_.Re()*pair_.Pr()*d_/dV_, SMALL));
    const volScalarField& Pe(tPe());

    volScalarField Nu
    (
        4.0/constant::mathematical::pi
        + (2.0/sqr(constant::mathematical::pi))*Pe
        + 60.0/(27.0*pow(constant::mathematical::pi, 3.0))*Pe*log(Pe)
    );

    return
        6.0
       *pair_.dispersed()*pos(pair_.dispersed() - cutoff)
       *pair_.continuous().kappa()
       *Nu
       *Xcorr_
       /sqr(pair_.dispersed().d());
}

// ************************************************************************* //
