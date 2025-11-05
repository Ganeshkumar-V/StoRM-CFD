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

#include "HolzerSommerfeld.H"
#include "mathematicalConstants.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(HolzerSommerfeld, 0);
    addToRunTimeSelectionTable(particleDragModel, HolzerSommerfeld, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::HolzerSommerfeld::HolzerSommerfeld
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    particleDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    phi_("phi", dimless, dict.get<scalar>("phi")),
    d_("d", dimLength, dict.get<scalar>("d")),
    t_("t", dimLength, dict.get<scalar>("t")),
    dV_("dV", dimLength, dict.get<scalar>("dV"))
{
      dimensionedScalar Across = constant::mathematical::pi*sqr(dV_)/4.0;
      dimensionedScalar Asurf = constant::mathematical::pi*sqr(d_)/2.0
                                + constant::mathematical::pi*d_*t_;

      // Length wise falling disks
      dimensionedScalar Aproj = d_*t_;
      dimensionedScalar AprojPerp = constant::mathematical::pi*sqr(d_)/4.0;
      length_.phiPara = Across/(0.5*Asurf - AprojPerp);
      length_.phiPerp = Across/Aproj;
      length_.A = 8.0/sqrt(length_.phiPara) + 16.0/sqrt(phi_);
      length_.B = 3.0/pow(phi_, 3.0/4.0);
      length_.C = pow(0.4210, 0.4*pow((-log(phi_)), 0.2))/length_.phiPerp;

      // Cross wise falling disks
      Aproj = constant::mathematical::pi*sqr(d_)/4.0;
      AprojPerp = d_*t_;
      cross_.phiPara = Across/(0.5*Asurf - AprojPerp);
      cross_.phiPerp = Across/Aproj;
      cross_.A = 8.0/sqrt(cross_.phiPara) + 16.0/sqrt(phi_);
      cross_.B = 3.0/pow(phi_, 3.0/4.0);
      cross_.C = pow(0.4210, 0.4*pow((-log(phi_)), 0.2))/cross_.phiPerp;

      // Print Coefficients
      Info << "Holzer-Sommerfeld Sphericity: " << endl;
      Info << "     Cross-wise falling disks ->  " << endl;
      Info << "         phi|| = " << cross_.phiPara.value() << endl;
      Info << "         phi|_ = " << cross_.phiPerp.value() << endl;
      Info << "     length-wise falling disks ->  " << endl;
      Info << "         phi|| = " << length_.phiPara.value() << endl;
      Info << "         phi|_ = " << length_.phiPerp.value() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModels::HolzerSommerfeld::~HolzerSommerfeld()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::HolzerSommerfeld::CdRe() const
{
    const tmp<volScalarField> tRe(max(pair_.Re(), SMALL));
    const volScalarField& Re(tRe());

    return 0.5*(length_.A + cross_.A)
           + 0.5*(length_.B + cross_.B)*sqrt(Re)
           + 0.5*(length_.C + cross_.C)*Re;
}

// ************************************************************************* //
