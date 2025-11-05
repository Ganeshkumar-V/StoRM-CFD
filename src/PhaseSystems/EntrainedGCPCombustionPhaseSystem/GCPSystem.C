/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation and Ganeshkumar V, IIT Gandhinagar
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
#include "GCPSystem.H"

namespace Foam
{
    // defineTypeNameAndDebug(GCPSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GCPSystem::GCPSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    rhoP("rhoParticle", dimDensity, 2700),
    alphaP("alphaP", dimless, dict.get<scalar>("alphaPReactant")),
    RT("RT", dimPressure/dimDensity, 0)
    {
      // Read Data from dictionary
        dictionary gasDict = mesh.lookupObject<dictionary>("thermophysicalProperties.gas");

        // Adiabatic Flame Temperature
        Tad = dict.get<scalar>("Tad");

        // Find mass transfer rate factors
        mtf.particles = dict.subDict("massTransferRates").get<scalar>("particles");
        mtf.gas = dict.subDict("massTransferRates").get<scalar>("gas");

        Info << "Mass Transfer Rate Factors: " << endl;
        Info << "       particles: " << mtf.particles << endl;
        Info << "             gas: " << mtf.gas << endl;

        // gas phase density
        scalar MWGas = gasDict.subDict("mixture").subDict("specie").get<scalar>("molWeight");
        RT = dimensionedScalar("RT", dimPressure/dimDensity, 8314.5*Tad/MWGas);
    }
