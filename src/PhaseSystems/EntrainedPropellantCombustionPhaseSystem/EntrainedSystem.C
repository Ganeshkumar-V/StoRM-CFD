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
#include "EntrainedSystem.H"

namespace Foam
{
    // defineTypeNameAndDebug(EntrainedSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EntrainedSystem::EntrainedSystem
(
    const dictionary& dict
)
:
    rhoPropellant("rhop", dimDensity, 0),
    rhoP("rhoParticle", dimDensity, 0),
    alphaRhoAlp("alphaRhoAl", dimDensity, 0),
    alphaP("alphaP", dimless, 0),
    RT("RT", dimPressure/dimDensity, 0)
    {
      // Read Data from dictionary

        // MolecularWeights
        molecularWeights = dict.subDict("molecularWeight");
        MW.Al = molecularWeights.get<scalar>("Al");
        MW.Al2O3 = molecularWeights.get<scalar>("Al2O3");
        MW.H2O = molecularWeights.get<scalar>("H2O");
        MW.H2 = molecularWeights.get<scalar>("H2");

        // Active Aluminum content
        rhoPropellant = dimensionedScalar("", dimDensity,
                            dict.get<scalar>("propellantRho"));
        AAlC = dict.get<scalar>("activeAlContent");
        K = (MW.Al/MW.Al2O3)*(1/AAlC - 1);

        // equivalenceRatio
        phi = dict.get<scalar>("equivalenceRatio");
        zeta = 1.5/phi;

        // Chemical efficiency
        eta = dict.get<scalar>("chemicalEfficiency");

        // number of Moles
        n.Al2O3 = (K + eta/2);
        n.Al = (1 - eta);
        n.H2O = (zeta - 1.5*eta);
        n.H2 = 1.5*eta;

        // Molecular Weight of Propellant and Product gas
        MW.Prop = MW.Al + K*MW.Al2O3 + zeta*MW.H2O;
        MW.gas = (n.H2*MW.H2 + n.H2O*MW.H2O)/(n.H2 + n.H2O);

        // Adiabatic Flame Temperature
        Tad = dict.get<scalar>("Tad");

        // Find mass transfer rate factors
        mtf.particles = (n.Al2O3*MW.Al2O3 + n.Al*MW.Al)/MW.Prop;
        mtf.H2 = n.H2*MW.H2/MW.Prop;
        mtf.H2O = n.H2O*MW.H2O/MW.Prop;
        Info << "Mass Transfer Rate Factors: " << endl;
        Info << "       particles: " << mtf.particles << endl;
        Info << "             H2O: " << mtf.H2O << endl;
        Info << "              H2: " << mtf.H2 << endl;
        Info << "Molecular Weight-> Gas: " << MW.gas << endl;

        // Particle Phase density
        rhoP = dimensionedScalar("", dimDensity, (n.Al2O3*MW.Al2O3 + n.Al*MW.Al)
                      *(2700*3950)/(2700*n.Al2O3*MW.Al2O3 + 3950*n.Al*MW.Al));
        Info << "Particle Phase Density: " << rhoP << endl;

        // Reactant apparent particle density
        alphaRhoAlp = dimensionedScalar("", dimDensity, (MW.Al + K*MW.Al2O3)
                                        *rhoPropellant.value()/MW.Prop);

        alphaP = alphaRhoAlp/rhoP;

        // gas phase density
        RT = dimensionedScalar("RT", dimPressure/dimDensity, 8314.5*Tad/MW.gas);
    }
