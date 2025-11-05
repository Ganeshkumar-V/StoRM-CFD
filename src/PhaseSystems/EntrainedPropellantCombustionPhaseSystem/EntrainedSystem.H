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

#ifndef EntrainedSystem_H
#define EntrainedSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

namespace Foam
{
  namespace Entrainment
  {
    struct molecularWeight
    {
      scalar Al;
      scalar Al2O3;
      scalar H2O;
      scalar H2;
      scalar Prop;
      scalar gas;
    };

    struct moles
    {
      scalar Al;
      scalar Al2O3;
      scalar H2O;
      scalar H2;
    };

    struct factors
    {
      scalar particles;
      scalar H2;
      scalar H2O;
      scalar gas;
    };
  }
/*---------------------------------------------------------------------------*\
                            Class Entrained System
\*---------------------------------------------------------------------------*/

class EntrainedSystem
{
protected:

  // Protected data

    // Molecular Weights of all the species
    dictionary molecularWeights;

    Entrainment::molecularWeight MW;

    Entrainment::moles n;

    // Propellant density and Active Aluminium content
    dimensionedScalar rhoPropellant;
    scalar AAlC;
    scalar K;

    // Equivalence Ratio
    scalar phi;
    scalar zeta;

    // Adiabatic Flame Temperature
    scalar Tad;

    // Chemical efficiency
    scalar eta;

    // Mass Transfer Rate factors
    Entrainment::factors mtf;

    // Particle Phase density
    dimensionedScalar rhoP;

    // Reactant particle phase apparent density
    dimensionedScalar alphaRhoAlp;
    dimensionedScalar alphaP;

    // gas phase RT
    dimensionedScalar RT;

public:

  //- Runtime type information
  // TypeName("EntrainedSystem");

  // Constructors

    //- Construct from dictionary
    EntrainedSystem(const dictionary& dict);

  // Methods

    // Public functions
    struct Entrainment::factors massTransfer() const
    {
        return mtf;
    }

    void gasVelocity
    (
        const volScalarField& P,
        const volScalarField& dmdt,
        const volScalarField& rb,
        const volScalarField& As,
        volVectorField& Ug
    )
    {
        // Gas phase density
        const tmp<volScalarField> trhoG(P/RT);
        const volScalarField& rhoG(trhoG());

        // Gas Velocity
        forAll(Ug, i)
        {
            if(As[i] != 0)
            {
                Ug[i] = ((rb[i] - (1 - mtf.particles)*dmdt[i]/
                          ((1 - alphaP)*As[i]*rhoG[i]))*vector(1, 0, 0)).value();
            }
        }
    }

    //- Destructor
    // virtual ~EntrainedSystem();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
