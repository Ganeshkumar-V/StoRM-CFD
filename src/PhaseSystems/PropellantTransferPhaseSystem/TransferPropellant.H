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
namespace Foam
{

  struct Mfactors
  {
    scalar particles;
    scalar gas;
  };

/*---------------------------------------------------------------------------*\
                            Class Chemical Efficiency
\*---------------------------------------------------------------------------*/

class TransferPropellant
{
protected:

  // Protected data

    // Propellant density and Active Aluminium content
    dimensionedScalar rhoPropellant;

    // Adiabatic Flame Temperature
    scalar Tad;

    // Mass Transfer Rate factors
    Mfactors mtf;

    // Particle Phase density
    dimensionedScalar rhoP;

    // Reactant particle phase apparent density
    dimensionedScalar alphaP;

    // gas phase RT
    dimensionedScalar RT;

    // Temperature Dependence
    scalar A_;
    scalar B_;
    scalar C_;

public:

  // Constructors

    //- Construct from dictionary
    TransferPropellant(const fvMesh& mesh, const dictionary& dict)
    :
      rhoPropellant("rhop", dimDensity, dict.get<scalar>("propellantRho")),
      rhoP("rhoParticle", dimDensity, 2700),
      alphaP("alphaP", dimless, dict.get<scalar>("alphaPReactant")),
      RT("RT", dimPressure/dimDensity, 0),
      A_(dict.get<scalar>("A")),
      B_(dict.get<scalar>("B")),
      C_(dict.get<scalar>("C"))
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

  // Methods

    // Public functions
    struct Mfactors massTransfer() const
    {
        return mtf;
    }

    void gasVelocity
    (
        const volVectorField& nHat,
        const volScalarField& P,
        const volScalarField& rb,
        volVectorField& Ug
    )
    {
        // Gas phase density
        const tmp<volScalarField> trhoG(P/RT);
        const volScalarField& rhoG(trhoG());

        // Gas Velocity
        Ug = rb*nHat*(1 - (rhoPropellant - alphaP*rhoP)/(rhoG*(1 - alphaP)));
    }

    // Enthalpy
    void getAdiabaticTemperature
    (
        const volScalarField& P,
        volScalarField& Tad
    )
    {
        forAll(Tad, i)
        {
           Tad[i] = A_*pow(P[i]/1e6, B_) + C_;
        }
    }

};

}
