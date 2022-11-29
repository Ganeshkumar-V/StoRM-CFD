/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

#include "PropellantRegressionPhaseSystem.H"
#include "interfaceTrackingModel.H"
#include "fvmSup.H"
#include "phaseSystem.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::rDmdt
(
    const phasePairKey& key
) const
{
    if (!rDmdt_.found(key))
    {
        return phaseSystem::dmdt(key);
    }

    const scalar rDmdtSign(Pair<word>::compare(rDmdt_.find(key).key(), key));
    return rDmdtSign**rDmdt_[key];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::PropellantRegressionPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    molecularWeights_(this->subDict("molecularWeight")),
    eqR_(this->template get<scalar>("equivalenceRatio")),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    )
{
    this->generatePairsAndSubModels
    (
        "interfaceTracking",
        interfaceTrackingModels_
    );

    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        this->rDmdt_.set
        (
            interfaceTrackingModelIter.key(),
            phaseSystem::dmdt(interfaceTrackingModelIter.key()).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::
~PropellantRegressionPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    return BasePhaseSystem::dmdt(key) + rDmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    // Fill the mass transfer rates with zero
    this->fillFields("dmdt", dimDensity/dimTime, dmdts);

    forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
        const volScalarField& rDmdt = *rDmdtIter();
        // Add massTransfer Rate to the Gas Phase only
        this->addField(pair.phase2(), "dmdt", rDmdt, dmdts);
    }
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    // (No Species Present)

    return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
  autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
          BasePhaseSystem::heatTransfer();

  // phaseSystem::heatTransferTable& eqns = eqnsPtr();

  forAllConstIter
  (
      interfaceTrackingModelTable,
      interfaceTrackingModels_,
      interfaceTrackingModelIter
  )
  {}

  return eqnsPtr;
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::solve()
{
  // Regress Propellant surface (Manipulate propellant volume fraction)
  forAllIter
  (
    interfaceTrackingModelTable,
    interfaceTrackingModels_,
    interfaceTrackingModelIter
  )
  {
    word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
    volScalarField& alpha = this->db().template lookupObjectRef<volScalarField>(propellant);
    interfaceTrackingModelIter()->regress(alpha);
  }

  // Solve other phase volume fraction equations
  BasePhaseSystem::solve();
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    //- Finds burning rate (rb = aP^n)
    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        interfaceTrackingModelIter()->correct();
        *rDmdt_[interfaceTrackingModelIter.key()]
                  = dimensionedScalar(dimDensity/dimTime);
    }

    //- return burning Rate, As and propellant density -> find mdot of propellant
    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        *rDmdt_[interfaceTrackingModelIter.key()]
              = interfaceTrackingModelIter()->rb()
                *interfaceTrackingModelIter()->As()
                *interfaceTrackingModelIter()->rhop();
    }

}

template<class BasePhaseSystem>
bool Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
