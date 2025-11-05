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

#include "EntrainedGCPCombustionPhaseSystem.H"
#include "interfaceTrackingModel.H"
#include "fvmSup.H"
#include "fvmDdt.H"
#include "phaseSystem.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::rDmdt
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
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::EntrainedGCPCombustionPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    phaseProperties(this),
    Tad
    (
      volScalarField
      (
        IOobject("Tadiabatic", mesh), mesh,
        dimensionedScalar("", dimTemperature, this->template get<scalar>("Tad"))
      )
    ),
    Hs1
    (
      volScalarField
      (
        IOobject("Hs1", mesh), mesh,
        dimensionedScalar("", dimVelocity*dimVelocity, 0)
      )
    ),
    Hs2
    (
      volScalarField
      (
        IOobject("Hs2", mesh), mesh,
        dimensionedScalar("", dimVelocity*dimVelocity, 0)
      )
    ),
    rb_
    (
      volScalarField
      (
        IOobject
        (
          "rb",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimVelocity,0)
      )
    ),
    rhoPropellant("rhoprop", dimDensity, this->template get<scalar>("propellantRho")),
    MR_("Mass Retained", dimless, this->template get<scalar>("Retained")),
    alphaOld
    (
      volScalarField
      (
        IOobject
        (
          "alphaOld",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless,0)
      )
    ),
    regressionAlpha
    (
      volScalarField
      (
        IOobject
        (
          "regressionAlpha",
          mesh.time().timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
        ),
        mesh
      )
    ),
    thinFlimModel_(this->template get<bool>("thinFlimModel")),
    Ug_
    (
      volVectorField
      (
        IOobject
        (
          "Ug",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector(0, 0, 0))
      )
    ),
    Up_
    (
      volVectorField
      (
        IOobject
        (
          "Up",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector(0, 0, 0))
      )
    ),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    ),
    eta(mesh, phaseProperties.subDict("GCP"))
{
    this->generatePairsAndSubModels
    (
        "interfaceTracking",
        interfaceTrackingModels_
    );
    regress_ = this->template get<bool>("regression");

    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[interfaceTrackingModelIter.key()];

        // Initially assume no mass transfer
        rDmdt_.set
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime)
            )
        );

        // Set Source terms
        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        Hs1 = phase1.thermo().he(phase1.thermo().p(), Tad);
        Hs2 = phase2.thermo().he(phase2.thermo().p(), Tad);
    }

    // clipping Regression Alpha
    regressionAlpha.clip(SMALL, 1-SMALL);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::
~EntrainedGCPCombustionPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    // return BasePhaseSystem::dmdt(key) + this->rDmdt(key);
    return BasePhaseSystem::dmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    // Fill the mass transfer rates with zero
    this->fillFields("dmdt", dimDensity/dimTime, dmdts);

    forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
        const volScalarField& rDmdt = *rDmdtIter();

        const Entrainment::factors mtf(eta.massTransfer());
        const scalar pcoeff = mtf.particles;
        const scalar gcoeff = (mtf.gas);

        this->addField(pair.phase1(), "dmdt", (pcoeff - MR_)*rDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", gcoeff*rDmdt, dmdts);

        // Subtract for Propellant Phase
        this->addField(this->phases()[2], "dmdt", -((pcoeff - MR_) + gcoeff)*rDmdt, dmdts);
    }
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    FatalErrorInFunction
        << "Mass Transfer() is not implemented"
        << exit(FatalError);
    return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
  autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
          BasePhaseSystem::heatTransfer();

  phaseSystem::heatTransferTable& eqns = eqnsPtr();

  forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
  {
    const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
    const volScalarField& rDmdt = *rDmdtIter();

    const Entrainment::factors mtf(eta.massTransfer());
    const scalar pcoeff = mtf.particles;
    const scalar gcoeff = (mtf.gas);

    const phaseModel& phase1 = pair.phase1();
    const phaseModel& phase2 = pair.phase2();

    // Enthalpy source (store this and don not calculate it again and again !!!)
    // const tmp<volScalarField> ths1(phase1.thermo().he(phase1.thermo().p(), Tad));
    // const volScalarField& hs1(ths1());
    //
    // const tmp<volScalarField> ths2(phase2.thermo().he(phase2.thermo().p(), Tad));
    // const volScalarField& hs2(ths2());

    // Equations
    fvScalarMatrix& eqn1 = *eqns[phase1.name()];
    fvScalarMatrix& eqn2 = *eqns[phase2.name()];

    // Implementation - 1
    eqn1 += - fvm::Sp((pcoeff - MR_)*rDmdt, eqn1.psi())
            + (pcoeff - MR_)*rDmdt*Hs1;
    eqn2 += - fvm::Sp(gcoeff*rDmdt, eqn2.psi())
            + gcoeff*rDmdt*Hs2;

    // Implementation - 2
    // Old Enthalpy
    // const tmp<volScalarField> th1(phase1.thermo().he());
    // const volScalarField& h1(th1());
    // const tmp<volScalarField> th2(phase2.thermo().he());
    // const volScalarField& h2(th2());
    //
    // eqn1 += pcoeff*rDmdt*(hs1 - h1);
    // eqn2 += gcoeff*rDmdt*(hs2 - h2);
  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
  autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
                BasePhaseSystem::momentumTransfer();
  if (thinFlimModel_)
  {
      phaseSystem::momentumTransferTable& eqns = eqnsPtr();
      forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
      {
          const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
          const volScalarField& rDmdt = *rDmdtIter();

          const Entrainment::factors mtf(eta.massTransfer());
          const scalar pcoeff = mtf.particles;
          const scalar gcoeff = (mtf.gas);

          const phaseModel& phase1 = pair.phase1();
          const phaseModel& phase2 = pair.phase2();

          // Equations
          fvVectorMatrix& eqn1 = *eqns[phase1.name()];
          fvVectorMatrix& eqn2 = *eqns[phase2.name()];

          // Momentum Source
          eqn1 += - fvm::Sp((pcoeff - MR_)*rDmdt, eqn1.psi())
                  + (pcoeff - MR_)*rDmdt*Up_;
          eqn2 += - fvm::Sp(gcoeff*rDmdt, eqn2.psi())
                  + gcoeff*rDmdt*Ug_;
      }
  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
void Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::solve()
{
  // Regress Propellant surface (Manipulate propellant volume fraction)
  if (regress_)
  {
      forAllIter
      (
          interfaceTrackingModelTable,
          interfaceTrackingModels_,
          interfaceTrackingModelIter
      )
      {
          const Entrainment::factors mtf(eta.massTransfer());

          word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
          volScalarField& alpha = this->db().template lookupObjectRef<volScalarField>(propellant);

          interfaceTrackingModelIter()->regress(mtf.particles, alpha);
      }
  }

  // Solve other phase volume fraction equations
  BasePhaseSystem::solve();
}

template<class BasePhaseSystem>
void Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        *rDmdt_[interfaceTrackingModelIter.key()]
                  = dimensionedScalar(dimDensity/dimTime);
    }

    //- return burning Rate, As and propellant density -> find mdot of propellant
    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        *rDmdt_[interfaceTrackingModelIter.key()]
        = interfaceTrackingModelIter()->dmdt()()*rhoPropellant;
    }

    // calculate velocity of the gas and particle source
    if (thinFlimModel_)
    {
        calculateVelocity();
    }
}

template<class BasePhaseSystem>
void Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::store()
{
  BasePhaseSystem::store();
  alphaOld = regressionAlpha;
  forAllIter
  (
      interfaceTrackingModelTable,
      interfaceTrackingModels_,
      interfaceTrackingModelIter
  )
  {
      interfaceTrackingModelIter()->store();
  }
}

template<class BasePhaseSystem>
void Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::calculateVelocity()
{
    //- Calculate velocity of the gas and particles comes
    //                into the combustion chamber

    // Velocity of gas and particle phase
    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
      const tmp<volScalarField> tAs = interfaceTrackingModelIter()->As();
      const tmp<volScalarField> trb = interfaceTrackingModelIter()->rb();
      const volScalarField& dmdt = *rDmdt_[interfaceTrackingModelIter.key()];

      const volScalarField& As(tAs());
      const volScalarField& rb(trb());

      eta.gasVelocity
      (
        (this->phases()[0].thermo().p()), dmdt, rb, As, Ug_
      );
    }
}

template<class BasePhaseSystem>
bool Foam::EntrainedGCPCombustionPhaseSystem<BasePhaseSystem>::read()
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
