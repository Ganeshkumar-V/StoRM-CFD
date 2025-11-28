#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<multiPhaseSystem> fluidPtr
(
    multiPhaseSystem::New(mesh)
);
multiPhaseSystem& fluid = fluidPtr();
multiPhaseSystem::phaseModelList& phases = fluid.phases();

// Phases
phaseModel& gasPhase(phases[0]);  // Gas phase model
phaseModel& particlePhase(phases[1]);  // Particle phase model
// phaseModel& propellantPhase(phases[2]);  // Propellant phase model

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = phases[0].thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());

// calculate Mach Number
volScalarField Mach
(
    IOobject
    (
      "Mach",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mag(phases[0].U())/sqrt(phases[0].thermo().gamma()*p/phases[0].thermo().rho())
);

volScalarField rhoPropellant
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phases[0].thermo().rho()
);
if (phases.size() > 2)
{
    rhoPropellant = phases[2].thermo().rho();
}
