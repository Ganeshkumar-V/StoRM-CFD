volVectorField Ugas
(
    IOobject
    (
      "U.gas",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volScalarField Tgas
(
    IOobject
    (
      "T.gas",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volVectorField Uparticles
(
    IOobject
    (
      "U.particles",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);
