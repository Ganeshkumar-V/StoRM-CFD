{
  PtrList<volScalarField> dmdts(fluid.dmdts());

  volScalarField dalpha(0*fluid.phases()[0]);

  forAll(fluid.phases(), phasei)
  {
    phaseModel& phase = fluid.phases()[phasei];
    volScalarField& alpha = phase;

    if (phase.index() == gasIndex) continue;
    // Alpha Equation solved only for dispersed Phase
    else if (phase.index() == particleIndex)
    {
      const volScalarField& rho = phase.rho();
      const surfaceScalarField rhoPhi(fvc::interpolate(rho)*phase.phi());

      fvScalarMatrix alphaEqn
      (
        fvm::ddt(rho, alpha)
        + fvm::div(rhoPhi, alpha)
        == dmdts[phasei]
      );
      if (propellantIndex != -1)
      {
        alphaEqn.setValues(purePropellantCells, setAlpha);
      }
      alphaEqn.solve();
    }
    else {}

    // Summation of dispersed phase to calculate gas alpha
    dalpha += alpha;
  }

  // Volume fraction of continuous phase
  forAll(fluid.phases(), phasei)
  {
    phaseModel& phase = fluid.phases()[phasei];
    if(phase.stationary()) continue;
    if (!phase.moving()) continue;

    volScalarField& alpha = phase;

    if (phase.name() == gasPhase)
    {
      alpha = 1.0 - dalpha;
    }
    phase.alphaPhiRef() = fvc::interpolate(alpha)*phase.phi();
    phase.alphaRhoPhiRef() = fvc::interpolate(phase.rho())*phase.alphaPhi();
    phase.clip(SMALL, 1 - SMALL);
  }

}
