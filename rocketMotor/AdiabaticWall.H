inline void ImposeWall(volScalarField& psi, const volScalarField& alpha)
{
  const labelList& Own(psi.mesh().owner());
  const labelList& Nei(psi.mesh().neighbour());

  scalar One(0.999); //1.0 - SMALL);

  forAll(Own, celli)
  {
      if ((alpha[Own[celli]] < One) && (alpha[Nei[celli]] >= One))
      {
          psi[Nei[celli]] = psi[Own[celli]];
      }
      else if ((alpha[Nei[celli]] < One) && (alpha[Own[celli]] >= One))
      {
          psi[Own[celli]] = psi[Nei[celli]];
      }
  }


  forAll(psi.mesh().boundary(),patchi)
  {
    if (isType<processorFvPatch>(psi.mesh().boundary()[patchi]))
    {
      const processorPolyPatch& pp
          = refCast<const processorPolyPatch>(psi.mesh().boundaryMesh()[patchi]);
      if (pp.neighbour())
      {
        const scalarField& psibF(psi.boundaryField()[patchi].patchNeighbourField());
        const scalarField& alphaF(alpha.boundaryField()[patchi].patchInternalField());
        const scalarField& alphabF(alpha.boundaryField()[patchi].patchNeighbourField());
        const labelList& fC(psi.mesh().boundary()[patchi].faceCells());
        forAll(fC, i)
        {
          if ((alphabF[i] < One) && (alphaF[i] >= One))
          {
            psi[fC[i]] = psibF[i];
          }
        }
      }
    }
  }

}

inline void correctFlux(surfaceScalarField& flux, const volScalarField& propellant)
{
  /***/
    const fvMesh& mesh(propellant.mesh());
    const labelList& Own(mesh.owner());
    const labelList& Nei(mesh.neighbour());
    scalar One(1 - SMALL);

    forAll(Own, i)
    {
      if(propellant[Own[i]] == One || propellant[Nei[i]] == One)
      {
         flux[i] = 0;
      }
      else continue;
    }
    /***/
}

inline void correctU(volVectorField& U, const volScalarField& propellant)
{
  scalar One(1 - SMALL);

  forAll(propellant, i)
  {
    if(propellant[i] == One)
    {
       U[i] = vector(0, 0, 0);
    }
    else continue;
  }
  // U = neg(propellant - (1 - SMALL))*U;
}
