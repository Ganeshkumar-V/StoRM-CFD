if (adjustTimeStep)
{
    scalar factor = 1.01;
    if (propellantIndex != -1)
    {
        // Propellant burning surface transition - Stabilization
        const scalarField iF(phases[propellantIndex].internalField());

        // Algorithm for finding burning surface transition
        scalar transition = 0.0;
        forAll(iF, i)
        {
            if ((iF[i] < 1 - SMALL) && (iF[i] > 0.98)) transition = 1.0;
            else if ((iF[i] > SMALL) && (iF[i] < 0.05)) transition = 1.0;
        }
        reduce(transition, sumOp<scalar>());

        if (transition == 0)
        {
            maxCo = pimple.dict().get<scalar>("maxCo");
        }
        else
        {
            maxCo = pimple.dict().get<scalar>("transCo");
        }
    }
    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), factor);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
