/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Application
    reactingMultiphaseEulerFoam

Description
    Solver for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentum, heat and mass transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "simpleMatrix.H"
#include "timeProfiler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void displayMatrix(fvScalarMatrix TEqn, volScalarField T)
{
  // Usage:- displayMatrix(YiEqn, Y[i]);

  label NC = T.mesh().nCells(); //Number of cells
  simpleMatrix<scalar> A(NC); //Coeff.matrix
  // Initialization of matrix
  for(label i=0; i<NC; i++)
  {
    A.source()[i] = 0.0;
    for(label j=0; j<NC; j++)
    {
      A[i][j] = 0.0;
    }
  }
  // Assigning diagonal coefficients
  for(label i=0; i<NC; i++)
  {
    A[i][i] = TEqn.diag()[i];
    A.source()[i] += TEqn.source()[i];
  }
  // Assigning off-diagonal coefficients
  for(label faceI=0; faceI<TEqn.lduAddr().lowerAddr().size(); faceI++)
  {
    label l = TEqn.lduAddr().lowerAddr()[faceI];
    label u = TEqn.lduAddr().upperAddr()[faceI];
    A[l][u] = TEqn.upper()[faceI];
    A[u][l] = TEqn.upper()[faceI];
  }
  // Assigning contribution from BC
  forAll(T.boundaryField(), patchI)
  {
    const fvPatch &pp =
    T.boundaryField()[patchI].patch();
    forAll(pp, faceI)
    {
      label cellI = pp.faceCells()[faceI];
      A[cellI][cellI]
      += TEqn.internalCoeffs()[patchI][faceI];
      A.source()[cellI]
      += TEqn.boundaryCoeffs()[patchI][faceI];
    }
  }
  // Info << "====Coefficients of Matrix " << T.name() << " ====" << endl;
  for(label i=0; i<NC; i++)
  {
    for(label j=0; j<NC; j++)
    {
      Info<< A[i][j] << " ";
    }
    Info<< A.source()[i] << endl;
  }
  // Info<< "\n==> Solution: " << A.solve() << endl;
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Switch faceMomentum
    (
        pimple.dict().getOrDefault<Switch>("faceMomentum", false)
    );
    Switch partialElimination
    (
        pimple.dict().getOrDefault<Switch>("partialElimination", false)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    StopWatch totalTime;
    StopWatch runLoopTime;
    StopWatch courantNoTime;
    StopWatch pimpleLoopTime;
    StopWatch alphaEqnTime;
    StopWatch EEqnTime;
    StopWatch PEqnTime;
    StopWatch UEqnTime;
    StopWatch infoTime;
    EEqntimeProfiler EEqnProfile;

    totalTime.start();
    while (!runTime.end())
    {
      runLoopTime.start();

        courantNoTime.start();
        #include "readTimeControls.H"

        int nEnergyCorrectors
        (
            pimple.dict().getOrDefault<int>("nEnergyCorrectors", 1)
        );

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        courantNoTime.stop();

        pimpleLoopTime.start();
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
          alphaEqnTime.start();
            fluid.solve();  // Just regress propellant surface

            // #include "findPropellantCells.H"
            label purePropellantSize = 0;
            if (propellantIndex != -1)
            {
              const volScalarField& propellant = phases[propellantIndex];

              forAll(propellant, i)
              {
                if (propellant[i] == 1)
                {
                  purePropellantSize++;
                }
              }
            }
            labelList purePropellantCells(purePropellantSize);
            // scalar Tad = fluid.get<scalar>("Tad");
            // scalarField setTemp(purePropellantSize, Tad);
            scalarField setAlpha(purePropellantSize, SMALL);
            vectorField setVelocity(purePropellantSize, vector(0, 0, 0));
            if (propellantIndex != -1)
            {
              const volScalarField& propellant = phases[propellantIndex];

              label j = 0;
              forAll(propellant, i)
              {
                if (propellant[i] == 1)
                {
                  purePropellantCells[j] = i;
                  j++;
                }
              }
            }

            #include "alphaEqn.H"
            fluid.correct();
          alphaEqnTime.stop();

            // #include "YEqns.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                #include "EEqns.H"
                #include "pUf/pEqn.H"
            }
            else
            {
              UEqnTime.start();
                #include "pU/UEqns.H"
              UEqnTime.stop();

              EEqnTime.start();
                #include "EEqns.H"
              EEqnTime.stop();

              PEqnTime.start();
                #include "pU/pEqn.H"
              PEqnTime.stop();
            }

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }
        pimpleLoopTime.stop();

        infoTime.start();
        forAll(phases, phasei)
        {
            phaseModel& phase = phases[phasei];
            if (phase.index() == propellantIndex) continue;
            Info<< phase.name() << " min/max T "
                << min(phase.thermo().T()).value()
                << " - "
                << max(phase.thermo().T()).value()
                << endl;
        }
        runTime.write();

        runTime.printExecutionTime(Info);
        infoTime.stop();

      runLoopTime.stop();
    }
    totalTime.stop();

    Info << "Time Profiling: " << endl;
    double totalTimed = totalTime.getTotalTime();
    Info << "TotalTime: " << totalTimed << endl;
    Info << "      runTimeLoop : " << runLoopTime.getTotalTime()/totalTimed*100.0 << " % ( " << runLoopTime.getTotalTime() << " s)"<< endl;
    Info << "       ->  courantNoTime : " << courantNoTime.getTotalTime()/totalTimed*100.0 << " % ( " << courantNoTime.getTotalTime() << " s)"<< endl;
    Info << "       ->     info/write : " << infoTime.getTotalTime()/totalTimed*100.0 << " % ( " << infoTime.getTotalTime() << " s)"<< endl;
    Info << "       -> pimpleLoopTime : " << pimpleLoopTime.getTotalTime()/totalTimed*100.0 << " % ( " << pimpleLoopTime.getTotalTime() << " s)"<< endl;
    Info << "               -> alphaEqnTime : " << alphaEqnTime.getTotalTime()/totalTimed*100.0 << " % ( " << alphaEqnTime.getTotalTime() << " s)"<< endl;
    Info << "               ->         EEqn : " << EEqnTime.getTotalTime()/totalTimed*100.0 << " % ( " << EEqnTime.getTotalTime() << " s)"<< endl;
    Info << "               ->         UEqn : " << UEqnTime.getTotalTime()/totalTimed*100.0 << " % ( " << UEqnTime.getTotalTime() << " s)"<< endl;
    Info << "               ->         PEqn : " << PEqnTime.getTotalTime()/totalTimed*100.0 << " % ( " << PEqnTime.getTotalTime() << " s)"<< endl;

    EEqnProfile.displayTime();
    // fluid.finalize();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
