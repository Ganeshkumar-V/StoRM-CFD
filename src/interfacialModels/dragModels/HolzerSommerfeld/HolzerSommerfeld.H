/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::particleDragModels::HolzerSommerfeld

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    HolzerSommerfeld.C

\*---------------------------------------------------------------------------*/

#ifndef HolzerSommerfeld_H
#define HolzerSommerfeld_H

#include "particleDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace particleDragModels
{

/*---------------------------------------------------------------------------*\
                           Class HolzerSommerfeld Declaration
\*---------------------------------------------------------------------------*/

class HolzerSommerfeld
:
    public particleDragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

    // typedefs
    typedef DimensionedField<scalar, fvMesh> Internal;
    typedef FieldField<fvPatchField, scalar> Boundary;


public:

    //- Runtime type information
    TypeName("HolzerSommerfeld");


    // Constructors

        //- Construct from a dictionary and a phase pair
        HolzerSommerfeld
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~HolzerSommerfeld();

    //
    struct fall
    {
      dimensionedScalar phiPara;
      dimensionedScalar phiPerp;
      dimensionedScalar A;
      dimensionedScalar B;
      dimensionedScalar C;
    };

    fall length_;
    fall cross_;

    dimensionedScalar phi_;
    dimensionedScalar d_;
    dimensionedScalar t_;
    dimensionedScalar dV_;

    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
