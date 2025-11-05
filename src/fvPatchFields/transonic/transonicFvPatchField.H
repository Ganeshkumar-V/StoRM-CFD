/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::transonicFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a phase fraction condition based on the
    local flow conditions, whereby the values are constrained to lay between
    user-specified upper and lower bounds.  The behaviour is described by:

    if alpha > upperBound:
    - apply a fixed value condition, with a uniform level of the upper bound

    if lower bound <= alpha <= upper bound:
    - apply a  zero-gradient condition

    if alpha < lowerBound:
    - apply a fixed value condition, with a uniform level of the lower bound

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        lowerBound   | lower bound for clipping | yes        |
        upperBound   | upper bound for clipping | yes        |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            transonic;
        lowerBound      0.0;
        upperBound      0.9;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    transonicFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef transonicFvPatchScalarField_H
#define transonicFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class transonicFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class transonicFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        //- Back Pressure
        scalar backPressure_;

        //- Gas name
        word gas_;

        //- gas phase properties
        scalar gamma_;
        scalar R_;

        //- Supersonic to subsonic switch
        bool supersonic;

public:

    //- Runtime scalar information
    TypeName("transonic");


    // Constructors

        //- Construct from patch and internal field
        transonicFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transonicFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  transonicFvPatchScalarField onto a new patch
        transonicFvPatchScalarField
        (
            const transonicFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        transonicFvPatchScalarField
        (
            const transonicFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new transonicFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        transonicFvPatchScalarField
        (
            const transonicFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new transonicFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Find Rankine Hugonoit Normal Shock downstream Pressure


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
