/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::vorticityStreamFunction

Group
    grpFieldFunctionObjects

Description
    Computes the vorticityStreamFunction, the curl of the velocity.

    \f[
        \vec \omega = \nabla \times \vec U
    \f]

    where
    \vartable
        \vec \omega   | vorticityStreamFunction vector    [1/s]
        \vec U        | Velocity vector     [m/s]
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volVectorField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    vorticityStreamFunction1
    {
        // Mandatory entries (unmodifiable)
        type        vorticityStreamFunction;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: vorticityStreamFunction               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func vorticityStreamFunction
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::vorticityStreamFunction

SourceFiles
    vorticityStreamFunction.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_vorticityStreamFunction_H
#define functionObjects_vorticityStreamFunction_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class vorticityStreamFunction Declaration
\*---------------------------------------------------------------------------*/

class vorticityStreamFunction
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the vorticityStreamFunction field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("vorticityStreamFunction");


    // Constructors

        //- Construct from Time and dictionary
        vorticityStreamFunction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        vorticityStreamFunction(const vorticityStreamFunction&) = delete;

        //- No copy assignment
        void operator=(const vorticityStreamFunction&) = delete;


    //- Destructor
    virtual ~vorticityStreamFunction() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
