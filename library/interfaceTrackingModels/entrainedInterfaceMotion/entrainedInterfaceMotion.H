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

Class
    Foam::interfaceTrackingModels::entrainedInterfaceMotion

Description

SourceFiles
    entrainedInterfaceMotion.C

\*---------------------------------------------------------------------------*/

#ifndef entrainedInterfaceMotion_H
#define entrainedInterfaceMotion_H

#include "interfaceTrackingModel.H"
#include "Surface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceTrackingModels
{

/*---------------------------------------------------------------------------*\
                          Class entrainedInterfaceMotion Declaration
\*---------------------------------------------------------------------------*/

class entrainedInterfaceMotion
:
    public interfaceTrackingModel
{
public:

    //- Runtime type information
    TypeName("entrainedInterfaceMotion");

    scalar n;
    scalar f;
    dimensionedScalar a;

    // Flame surface
    Surface flame_;

    // Propellant Bed surface
    Surface bed_;

    // mass transfer rate
    tmp<volScalarField> dmdt_;

    // Retainer
    scalar MR;

    // Constructors

        //- Construct from components
        entrainedInterfaceMotion
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~entrainedInterfaceMotion();

    // Member Functions
    virtual void correct();

    //- The burning rate
    virtual tmp<volScalarField> rb() const;

    //- Interface Area / Volume
    virtual tmp<volScalarField> As() const;

    //- Mass Transfer Rate
    virtual tmp<volScalarField> dmdt() const;

    // These virtual functions are not used in this tracking models
    virtual void regress(const scalar fp, volScalarField& alpha);
    virtual void findInterface(const volScalarField& alpha){}

    virtual void store();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTrackingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
