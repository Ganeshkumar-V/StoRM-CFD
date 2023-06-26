/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "transonicFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transonicFvPatchScalarField
::transonicFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    backPressure_(101325),
    gas_("gas"),
    gamma_(0.0),
    R_(1.0),
    supersonic(false)
{
    this->refValue() = 101325;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::transonicFvPatchScalarField
::transonicFvPatchScalarField
(
    const transonicFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    backPressure_(ptf.backPressure_),
    gas_(ptf.gas_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    supersonic(ptf.supersonic)
{}


Foam::transonicFvPatchScalarField
::transonicFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    backPressure_(dict.get<scalar>("backPressure")),
    gas_(dict.get<word>("gas")),
    gamma_(dict.get<scalar>("gamma")),
    R_(dict.get<scalar>("R")),
    supersonic(dict.getOrDefault<bool>("supersonic", false))
{
    this->refValue() = backPressure_;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->patchInternalField());
    }

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::transonicFvPatchScalarField
    ::transonicFvPatchScalarField
(
    const transonicFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    backPressure_(ptf.backPressure_),
    gas_(ptf.gas_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    supersonic(ptf.supersonic)
{}


Foam::transonicFvPatchScalarField
    ::transonicFvPatchScalarField
(
    const transonicFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    backPressure_(ptf.backPressure_),
    gas_(ptf.gas_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    supersonic(ptf.supersonic)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::transonicFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Gas phase density
    const fvPatchScalarField& Tg
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("T", gas_)
        )
    );

    // Gas phase Velocity
    const fvPatchVectorField& Ug
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            IOobject::groupName("U", gas_)
        )
    );

    if (supersonic)
    {
        if (this->size() > 0)   // To  stop other processors from calculation
        {
            // Check for normal shock possibility at the exit
            const scalarField M(mag(Ug)/sqrt(gamma_*R_*Tg));
            const scalarField P(*this);

            scalar Sp(sum(this->patch().magSf()));
            scalar Mavg = sum(M*this->patch().magSf())/Sp;
            scalar Pavg = sum(P*this->patch().magSf())/Sp;

            scalar Pp = (2*gamma_*sqr(Mavg) - (gamma_ - 1.0))*Pavg/(gamma_ + 1);
            if (Pp > backPressure_)
            {
                // Supersonic Flow - Normal Shock didn't form
                this->valueFraction() = 0.0;
            }
            else
            {
                // Subsonic Flow - Normal Shock at the exit
                this->valueFraction() = 1.0;
            }
        }
    }
    else
    {
        if (this->size() > 0) // To  stop other processors from calculation
        {
            // Mach Number internalField
            const scalarField Mc(mag(Ug.patchInternalField())/sqrt(gamma_*R_*Tg.patchInternalField()));

            if (max(Mc) > 1.0)
            {
                this->valueFraction() = 0.0;
            }
            else
            {
                this->valueFraction() = 1.0;
            }
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::transonicFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry<scalar>("backPressure", backPressure_);
    os.writeEntry<word>("gas", gas_);
    os.writeEntry<bool>("supersonic", supersonic);
    os.writeEntry<scalar>("gamma", gamma_);
    os.writeEntry<scalar>("R", R_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        transonicFvPatchScalarField
    );
}

// ************************************************************************* //
