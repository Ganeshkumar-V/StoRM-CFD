/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "waveTransmissiveVelocityFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::waveTransmissiveVelocityFvPatchField<Type>::waveTransmissiveVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("thermo:psi"),
    gamma_(0.0),
    inletValue_(Zero)
{}


template<class Type>
Foam::waveTransmissiveVelocityFvPatchField<Type>::waveTransmissiveVelocityFvPatchField
(
    const waveTransmissiveVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    inletValue_(ptf.inletValue_)
{}


template<class Type>
Foam::waveTransmissiveVelocityFvPatchField<Type>::waveTransmissiveVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    gamma_(dict.get<scalar>("gamma")),
    inletValue_("inletValue", dict, p.size())
{}


template<class Type>
Foam::waveTransmissiveVelocityFvPatchField<Type>::waveTransmissiveVelocityFvPatchField
(
    const waveTransmissiveVelocityFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_),
    inletValue_(ptpsf.inletValue_)
{}


template<class Type>
Foam::waveTransmissiveVelocityFvPatchField<Type>::waveTransmissiveVelocityFvPatchField
(
    const waveTransmissiveVelocityFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_),
    inletValue_(ptpsf.inletValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::waveTransmissiveVelocityFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip =
        this->patch().template
            lookupPatchField<volScalarField, scalar>(psiName_);

    const surfaceScalarField& phi =
        this->db().template lookupObject<surfaceScalarField>(this->phiName_);

    fvsPatchField<scalar> phip =
        this->patch().template
            lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop =
            this->patch().template
                lookupPatchField<volScalarField, scalar>(this->rhoName_);

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf() + sqrt(gamma_/psip);
}

template<class Type>
void Foam::waveTransmissiveVelocityFvPatchField<Type>::updateCoeffs()
{
  if (this->updated())
  {
      return;
  }

  const fvMesh& mesh = this->internalField().mesh();

  word ddtScheme
  (
      mesh.ddtScheme(this->internalField().name())
  );
  scalar deltaT = this->db().time().deltaTValue();

  const GeometricField<Type, fvPatchField, volMesh>& field =
      this->db().objectRegistry::template
      lookupObject<GeometricField<Type, fvPatchField, volMesh>>
      (
          this->internalField().name()
      );

  // Calculate the advection speed of the field wave
  // If the wave is incoming set the speed to 0.
  const scalarField w(Foam::max(advectionSpeed(), scalar(0)));

  // Calculate the field wave coefficient alpha (See notes)
  const scalarField alpha(w*deltaT*this->patch().deltaCoeffs());

  label patchi = this->patch().index();

  const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            this->phiName_
        );

  // Non-reflecting outflow boundary
  // If this->lInf() defined setup relaxation to the value this->fieldInf().
  if (this->lInf() > 0)
  {
      // Calculate the field relaxation coefficient k (See notes)
      const scalarField k(w*deltaT/this->lInf());

      if
      (
          ddtScheme == fv::EulerDdtScheme<scalar>::typeName
       || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
      )
      {
          this->refValue() = pos0(phip)*
          (
              field.oldTime().boundaryField()[patchi] + k*this->fieldInf()
          )/(1.0 + k) + (1.0 - pos0(phip))*inletValue_;

          this->valueFraction() = pos0(phip)*(1.0 + k)/(1.0 + alpha + k) + (1.0 - pos0(phip));
      }
      else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
      {
          this->refValue() = pos0(phip)*
          (
              2.0*field.oldTime().boundaryField()[patchi]
            - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
            + k*this->fieldInf()
          )/(1.5 + k) + (1.0 - pos0(phip))*inletValue_;

          this->valueFraction() = pos0(phip)*(1.5 + k)/(1.5 + alpha + k) + (1.0 - pos0(phip));
      }
      else if
      (
          ddtScheme == fv::localEulerDdtScheme<scalar>::typeName
      )
      {
          const volScalarField& rDeltaT =
              fv::localEulerDdt::localRDeltaT(mesh);

          // Calculate the field wave coefficient alpha (See notes)
          const scalarField alpha
          (
              w*this->patch().deltaCoeffs()/rDeltaT.boundaryField()[patchi]
          );

          // Calculate the field relaxation coefficient k (See notes)
          const scalarField k(w/(rDeltaT.boundaryField()[patchi]*this->lInf()));

          this->refValue() = pos0(phip)*
          (
              field.oldTime().boundaryField()[patchi] + k*this->fieldInf()
          )/(1.0 + k) + (1.0 - pos0(phip))*inletValue_;

          this->valueFraction() = pos0(phip)*(1.0 + k)/(1.0 + alpha + k) + (1.0 - pos0(phip));
      }
      else
      {
          FatalErrorInFunction
              << ddtScheme << nl
              << "    on patch " << this->patch().name()
              << " of field " << this->internalField().name()
              << " in file " << this->internalField().objectPath()
              << exit(FatalError);
      }
  }
  else
  {
      if
      (
          ddtScheme == fv::EulerDdtScheme<scalar>::typeName
       || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
      )
      {
          this->refValue() = field.oldTime().boundaryField()[patchi];

          this->valueFraction() = 1.0/(1.0 + alpha);
      }
      else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
      {
          this->refValue() =
          (
              2.0*field.oldTime().boundaryField()[patchi]
            - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
          )/1.5;

          this->valueFraction() = 1.5/(1.5 + alpha);
      }
      else if
      (
          ddtScheme == fv::localEulerDdtScheme<scalar>::typeName
      )
      {
          const volScalarField& rDeltaT =
              fv::localEulerDdt::localRDeltaT(mesh);

          // Calculate the field wave coefficient alpha (See notes)
          const scalarField alpha
          (
              w*this->patch().deltaCoeffs()/rDeltaT.boundaryField()[patchi]
          );

          this->refValue() = field.oldTime().boundaryField()[patchi];

          this->valueFraction() = 1.0/(1.0 + alpha);
      }
      else
      {
          FatalErrorInFunction
              << ddtScheme
              << "\n    on patch " << this->patch().name()
              << " of field " << this->internalField().name()
              << " in file " << this->internalField().objectPath()
              << exit(FatalError);
      }
  }

  mixedFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void Foam::waveTransmissiveVelocityFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", this->rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);

    os.writeEntry("gamma", gamma_);
    inletValue_.writeEntry("inletValue", os);
    // os.writeEntry("inletValue", inletValue_);

    if (this->lInf() > SMALL)
    {
        os.writeEntry("fieldInf", this->fieldInf());
        os.writeEntry("lInf", this->lInf());
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
