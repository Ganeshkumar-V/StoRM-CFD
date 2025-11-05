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

#include "waveTransmissiveShockFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::waveTransmissiveShockFvPatchField<Type>::waveTransmissiveShockFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("thermo:psi"),
    gas_("gas"),
    gamma_(0.0),
    R_(1.0),
    supersonic(false)
{}


template<class Type>
Foam::waveTransmissiveShockFvPatchField<Type>::waveTransmissiveShockFvPatchField
(
    const waveTransmissiveShockFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gas_(ptf.gas_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    supersonic(ptf.supersonic)
{}


template<class Type>
Foam::waveTransmissiveShockFvPatchField<Type>::waveTransmissiveShockFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    gas_(dict.get<word>("gas")),
    gamma_(dict.get<scalar>("gamma")),
    R_(dict.get<scalar>("R")),
    supersonic(dict.getOrDefault<bool>("supersonic", false)),
    transition(dict.getOrDefault<bool>("transition", false))
{}


template<class Type>
Foam::waveTransmissiveShockFvPatchField<Type>::waveTransmissiveShockFvPatchField
(
    const waveTransmissiveShockFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    gas_(ptpsf.gas_),
    gamma_(ptpsf.gamma_),
    R_(ptpsf.R_),
    supersonic(ptpsf.supersonic)
{}


template<class Type>
Foam::waveTransmissiveShockFvPatchField<Type>::waveTransmissiveShockFvPatchField
(
    const waveTransmissiveShockFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    gas_(ptpsf.gas_),
    gamma_(ptpsf.gamma_),
    R_(ptpsf.R_),
    supersonic(ptpsf.supersonic),
    transition(ptpsf.transition)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::waveTransmissiveShockFvPatchField<Type>::advectionSpeed() const
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
void Foam::waveTransmissiveShockFvPatchField<Type>::updateCoeffs()
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
          this->refValue() =
          (
              field.oldTime().boundaryField()[patchi] + k*this->fieldInf()
          )/(1.0 + k);

          this->valueFraction() = (1.0 + k)/(1.0 + alpha + k);
      }
      else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
      {
          this->refValue() =
          (
              2.0*field.oldTime().boundaryField()[patchi]
            - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
            + k*this->fieldInf()
          )/(1.5 + k);

          this->valueFraction() = (1.5 + k)/(1.5 + alpha + k);
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

          this->refValue() =
          (
              field.oldTime().boundaryField()[patchi] + k*this->fieldInf()
          )/(1.0 + k);

          this->valueFraction() = (1.0 + k)/(1.0 + alpha + k);
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
     FatalErrorInFunction <<
            "lInf cannot be less than or equal to zero!"
            << exit(FatalError);
  }

  // Check for the need of subsonic transition
  // Gas phase density
  const fvPatchScalarField& Tg
  (
      this->patch().template lookupPatchField<volScalarField, scalar>
      (
          IOobject::groupName("T", gas_)
      )
  );

  // Gas phase Velocity
  const fvPatchVectorField& Ug
  (
      this->patch().template lookupPatchField<volVectorField, vector>
      (
          IOobject::groupName("U", gas_)
      )
  );

  if ((supersonic && (this->size() > 0)) && !transition)
  {
      // Check for normal shock possibility at the exit
      const scalarField M(mag(Ug)/sqrt(gamma_*R_*Tg));
      const scalarField P(mag(*this));

      scalar Sp(sum(this->patch().magSf()));
      scalar Mavg = sum(M*this->patch().magSf())/Sp;
      scalar Pavg = sum(P*this->patch().magSf())/Sp;

      scalar Pp = (2*gamma_*sqr(min(M)) - (gamma_ - 1.0))*min(P)/(gamma_ + 1);
      if (Pp <= 101325)
      {
	  transition = true;
	  Info << "Normal Shock at exit! -> M1 = " << Mavg << " P1 = " << Pavg << endl;
          // Subsonic Flow - Normal Shock at the exit
          this->valueFraction() = 1.0;
          this->refValue() = this->fieldInf();
      }
  }
  
  if (transition)
  {
          this->valueFraction() = 1.0;
          this->refValue() = this->fieldInf();
  }

  mixedFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void Foam::waveTransmissiveShockFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", this->rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);

    os.writeEntry<word>("gas", gas_);
    os.writeEntry<bool>("supersonic", supersonic);
    os.writeEntry<bool>("transition", transition);
    os.writeEntry<scalar>("gamma", gamma_);
    os.writeEntry<scalar>("R", R_);

    if (this->lInf() > SMALL)
    {
        os.writeEntry("fieldInf", this->fieldInf());
        os.writeEntry("lInf", this->lInf());
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
