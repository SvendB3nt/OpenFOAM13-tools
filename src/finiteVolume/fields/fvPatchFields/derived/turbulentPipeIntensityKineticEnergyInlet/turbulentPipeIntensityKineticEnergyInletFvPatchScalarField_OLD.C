/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "turbulentPipeIntensityKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    calcHydDia_(!dict.found("hydraulicDiameter")),
    hydraulicDiameter_(dict.lookupOrDefault<scalar>("hydraulicDiameter", 0.0)),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    if (!calcHydDia_ && hydraulicDiameter_ < VSMALL)
    {
        FatalErrorInFunction
            << "Hydraulic Diameter should be larger than 0\n"
               " value given is " << hydraulicDiameter_ << nl
            << " on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    if (calcHydDia_)
    {
        updateHydDia();
    }

    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    hydraulicDiameter_(ptf.hydraulicDiameter_),
    UName_(ptf.UName_)
{}


Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    hydraulicDiameter_(ptf.hydraulicDiameter_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (patch().boundaryMesh().mesh().dynamic())
    {
        updateHydDia();
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvPatchScalarField& nup = 
        patch().lookupPatchField<volScalarField, scalar>("nu");
        
    const scalarField& Sf = patch().magSf();

    const scalar Re = gSum(mag(Up)*Sf/nup)*hydraulicDiameter_ / max(VSMALL,gSum(Sf));

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = 1.5 * sqr(0.16 * 1/pow4(sqr(Re)))*magSqr(Up);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    if (calcHydDia_)
    {
    writeEntry(os, "hydraulicDiameter", hydraulicDiameter_);
    }
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "value", *this);
}

void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::updateHydDia()
{
    hydraulicDiameter_ = sqrt(gSum(patch().magSf())*4/constant::mathematical::pi);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
    );
}

// ************************************************************************* //
