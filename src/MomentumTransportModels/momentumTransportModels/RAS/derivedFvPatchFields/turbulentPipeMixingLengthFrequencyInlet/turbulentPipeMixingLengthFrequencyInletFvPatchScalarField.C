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

#include "turbulentPipeMixingLengthFrequencyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPipeMixingLengthFrequencyInletFvPatchScalarField::
turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletPipeFvPatchScalarField(p, iF),
    turbulentLengthScale_(dict.lookupOrDefault<scalar>("L", 0.07)),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    this->calcHydDia_= !dict.found("hydraulicDiameter");
    this->hydraulicDiameter_ = dict.lookupOrDefault<scalar>("hydraulicDiameter", 0.0);
    this->shape_ = dict.lookupOrDefault<word>("shape", "circle");
    this->shapeRatio_ = 1.0; 
    this->hydDiaFnc_ = nullptr;

    this->checkInput(dict);

    this->selectHydDiaFnc();

    if (this->calcHydDia_)
    {
        this->updateHydDia();
    }

    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


turbulentPipeMixingLengthFrequencyInletFvPatchScalarField::
turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
(
    const turbulentPipeMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    inletOutletPipeFvPatchScalarField(ptf, p, iF, mapper),
    turbulentLengthScale_(ptf.turbulentLengthScale_),
    kName_(ptf.kName_),
    Cmu_(ptf.Cmu_)
{}


turbulentPipeMixingLengthFrequencyInletFvPatchScalarField::
turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
(
    const turbulentPipeMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletPipeFvPatchScalarField(ptf, iF),
    turbulentLengthScale_(ptf.turbulentLengthScale_),
    kName_(ptf.kName_),
    Cmu_(ptf.Cmu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentPipeMixingLengthFrequencyInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (patch().boundaryMesh().mesh().dynamic())
    {
        this->updateHydDia();
    }
    
    const scalar Cmu25 = pow(Cmu_, 0.25);

    const fvPatchScalarField& kp =
        patch().lookupPatchField<volScalarField, scalar>(kName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = 
        sqrt(kp)/(Cmu25 * turbulentLengthScale_ * this->hydraulicDiameter_);
    
    this->valueFraction() = neg(phip);

    inletOutletPipeFvPatchScalarField::updateCoeffs();
}


void turbulentPipeMixingLengthFrequencyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "L", turbulentLengthScale_);
    if (!this->calcHydDia_)
    {
        writeEntry(os, "hydraulicDiameter", this->hydraulicDiameter_);
    }

    if (this->shape_ == "rectangle")
    {
        writeEntry(os, "shapeRatio", this->shapeRatio_);
    }
    
    writeEntry(os, "Cmu", Cmu_);
    writeEntryIfDifferent<word>(os, "shape", "circle",this->shape_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
