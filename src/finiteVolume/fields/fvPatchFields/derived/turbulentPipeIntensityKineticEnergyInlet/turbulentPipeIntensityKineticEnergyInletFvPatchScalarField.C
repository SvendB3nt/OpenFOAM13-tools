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
    inletOutletPipeFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    this->calcHydDia_ = !dict.found("hydraulicDiameter");
    this->shape_ = dict.lookupOrDefault<word>("shape", "circle");
    this->shapeRatio_ = 1.0;
    this->hydDiaFnc_ = nullptr;
    
    checkInput(dict);

    selectHydDiaFnc();

    if (this->calcHydDia_)
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
    inletOutletPipeFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    nuName_(ptf.nuName_)
{

    this->selectHydDiaFnc();
}


Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletPipeFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    nuName_(ptf.nuName_)
{
    this->selectHydDiaFnc();
}

Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletPipeFvPatchScalarField(p, iF),
    UName_("U"),
    nuName_("nu")   
{
    this->selectHydDiaFnc();
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (patch().boundaryMesh().mesh().dynamic())
    {
        this->updateHydDia();
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvPatchScalarField& nup = 
        patch().lookupPatchField<volScalarField, scalar>(nuName_);
        
    const scalarField& Sf = patch().magSf();

    const scalar Re = gSum(mag(Up)*Sf/nup)*this->hydraulicDiameter_ / max(VSMALL,gSum(Sf));

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = 1.5 * sqr(0.16 * 1/sqrt(sqrt(sqrt(Re))))*magSqr(Up);
    this->valueFraction() = neg(phip);

    inletOutletPipeFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    if (!this->calcHydDia_)
    {
        writeEntry(os, "hydraulicDiameter", this->hydraulicDiameter_);
    }

    if (this->shape_ == "rectangle")
    {
        writeEntry(os, "shapeRatio", this->shapeRatio_);
    }
    
    writeEntryIfDifferent<word>(os, "shape", "circle", this->shape_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "value", *this);
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
