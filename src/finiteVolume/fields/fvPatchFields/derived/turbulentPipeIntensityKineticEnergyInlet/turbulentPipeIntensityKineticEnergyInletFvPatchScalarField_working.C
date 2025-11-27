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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::selectHydDiaFnc()
{
    if (shape_ == "circle")
    {
        hydDiaFnc_ = &turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::hydDiaCircle;
    }
    else if (shape_ == "square")
    {
        hydDiaFnc_ = &turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::hydDiaSquare;
    }
    else if (shape_ == "rectangle")
    {
        hydDiaFnc_ = &turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::hydDiaRectangle;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown hydraulic-diameter shape: " << shape_
            << ". Valid shapes: circle, square, rectangle."
            << exit(FatalError);
    }
}


// Shape-specific hydraulic diameter formulas

Foam::scalar Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::hydDiaCircle() const
{
    // Dh = sqrt(4*A/pi)
    const scalar A = gSum(patch().magSf());
    return sqrt(4.0 * A / constant::mathematical::pi);
}

Foam::scalar Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::hydDiaSquare() const
{
    // Square of side s: A = s^2, Dh = 4A/P = s
    const scalar A = gSum(patch().magSf());
    return sqrt(A);
}

Foam::scalar Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::hydDiaRectangle() const
{
    // Rectangle with sides a (short) and b (long), shapeRatio = a/b (0<r<=1)
    // A = a*b
    // a = sqrt(A * r)
    // b = sqrt(A / r)
    // Dh = 2ab/(a+b)
    const scalar A = gSum(patch().magSf());
    
    const scalar a = sqrt(A * shapeRatio_);
    const scalar b = sqrt(A / shapeRatio_);
    return 2.0 * a * b / (a + b);
}


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
    UName_(dict.lookupOrDefault<word>("U", "U")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    shape_(dict.lookupOrDefault<word>("shape", "circle")),
    shapeRatio_(1.0), 
    hydDiaFnc_(nullptr)
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

    if (shape_ == "rectangle")
    {
        if (!dict.found("shapeRatio"))
        {
            FatalErrorInFunction
                << "shapeRatio must be provided for shape rectangle on patch "
                << this->patch().name() << ".\n"
                << "Provide shapeRatio = a/b where a is the shortest side (0 < shapeRatio <= 1)."
                << exit(FatalError);
        }
        shapeRatio_ = dict.lookup<scalar>("shapeRatio");

        if (shapeRatio_ <= VSMALL || shapeRatio_ > 1.0)
        {
            FatalErrorInFunction
                << "Invalid shapeRatio: " << shapeRatio_ << " on patch "
                << this->patch().name() << ".\n"
                << "Must satisfy 0 < shapeRatio <= 1 (a/b where a is shortest side)."
                << exit(FatalError);
        }
    }

    selectHydDiaFnc();

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
    calcHydDia_(ptf.calcHydDia_),
    hydraulicDiameter_(ptf.hydraulicDiameter_),
    UName_(ptf.UName_),
    nuName_(ptf.nuName_),
    shape_(ptf.shape_),
    shapeRatio_(ptf.shapeRatio_),
    hydDiaFnc_(nullptr)
{
    // set function pointer according to copied shape
    selectHydDiaFnc();
}


Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    calcHydDia_(ptf.calcHydDia_),
    hydraulicDiameter_(ptf.hydraulicDiameter_),
    UName_(ptf.UName_),
    nuName_(ptf.nuName_),
    shape_(ptf.shape_),
    shapeRatio_(ptf.shapeRatio_),
    hydDiaFnc_(nullptr)
{
    selectHydDiaFnc();
}

Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::
turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    calcHydDia_(true),
    hydraulicDiameter_(0.0),
    UName_("U"),
    nuName_("nu"),
    shape_("circle"),
    shapeRatio_(1.0), 
    hydDiaFnc_(nullptr)    
{
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
        updateHydDia();
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvPatchScalarField& nup = 
        patch().lookupPatchField<volScalarField, scalar>(nuName_);
        
    const scalarField& Sf = patch().magSf();

    const scalar Re = gSum(mag(Up)*Sf/nup)*hydraulicDiameter_ / max(VSMALL,gSum(Sf));

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = 1.5 * sqr(0.16 * 1/sqrt(sqrt(sqrt(Re))))*magSqr(Up);
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

    if (shape_ == "rectangle")
    {
        writeEntry(os, "shapeRatio", shapeRatio_);
    }
    
    writeEntry(os, "shape", shape_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "value", *this);
}

void Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField::updateHydDia()
{
    // no conditionals here — hydDiaFnc_ chosen in constructor / copy ctors
    if (!hydDiaFnc_)
    {
        FatalErrorInFunction
            << "hydDiaFnc_ not set for patch " << this->patch().name()
            << exit(FatalError);
    }

    hydraulicDiameter_ = (this->*hydDiaFnc_)();
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
