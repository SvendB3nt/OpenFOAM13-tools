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

#include "inletOutletPipeFvPatchField.H"

template<class Type>
void Foam::inletOutletPipeFvPatchField<Type>::selectHydDiaFnc()
{
    if (shape_ == "circle")
    {
        hydDiaFnc_ = &inletOutletPipeFvPatchField<Type>::hydDiaCircle;
    }
    else if (shape_ == "square")
    {
        hydDiaFnc_ = &inletOutletPipeFvPatchField<Type>::hydDiaSquare;
    }
    else if (shape_ == "rectangle")
    {
        hydDiaFnc_ = &inletOutletPipeFvPatchField<Type>::hydDiaRectangle;
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
template<class Type>
Foam::scalar Foam::inletOutletPipeFvPatchField<Type>::hydDiaCircle() const
{
    // Dh = sqrt(4*A/pi)
    const scalar A = gSum(this->patch().magSf());
    return sqrt(4.0 * A / constant::mathematical::pi);
}


template<class Type>
Foam::scalar Foam::inletOutletPipeFvPatchField<Type>::hydDiaSquare() const
{
    // Square of side s: A = s^2, Dh = 4A/P = s
    const scalar A = gSum(this->patch().magSf());
    return sqrt(A);
}


template<class Type>
Foam::scalar Foam::inletOutletPipeFvPatchField<Type>::hydDiaRectangle() const
{
    // Rectangle with sides a (short) and b (long), shapeRatio = a/b (0<r<=1)
    // A = a*b
    // a = sqrt(A * r)
    // b = sqrt(A / r)
    // Dh = 2ab/(a+b)
    const scalar A = gSum(this->patch().magSf());
    
    const scalar a = sqrt(A * shapeRatio_);
    const scalar b = sqrt(A / shapeRatio_);
    return 2.0 * a * b / (a + b);
}


template<class Type>
void Foam::inletOutletPipeFvPatchField<Type>::updateHydDia()
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

template<class Type>
void Foam::inletOutletPipeFvPatchField<Type>::checkInput(const dictionary& dict)
{
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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::inletOutletPipeFvPatchField<Type>::inletOutletPipeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(p, iF),
    calcHydDia_(true),
    hydraulicDiameter_(0.0),
    shape_("circle"),
    shapeRatio_(1.0),
    hydDiaFnc_(nullptr)
{
    this->phiName_ = "phi";
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0;
    fvPatchField<Type>::operator=(Zero);
}


template<class Type>
Foam::inletOutletPipeFvPatchField<Type>::inletOutletPipeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<Type>(p, iF, dict),
    calcHydDia_(!dict.found("hydraulicDiameter")),
    hydraulicDiameter_(dict.lookupOrDefault<scalar>("hydraulicDiameter", 0.0)),
    shape_(dict.lookupOrDefault<word>("shape", "circle")),
    shapeRatio_(1.0),
    hydDiaFnc_(nullptr)
{
 
    checkInput(dict); 

    selectHydDiaFnc();

    if (calcHydDia_)
    {
        updateHydDia();
    }
}


template<class Type>
Foam::inletOutletPipeFvPatchField<Type>::inletOutletPipeFvPatchField
(
    const inletOutletPipeFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    inletOutletFvPatchField<Type>(ptf, p, iF, mapper),
    calcHydDia_(ptf.calcHydDia_),
    hydraulicDiameter_(ptf.hydraulicDiameter_),
    shape_(ptf.shape_),
    shapeRatio_(ptf.shapeRatio_),
    hydDiaFnc_(nullptr)
{
    // set function pointer according to copied shape
    selectHydDiaFnc();
}


template<class Type>
Foam::inletOutletPipeFvPatchField<Type>::inletOutletPipeFvPatchField
(
    const inletOutletPipeFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(ptf, iF),
    calcHydDia_(ptf.calcHydDia_),
    hydraulicDiameter_(ptf.hydraulicDiameter_),
    shape_(ptf.shape_),
    shapeRatio_(ptf.shapeRatio_),
    hydDiaFnc_(nullptr)
{
    selectHydDiaFnc();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::inletOutletPipeFvPatchField<Type>::updateCoeffs()
{
    inletOutletFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::inletOutletPipeFvPatchField<Type>::write(Ostream& os) const
{
    if (!calcHydDia_)
    {
        writeEntry(os, "hydraulicDiameter", hydraulicDiameter_);
    }

    if (shape_ == "rectangle")
    {
        writeEntry(os, "shapeRatio", shapeRatio_);
    }
    
    writeEntryIfDifferent<word>(os, "shape", "circle", shape_);
    inletOutletFvPatchField<Type>::write(os);
}


// ************************************************************************* //
