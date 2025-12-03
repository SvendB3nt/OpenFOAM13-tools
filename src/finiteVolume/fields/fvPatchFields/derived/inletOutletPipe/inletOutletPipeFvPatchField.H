/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Class
    Foam::inletOutletPipeFvPatchField

Description
<<<<<<< HEAD
    This boundary condition provides a generic outflow condition, with
    specified inflow for the case of return flow.
=======
    This boundary condition is an extention of the inletOutlet. 
    It provides child classes with access to a user defined hydraulic 
    diameter or one that is calculatedfor simple cross sectional 
    geometries.
>>>>>>> 6da4759 (Added decriptions to header files)

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | Flux field name         | no          | phi
        inletValue   | Inlet value for reverse flow | yes    |
<<<<<<< HEAD
=======
        hydraulicDiameter | hydraulic diameter for the patch | no | - 
        shape        | geometry for hydraulic diameter (circle|square|rectangle) | no | circle
        shapeRatio   | (for rectangle only) a/b where a<=b (0<ratio<=1) | required if shape rectangle | -
>>>>>>> 6da4759 (Added decriptions to header files)
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
<<<<<<< HEAD
        type            inletOutlet;
=======
        type            inletOutletPipe;
>>>>>>> 6da4759 (Added decriptions to header files)
        phi             phi;
        inletValue      uniform 0;
        value           uniform 0;
    }
<<<<<<< HEAD
=======

    <patchName>
    {
        type            inletOutletPipe;
        phi             phi;
        inletValue      uniform 0;
        value           uniform 0;
        hydraulicDiameter    1;
    }

    <patchName>
    {
        type            inletOutletPipe;
        phi             phi;
        inletValue      uniform 0;
        value           uniform 0;
        shape           square;
    }

    <patchName>
    {
        type            inletOutletPipe;
        phi             phi;
        inletValue      uniform 0;
        value           uniform 0;
        shape           rectangle;
        shapeRatio      0.6;
    }

>>>>>>> 6da4759 (Added decriptions to header files)
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch faces.

    Note:
        Sign conventions:
      - Positive flux (out of domain): apply zero-gradient condition
      - Negative flux (into of domain): apply the "inletValue" fixed-value

See also
    Foam::mixedFvPatchField
    Foam::zeroGradientFvPatchField
    Foam::outletInletFvPatchField

SourceFiles
    inletOutletPipeFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef inletOutletPipeFvPatchField_H
#define inletOutletPipeFvPatchField_H

#include "inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class inletOutletPipeFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class inletOutletPipeFvPatchField
:
    public inletOutletFvPatchField<Type>
{

protected:

    // Protected data

    //- Switch for calculating the hydraulic diameter
    bool calcHydDia_;

    //- Hydraulic diameter of patch
    scalar hydraulicDiameter_;

    //- Shape specifier: "circle", "square", "rectangle", ...
    word shape_;

    //- For rectangles: shapeRatio = a/b (a = shortest side, b = longest side)
    scalar shapeRatio_;

    //- Function pointer to the hydraulic diameter calculation selected at construction
    scalar (inletOutletPipeFvPatchField::*hydDiaFnc_)() const;

    //- Helper to set hydDiaFnc_ based on shape_
    void selectHydDiaFnc();

    //- Shape-specific hydraulic diameter calculators
    scalar hydDiaCircle() const;
    scalar hydDiaSquare() const;
    scalar hydDiaRectangle() const;

    //- Update the value of the hydraulic diameter
    void updateHydDia();

    //- check user input from dict
    void checkInput(const dictionary&);


public:

    //- Runtime type information
    TypeName("inletOutletPipe");


    // Constructors

        //- Construct from patch and internal field
        inletOutletPipeFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inletOutletPipeFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inletOutletPipeFvPatchField onto a new patch
        inletOutletPipeFvPatchField
        (
            const inletOutletPipeFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        inletOutletPipeFvPatchField(const inletOutletPipeFvPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        inletOutletPipeFvPatchField
        (
            const inletOutletPipeFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new inletOutletPipeFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Attributes


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "inletOutletPipeFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
