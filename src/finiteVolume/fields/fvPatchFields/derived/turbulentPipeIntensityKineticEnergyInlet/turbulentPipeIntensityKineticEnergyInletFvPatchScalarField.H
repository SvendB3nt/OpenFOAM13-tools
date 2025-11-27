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
    Foam::turbulentPipeIntensityKineticEnergyInletFvPatchScalarField

Description
    This boundary condition provides a turbulent kinetic energy condition
    for pipe flows. The turbulence intensity is based on either 
    user-supplied hydraulic diameter of the patch, or shape of the patch.
    The Reynolds number for the patched is based on area averaged values
    of U and nu. 

        \f[
            I = 0.16 Re^(-1/8)
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        Re      | reynolds number for the patch
        I       | turbulence intensity
        k_p     | kinetic energy at the patch
        U       | velocity field
    \endvartable

    In the event of reverse flow, a zero-gradient condition is applied.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        nu           | kinematic viscosity field name     | no          | nu
        phi          | flux field name         | no          | phi
        hydraulicDiameter | hydraulic diameter for the patch | no | - 
        shape        | geometry for hydraulic diameter (circle|square|rectangle) | no | circle
        shapeRatio   | (for rectangle only) a/b where a<=b (0<ratio<=1) | required if shape rectangle | -
    \endtable

    Examples of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        turbulentPipeIntensityKineticEnergyInlet;
        value       uniform 1;      // placeholder 
    }

    <patchName>
    {
        type        turbulentPipeIntensityKineticEnergyInlet;
        value       uniform 1;      // placeholder
        shape       rectangle;
        shapeRatio  0.5;            // a/b (a is shortest side)
    }

    <patchName>
    {
        type        turbulentPipeIntensityKineticEnergyInlet;
        value       uniform 1;      // placeholder
        hydraulicDiameter 1; 
    }

    \endverbatim

See also
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentPipeIntensityKineticEnergyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentPipeIntensityKineticEnergyInletFvPatchScalarField_H
#define turbulentPipeIntensityKineticEnergyInletFvPatchScalarField_H

#include "inletOutletPipeFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class turbulentPipeIntensityKineticEnergyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
:
    public inletOutletPipeFvPatchScalarField
{
    // Private Data
        //- Name of the velocity field
        word UName_;

        //- Name of the kinematic viscosity field
        word nuName_;

 

public:

    //- Runtime type information
    TypeName("turbulentPipeIntensityKineticEnergyInlet");


    // Constructors

        //- Construct from patch, internal field and dictionary
        turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
        //  onto a new patch
        turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
        (
            const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
        (
            const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
        (
            const turbulentPipeIntensityKineticEnergyInletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentPipeIntensityKineticEnergyInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
