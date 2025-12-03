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

Class
    Foam::turbulentPipeMixingLengthFrequencyInletFvPatchScalarField

Description
    This boundary condition provides a turbulence specific dissipation,
    \f$\omega\f$ (omega) inlet condition based on the turbulent mixing length.
    The 
    The patch values are calculated using:

        \f[
            \omega_p = \frac{k^{0.5}}{C_{\mu}^{0.25} L D_{h}}
        \f]

    where

    \vartable
        \omega_p | patch omega values
        C_{\mu} | Model coefficient, set to 0.09
        k       | turbulence kinetic energy
        L       | turbulent length scale, set to 0.07
        D_{h}   | hydraulic diameter
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
        Cmu          | Turbulence model coefficient | no | 0.09
        L            | turbulent length scale  | no | 0.07
        hydraulicDiameter | hydraulic diameter for the patch | no | - 
        shape        | geometry for hydraulic diameter (circle|square|rectangle) | no | circle
        shapeRatio   | (for rectangle only) a/b where a<=b (0<ratio<=1) | required if shape rectangle | -
    \endtable

    Examples of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentPipeMixingLengthFrequencyInlet;
        value           uniform 200;   // placeholder
    }

    \endverbatim

    Note:
        In the event of reverse flow, a zero-gradient condition is applied

See also
    Foam::inletOutletPipeFvPatchField

SourceFiles
    turbulentPipeMixingLengthFrequencyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentPipeMixingLengthFrequencyInletFvPatchScalarField_H
#define turbulentPipeMixingLengthFrequencyInletFvPatchScalarField_H

#include "inletOutletPipeFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class turbulentPipeMixingLengthFrequencyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
:
    public inletOutletPipeFvPatchScalarField
{
    // Private Data

        
        //- Turbulent length scale
        scalar turbulentLengthScale_;

        //- Name of the turbulent kinetic energy field
        word kName_;

        //- Turbulence model coefficient
        scalar Cmu_;


public:

    //- Runtime type information
    TypeName("turbulentPipeMixingLengthFrequencyInlet");


    // Constructors

        //- Construct from patch, internal field and dictionary
        turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
        //  onto a new patch
        turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
        (
            const turbulentPipeMixingLengthFrequencyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
        (
            const turbulentPipeMixingLengthFrequencyInletFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
        (
            const turbulentPipeMixingLengthFrequencyInletFvPatchScalarField&,
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
                new turbulentPipeMixingLengthFrequencyInletFvPatchScalarField
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
