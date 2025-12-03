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
    Foam::turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField

Description
    This boundary condition provides a turbulence dissipation, \f$\epsilon\f$
    (epsilon) inlet condition based on a specified mixing length.  The patch
    values are calculated using:

        \f[
            \epsilon_p = \frac{C_{\mu}^{0.75} k^{1.5}}{L D_{h}}
        \f]

    where

    \vartable
        \epsilon_p | patch epsilon values
        C_{\mu} | Model coefficient, set to 0.09
        k       | turbulence kinetic energy
        L       | length scale
        D_h     | hydraulic Diameter
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        mixingLength | Length scale [m]        | yes         |
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
        Cmu          | Turbulence model coefficient | no | 0.09
        hydraulicDiameter | hydraulic diameter for the patch | no | - 
        shape        | geometry for hydraulic diameter (circle|square|rectangle) | no | circle
        shapeRatio   | (for rectangle only) a/b where a<=b (0<ratio<=1) | required if shape rectangle | -
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentPipeMixingLengthDissipationRateInlet;
        mixingLength    0.005;
        value           uniform 200;   // placeholder
    }
    \endverbatim

    Note:
        In the event of reverse flow, a zero-gradient condition is applied

See also
    Foam::inletOutletPipeFvPatchField

SourceFiles
    turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField_H
#define turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField_H

#include "inletOutletPipeFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
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
    TypeName("turbulentPipeMixingLengthDissipationRateInlet");


    // Constructors

        //- Construct from patch, internal field and dictionary
        turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
        //  onto a new patch
        turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField&,
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
                new turbulentPipeMixingLengthDissipationRateInletFvPatchScalarField
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
