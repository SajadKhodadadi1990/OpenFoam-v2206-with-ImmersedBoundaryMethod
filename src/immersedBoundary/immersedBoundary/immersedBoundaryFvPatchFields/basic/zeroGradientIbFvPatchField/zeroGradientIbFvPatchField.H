/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::zeroGradientIbFvPatchField

Description
    Foam::zeroGradientIbFvPatchField

Author
    Hrvoje Jasak

SourceFiles
    zeroGradientIbFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroGradientIbFvPatchField_H
#define zeroGradientIbFvPatchField_H

#include "zeroGradientFvPatchFields.H"
#include "immersedBoundaryFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class zeroGradientIbFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroGradientIbFvPatchField
:
    public zeroGradientFvPatchField<Type>,
    public immersedBoundaryFieldBase<Type>
{
    // Private Member Functions

        //- Update IB value and gradient
        void updateIbValues();


public:

    //- Runtime type information
    TypeName("zeroGradientIb");


    // Constructors

        //- Construct from patch and internal field
        zeroGradientIbFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroGradientIbFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given zeroGradientIbFvPatchField
        //  onto a new patch
        zeroGradientIbFvPatchField
        (
            const zeroGradientIbFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroGradientIbFvPatchField
        (
            const zeroGradientIbFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new zeroGradientIbFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroGradientIbFvPatchField
        (
            const zeroGradientIbFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new zeroGradientIbFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~zeroGradientIbFvPatchField()
    {}


    // Member functions

        // Mapping functions

            // //- Map (and resize as needed) from self given a mapping object
            // virtual void autoMap
            // (
            //     const fvPatchFieldMapper&
            // );

            // //- Reverse map the given fvPatchField onto this fvPatchField
            // virtual void rmap
            // (
            //     const fvPatchField<Type>&,
            //     const labelList&
            // );

            //- Update on mesh motion
            virtual void updateOnMotion();


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =  Pstream::commsTypes::blocking
            );

            //- Manipulate a matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "zeroGradientIbFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
