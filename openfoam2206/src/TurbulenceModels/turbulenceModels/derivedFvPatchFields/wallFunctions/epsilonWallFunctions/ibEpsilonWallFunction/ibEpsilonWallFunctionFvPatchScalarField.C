/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "ibEpsilonWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedInternalValueFvPatchField.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::ibEpsilonWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::ibEpsilonWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<ibEpsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            ibEpsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<ibEpsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            for (const auto& faceCell : faceCells)
            {
                ++weights[faceCell];
            }
        }
    }

    cornerWeights_.setSize(bf.size());

    for (const auto& patchi : epsilonPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), Zero);
    epsilon_.setSize(internalField().size(), Zero);

    initialised_ = true;
}


Foam::ibEpsilonWallFunctionFvPatchScalarField&
Foam::ibEpsilonWallFunctionFvPatchScalarField::epsilonPatch
(
    const label patchi
)
{
    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const auto& epf =
        refCast<const ibEpsilonWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<ibEpsilonWallFunctionFvPatchScalarField&>(epf);
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            ibEpsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            ibEpsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar Cmu75 = pow(wallCoeffs_.Cmu(), 0.75);
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    // Set epsilon and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar w = cornerWeights[facei];

        // Contribution from the viscous sublayer
        const scalar epsilonVis = w*2.0*k[celli]*nuw[facei]/sqr(y[facei]);

        // Contribution from the inertial sublayer
        const scalar epsilonLog = w*Cmu75*pow(k[celli], 1.5)/(kappa*y[facei]);

        switch (blender_)
        {
            case blenderType::STEPWISE:
            {
                if (lowReCorrection_ && yPlus < yPlusLam)
                {
                    epsilon0[celli] += epsilonVis;
                }
                else
                {
                    epsilon0[celli] += epsilonLog;
                }
                break;
            }

            case blenderType::BINOMIAL:
            {
                // (ME:Eqs. 15-16)
                epsilon0[celli] +=
                    pow
                    (
                        pow(epsilonVis, n_) + pow(epsilonLog, n_),
                        scalar(1)/n_
                    );
                break;
            }

            case blenderType::MAX:
            {
                // (PH:Eq. 27)
                epsilon0[celli] += max(epsilonVis, epsilonLog);
                break;
            }

            case blenderType::EXPONENTIAL:
            {
                // (PH:p. 193)
                const scalar Gamma = 0.001*pow4(yPlus)/(scalar(1) + yPlus);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);
                epsilon0[celli] +=
                    epsilonVis*exp(-Gamma) + epsilonLog*exp(-invGamma);
                break;
            }

            case blenderType::TANH:
            {
                // (KAS:Eqs. 33-34)
                const scalar phiTanh = tanh(pow4(0.1*yPlus));
                const scalar b1 = epsilonVis + epsilonLog;
                const scalar b2 =
                    pow(pow(epsilonVis, 1.2) + pow(epsilonLog, 1.2), 1.0/1.2);

                epsilon0[celli] += phiTanh*b1 + (1 - phiTanh)*b2;
                break;
            }
        }

        if (!lowReCorrection_ || (yPlus > yPlusLam))
        {
            G0[celli] +=
                w
               *(nutw[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(k[celli])
               /(kappa*y[facei]);
        }
    }
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
    os.writeEntryIfDifferent<bool>("lowReCorrection", false, lowReCorrection_);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibEpsilonWallFunctionFvPatchScalarField::
ibEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    wallFunctionBlenders(),
    lowReCorrection_(false),
    initialised_(false),
    master_(-1),
    wallCoeffs_(),
    G_(),
    epsilon_(),
    cornerWeights_()
{

}


Foam::ibEpsilonWallFunctionFvPatchScalarField::
ibEpsilonWallFunctionFvPatchScalarField
(
    const ibEpsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    lowReCorrection_(ptf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ptf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{

}


Foam::ibEpsilonWallFunctionFvPatchScalarField::
ibEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::STEPWISE, scalar(2)),
    lowReCorrection_(dict.getOrDefault("lowReCorrection", false)),
    initialised_(false),
    master_(-1),
    wallCoeffs_(dict),
    G_(),
    epsilon_(),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::ibEpsilonWallFunctionFvPatchScalarField::
ibEpsilonWallFunctionFvPatchScalarField
(
    const ibEpsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf),
    wallFunctionBlenders(ewfpsf),
    lowReCorrection_(ewfpsf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ewfpsf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{

}


Foam::ibEpsilonWallFunctionFvPatchScalarField::
ibEpsilonWallFunctionFvPatchScalarField
(
    const ibEpsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf, iF),
    wallFunctionBlenders(ewfpsf),
    lowReCorrection_(ewfpsf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ewfpsf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::ibEpsilonWallFunctionFvPatchScalarField::G
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


Foam::scalarField& Foam::ibEpsilonWallFunctionFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const label patchi = patch().index();
    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar Cmu75 = pow(wallCoeffs_.Cmu(), 0.75);
    const scalar kappa_ = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    typedef DimensionedField<scalar, volMesh> FieldType;

    //~ volScalarField& G = const_cast<volScalarField&>
        //~ (db().lookupObject<volScalarField>(turbModel.GName()));
    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    // Note: epsilon is now a refValue and set in
    // fixedInternalValueFvPatchField
    // HJ, 3/Aug/2011
    scalarField& epsilon = refValue();
    // Set epsilon and G
    forAll (nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus = Cmu25*y[faceI]*sqrt(k[faceCellI])/nuw[faceI];

        // Note: epsilon is now a refValue and set in
        // fixedInternalValueFvPatchField
        // HJ, 3/Aug/2011
        epsilon[faceI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);
;
        if (yPlus > yPlusLam)
        {
            // Original OpenFOAM implementation
//             G[faceCellI] =
//                 (nutw[faceI] + nuw[faceI])*magGradUw[faceI]
//                *Cmu25*sqrt(k[faceCellI])/(kappa_*y[faceI]);

            // Change for consistency with Fluent implementation.
            // Emil Baric, NUMAP-FOAM 2011
            // HJ, 13/Dec/2011
            G[faceCellI] =
                sqr((nutw[faceI] + nuw[faceI])*magGradUw[faceI])/
                (Cmu25*sqrt(k[faceCellI])*kappa_*y[faceI]);
        }
        else
        {
            G[faceCellI] = 0.0;
        }
    }
    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    scalarField& epsilonf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            epsilon[celli] = (1.0 - w)*epsilon[celli] + w*epsilon0[celli];
            epsilonf[facei] = epsilon[celli];
        }
    }

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}
/*

void Foam::ibEpsilonWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

   // matrix.setValues(patch().faceCells(), patchInternalField());

    fixedInternalValueFvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::ibEpsilonWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintValues(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& fld = internalField();

    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            const label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << constraintCells.size()
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues(constraintCells, constraintValues);

    fixedInternalValueFvPatchField<scalar>::manipulateMatrix(matrix);
}

*/
void Foam::ibEpsilonWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        ibEpsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
