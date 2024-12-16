/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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

#include "ibOmegaWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedInternalValueFvPatchField.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::ibOmegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::ibOmegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<ibOmegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            ibOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void Foam::ibOmegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

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

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<ibOmegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            for (const auto& celli : faceCells)
            {
                ++weights[celli];
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    for (const auto& patchi : omegaPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::ibOmegaWallFunctionFvPatchScalarField&
Foam::ibOmegaWallFunctionFvPatchScalarField::omegaPatch
(
    const label patchi
)
{
    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const auto& opf =
        refCast<const ibOmegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<ibOmegaWallFunctionFvPatchScalarField&>(opf);
}


void Foam::ibOmegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            ibOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            ibOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void Foam::ibOmegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
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
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    // Set omega and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];
        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar w = cornerWeights[facei];

        // Contribution from the viscous sublayer
        const scalar omegaVis = 6.0*nuw[facei]/(beta1_*sqr(y[facei]));

        // Contribution from the inertial sublayer
        const scalar omegaLog = sqrt(k[celli])/(Cmu25*kappa*y[facei]);

        switch (blender_)
        {
            case blenderType::STEPWISE:
            {
                if (yPlus > yPlusLam)
                {
                    omega0[celli] += w*omegaLog;
                }
                else
                {
                    omega0[celli] += w*omegaVis;
                }
                break;
            }

            case blenderType::BINOMIAL:
            {
                omega0[celli] +=
                    w*pow
                    (
                        pow(omegaVis, n_) + pow(omegaLog, n_),
                        scalar(1)/n_
                    );
                break;
            }

            case blenderType::MAX:
            {
                // (PH:Eq. 27)
                omega0[celli] += max(omegaVis, omegaLog);
                break;
            }

            case blenderType::EXPONENTIAL:
            {
                // (PH:Eq. 31)
                const scalar Gamma = 0.01*pow4(yPlus)/(1 + 5*yPlus);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);

                omega0[celli] +=
                    w*(omegaVis*exp(-Gamma) + omegaLog*exp(-invGamma));
                break;
            }

            case blenderType::TANH:
            {
                // (KAS:Eqs. 33-34)
                const scalar phiTanh = tanh(pow4(0.1*yPlus));
                const scalar b1 = omegaVis + omegaLog;
                const scalar b2 =
                    pow(pow(omegaVis, 1.2) + pow(omegaLog, 1.2), 1.0/1.2);

                omega0[celli] += phiTanh*b1 + (1 - phiTanh)*b2;
                break;
            }
        }

        if (!(blender_ == blenderType::STEPWISE) || yPlus > yPlusLam)
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


void Foam::ibOmegaWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
    os.writeEntryIfDifferent<scalar>("beta1", 0.075, beta1_);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibOmegaWallFunctionFvPatchScalarField::ibOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    wallFunctionBlenders(),
    initialised_(false),
    master_(-1),
    beta1_(0.075),
    wallCoeffs_(),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::ibOmegaWallFunctionFvPatchScalarField::ibOmegaWallFunctionFvPatchScalarField
(
    const ibOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    initialised_(false),
    master_(-1),
    beta1_(ptf.beta1_),
    wallCoeffs_(ptf.wallCoeffs_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::ibOmegaWallFunctionFvPatchScalarField::ibOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::BINOMIAL, scalar(2)),
    initialised_(false),
    master_(-1),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    wallCoeffs_(dict),
    G_(),
    omega_(),
    cornerWeights_()
{
    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::ibOmegaWallFunctionFvPatchScalarField::ibOmegaWallFunctionFvPatchScalarField
(
    const ibOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf),
    wallFunctionBlenders(owfpsf),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    wallCoeffs_(owfpsf.wallCoeffs_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::ibOmegaWallFunctionFvPatchScalarField::ibOmegaWallFunctionFvPatchScalarField
(
    const ibOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf, iF),
    wallFunctionBlenders(owfpsf),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    wallCoeffs_(owfpsf.wallCoeffs_),
    G_(),
    omega_(),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::ibOmegaWallFunctionFvPatchScalarField::G
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

    return omegaPatch(master_).G();
}


Foam::scalarField& Foam::ibOmegaWallFunctionFvPatchScalarField::omega
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void Foam::ibOmegaWallFunctionFvPatchScalarField::updateCoeffs()
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
    vectorField n = patch().nf();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa_ = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    typedef DimensionedField<scalar, volMesh> FieldType;
    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());
    scalarField& omega = refValue();
    forAll (nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus = Cmu25*y[faceI]*sqrt(k[faceCellI])/nuw[faceI];

        scalar omegaVis = 6.0*nuw[faceI]/(beta1_*sqr(y[faceI]));

        scalar omegaLog = sqrt(k[faceCellI])/(Cmu25*kappa_*y[faceI]);
        
        omega[faceI] = sqrt(sqr(omegaVis) + sqr(omegaLog));

        if (yPlus > yPlusLam)
        {
            G[faceCellI] =
                (nutw[faceI] + nuw[faceI])
               *magGradUw[faceI]
               *Cmu25*sqrt(k[faceCellI])
               /(kappa_*y[faceI]);
        }
        else
        {
            G[faceCellI] = 0.0;
        }
    }
    fixedInternalValueFvPatchField<scalar>::updateCoeffs();    
}


void Foam::ibOmegaWallFunctionFvPatchScalarField::updateWeightedCoeffs
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
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}

void Foam::ibOmegaWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        ibOmegaWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
