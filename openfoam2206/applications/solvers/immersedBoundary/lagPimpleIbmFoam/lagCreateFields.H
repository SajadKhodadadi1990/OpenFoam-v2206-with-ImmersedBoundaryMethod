dimensionedScalar nu
(
    "nu",
    dimViscosity,
    laminarTransport
);

dimensionedScalar rhoVal
(
    "rho",
    dimDensity,
    laminarTransport
);

    // dummy-fields to satisfy the requirements of the kinematicCloud
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoVal,
        "zeroGradient"
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoVal*nu,
        "zeroGradient"
    );

#   include "readGravitationalAcceleration.H"

    Info<< "Constructing kinematicCloud" << endl;
    basicKinematicCloud kinematicCloud
    (
        "kinematicCloud",
        rho,
        U,
        mu,
        g
    );
    volScalarField ParticleMassInIbPatchCell
    (
                IOobject
                (
                    "ParticleMassInIbPatchCell",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass, 0.0)
    );
    label particleCounter=0;

