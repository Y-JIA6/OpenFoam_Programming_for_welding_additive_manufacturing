Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// --- physicalProperties ---
Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading density rho\n" << endl;
dimensionedScalar rho( physicalProperties.lookup("rho") );

Info<< "Reading thermal conductivity k\n" << endl;
dimensionedScalar k( physicalProperties.lookup("k") );

Info<< "Reading specific heat Cp\n" << endl;
dimensionedScalar Cp( physicalProperties.lookup("Cp") );


volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCp", rho.dimensions()*Cp.dimensions(), rho.value()*Cp.value())
);


#include "createFvModels.H"
#include "createFvConstraints.H"
