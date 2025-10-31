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

// --- Temperature-dependent property fields
volScalarField k
(
    IOobject
    (
        "k", 
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.0)
);

volScalarField Cp
(
    IOobject
    (
        "Cp", 
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,2,-2,-1,0,0,0), 0.0)
);

volScalarField rho
(
    IOobject
    (
        "rho", 
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-3,0,0,0,0,0), 0.0)
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp", 
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-1,-2,-1,0,0,0), 0.0)
);

// --- Reading material properties dictionary ---
IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// --- Laser source term ---
Info<< "Creating volumetric heat source field Q\n" << endl;

volScalarField Q
(
    IOobject
    (
        "Q", 
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Q", dimEnergy/dimTime/dimVolume, 0.0)
);

// --- Laser properties ---
Info<< "Reading laserProperties\n" << endl;

IOdictionary laserProperties
(
    IOobject("laserProperties", runTime.constant(), mesh, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE)
);

const scalar P          = readScalar(laserProperties.lookup("P"));
const scalar r0         = readScalar(laserProperties.lookup("r0"));
const scalar absorption = readScalar(laserProperties.lookup("absorption"));
const vector vLaser     (laserProperties.lookup("vLaser"));
const vector laserOrigin(laserProperties.lookup("laserOrigin"));

const scalar startTime =
    laserProperties.found("startTime") ? readScalar(laserProperties.lookup("startTime")) : 0.0;
const scalar endTime =
    laserProperties.found("endTime")   ? readScalar(laserProperties.lookup("endTime"))   : GREAT;

const word laserPatchName(laserProperties.lookup("patchName"));

#include "createFvModels.H"
#include "createFvConstraints.H"
