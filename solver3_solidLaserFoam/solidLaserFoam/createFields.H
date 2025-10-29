Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",                     // Field name: "T"
        runTime.timeName(),      // Time directory (e.g., "0", "1", etc.)
        mesh,                    // Associated mesh
        IOobject::MUST_READ,     // Must be read from file
        IOobject::AUTO_WRITE     // Automatically written at each time step
    ),
    mesh
);


// --- physicalProperties ---
Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",    // Dictionary name
        runTime.constant(),      // File located in "constant/"
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,  // Re-read if modified during runtime
        IOobject::NO_WRITE
    )
);


// Read while preserving physical dimensions (do not use .value())
Info<< "Reading density rho\n" << endl;
dimensionedScalar rho( physicalProperties.lookup("rho") );

Info<< "Reading thermal conductivity k\n" << endl;
dimensionedScalar k( physicalProperties.lookup("k") );

Info<< "Reading specific heat Cp\n" << endl;
dimensionedScalar Cp( physicalProperties.lookup("Cp") );


// Optionally create a combined field (choose either scalar or field version, not both)
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
    dimensionedScalar
    (
        "rhoCp",
        rho.dimensions()*Cp.dimensions(),  // Combined dimensions
        rho.value()*Cp.value()             // Combined numerical value
    )
);


// --- volumetric heat source Q: W/m^3 ---
Info<< "Creating volumetric heat source field Q\n" << endl;
volScalarField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ ,   // ⚠️ Usually set to NO_READ if Q is generated internally
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Q", dimEnergy/dimTime/dimVolume, 0.0) // Units: [W/m^3]
);


// --- laserProperties ---
Info<< "Reading laserProperties\n" << endl;

IOdictionary laserProperties
(
    IOobject
    (
        "laserProperties",       // File name: constant/laserProperties
        runTime.constant(),      // Located in constant/
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,     // Must exist
        IOobject::NO_WRITE
    )
);


// Laser parameters
const scalar P          = readScalar(laserProperties.lookup("P"));          // Laser power [W]
const scalar r0         = readScalar(laserProperties.lookup("r0"));         // Beam radius [m]
const scalar absorption = readScalar(laserProperties.lookup("absorption")); // Absorption coefficient [-]
const vector vLaser     (laserProperties.lookup("vLaser"));                 // Laser motion velocity [m/s]
const vector laserOrigin(laserProperties.lookup("laserOrigin"));            // Laser starting position [m]


// Time window during which the laser is active
const scalar startTime =
    laserProperties.found("startTime") ? readScalar(laserProperties.lookup("startTime")) : 0.0;
const scalar endTime =
    laserProperties.found("endTime")   ? readScalar(laserProperties.lookup("endTime"))   : GREAT;


// Name of the patch that receives laser heating
const word laserPatchName
(
    laserProperties.lookup("patchName")   // Read the laser heating patch name
);


// --- Dimensional consistency check (printed once for verification) ---
Info<< "dims rho="   << rho.dimensions()
    << " Cp="        << Cp.dimensions()
    << " k="         << k.dimensions()
    << " rhoCp="     << rhoCp.dimensions()
    << " Q="         << Q.dimensions() << nl;


// Initialize fvConstraints (e.g., for boundary conditions or fixed fields)


#include "createFvModels.H"
#include "createFvConstraints.H"
