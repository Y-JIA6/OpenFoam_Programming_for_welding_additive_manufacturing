Info<< "Calculating laser volumetric heat source Q\n" << endl;

// Current physical time [s]
const scalar timeNow = runTime.value();

// Find the patch where the laser is applied
const label patchID = mesh.boundaryMesh().findPatchID(laserPatchName);
if (patchID == -1)
{
    FatalErrorInFunction
        << "Patch \"" << laserPatchName << "\" not found in mesh!" << nl
        << exit(FatalError);
}

// --- Geometric information of this patch ---
const labelUList& patchFaceCells = mesh.boundary()[patchID].faceCells();  // Adjacent cells for each boundary face
const vectorField& patchCf = mesh.boundary()[patchID].Cf();               // Face centers
const scalarField& patchmagSf = mesh.boundary()[patchID].magSf();         // Face areas
// const vectorField& patchn = mesh.boundary()[patchID].nf();             // Face outward normals (if needed)

// --- Global cell geometric information ---
// const vectorField& cellC = mesh.C();                                   // Cell centers
const scalarField& cellVol = mesh.V();                                   // Cell volumes

// --- Total power and accumulated energy ---
scalar Q_total = 0.0;        // Instantaneous laser power at current time step [W]
static scalar E_accum = 0.0; // Accumulated deposited energy [J]

// --- Time control: deactivate laser outside defined time window ---
if (timeNow < startTime || timeNow > endTime)
{
    // If laser is off, reset Q field to zero while keeping correct dimensions
    Q = dimensionedScalar("Q", Q.dimensions(), 0.0);
}
else
{
    // Current laser center position = initial position + velocity × (t - startTime)
    const vector laserCenter = laserOrigin + vLaser * (timeNow - startTime);

    // Gaussian beam peak surface heat flux [W/m^2]
    const scalar q0 = 2.0 * absorption * P
                    / (Foam::constant::mathematical::pi * sqr(r0));

    // Convert to volumetric heat source [W/m^3]
    const scalar Q0 = q0;

    // Reset Q to zero to avoid residual values from previous time steps
    Q = dimensionedScalar("Q", Q.dimensions(), 0.0);

    // Loop over all faces of the selected laser patch
    forAll(patchFaceCells, fI)
    {
        const label cellI = patchFaceCells[fI];
        // const vector cellCenter = cellC[cellI];
        const vector faceCenter = patchCf[fI];
        // const vector normal = patchn[fI];

        // --- Compute lateral Gaussian distribution ---
        // r² = (x - x₀)² + (y - y₀)² assuming laser moves in x-y plane
        const scalar r2 = sqr(faceCenter.x() - laserCenter.x())
                        + sqr(faceCenter.y() - laserCenter.y());

        // Volumetric heat source for the cell [W/m^3]
        // Qval = q0 * (beam profile) * (irradiated face area) / (cell volume)
        const scalar Qval =
            Q0 * patchmagSf[fI] * Foam::exp(-2.0 * r2 / sqr(r0)) / cellVol[fI];

        Q[cellI] = Qval;

        // Power accumulation: Qval [W/m^3] × volume [m^3] = [W]
        Q_total += Qval * cellVol[cellI];
    }

    // --- Compute deposited energy for this time step ---
    const scalar deltaT = runTime.deltaTValue();
    const scalar E_step = Q_total * deltaT;  // Energy in Joules [J]
    E_accum += E_step;

    Info<< "  Laser power at t=" << timeNow << " s : " << Q_total << " W" << nl
        << "  Energy deposited in this step: " << E_step << " J" << nl
        << "  Accumulated total energy     : " << E_accum << " J" << endl;
}
