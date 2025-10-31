// -----------------------------------------------------------------------------
// updateMaterialProperties.H
// Temperature-dependent material properties via table interpolation
// Handles out-of-range by clamping to min/max table values.
// -----------------------------------------------------------------------------

// --- Load tables once (static) ---
static List<Tuple2<scalar, scalar>> kTable, CpTable, rhoTable;
static word outOfRange = "clamp";
static bool tablesLoaded = false;

if (!tablesLoaded)
{
    Info<< "Loading temperature-dependent material properties..." << endl;
    physicalProperties.lookup("kTable")   >> kTable;
    physicalProperties.lookup("CpTable")  >> CpTable;
    physicalProperties.lookup("rhoTable") >> rhoTable;
    outOfRange = physicalProperties.lookupOrDefault<word>("outOfRange", "clamp");

    // Sort once
    Foam::sort(kTable,  [](const auto& a, const auto& b){ return a.first() < b.first(); });
    Foam::sort(CpTable, [](const auto& a, const auto& b){ return a.first() < b.first(); });
    Foam::sort(rhoTable,[](const auto& a, const auto& b){ return a.first() < b.first(); });

    tablesLoaded = true;

    if (Pstream::master()) 
{
    Info<< "\n--- Loaded material property tables ---" << nl;

    Info<< "kTable (Temperature [K] , k [W/m/K])" << endl;
    forAll(kTable, i)
    {
        Info<< "  " << kTable[i].first() << "  " << kTable[i].second() << endl;
    }

    Info<< "CpTable (Temperature [K] , Cp [J/kg/K])" << endl;
    forAll(CpTable, i)
    {
        Info<< "  " << CpTable[i].first() << "  " << CpTable[i].second() << endl;
    }

    Info<< "rhoTable (Temperature [K] , rho [kg/m^3])" << endl;
    forAll(rhoTable, i)
    {
        Info<< "  " << rhoTable[i].first() << "  " << rhoTable[i].second() << endl;
    }

    Info<< "----------------------------------------" << nl << endl;
}
}




// --- Linear interpolation helper ---
auto linearInterp = [](scalar x0, scalar y0, scalar x1, scalar y1, scalar x)
{
    if (mag(x1 - x0) < SMALL) return 0.5 * (y0 + y1);
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
};

// --- Sample property with out-of-range handling ---
auto sampleProperty = [&](const List<Tuple2<scalar, scalar>>& t, scalar Tval) -> scalar
{
    if (t.empty()) return 0.0;

    // Below min
    if (Tval <= t.first().first())
    {
        return t.first().second(); // clamp to min value
    }

    // Above max
    if (Tval >= t.last().first())
    {
        return t.last().second(); // clamp to max value
    }

    // Interpolate between table points
    for (label i = 1; i < t.size(); ++i)
    {
        if (Tval < t[i].first())
        {
            return linearInterp
            (
                t[i-1].first(), t[i-1].second(),
                t[i].first(),   t[i].second(),
                Tval
            );
        }
    }

    return t.last().second(); // fallback
};

// --- Update fields per cell ---
forAll(T, i)
{
    const scalar Ti = T[i];

    k[i]   = max(sampleProperty(kTable,  Ti), SMALL);
    Cp[i]  = max(sampleProperty(CpTable, Ti), SMALL);
    rho[i] = max(sampleProperty(rhoTable,Ti), SMALL);
}

forAll(T.boundaryField(), patchI)
{
    const fvPatchScalarField& Tp = T.boundaryField()[patchI];
    fvPatchScalarField& kp  = k.boundaryFieldRef()[patchI];
    fvPatchScalarField& cpp = Cp.boundaryFieldRef()[patchI];
    fvPatchScalarField& rhop= rho.boundaryFieldRef()[patchI];

    forAll(Tp, faceI)
    {
        const scalar Ti = Tp[faceI];
        kp[faceI]   = max(sampleProperty(kTable,  Ti), SMALL);
        cpp[faceI]  = max(sampleProperty(CpTable, Ti), SMALL);
        rhop[faceI] = max(sampleProperty(rhoTable,Ti), SMALL);
    }
}

// --- Combine fields
rhoCp = rho * Cp;
