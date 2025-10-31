/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    solidLaserFoam

Description
    Transient heat conduction solver:
    rho*Cp*dT/dt = div(k*grad(T)) + source
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"


    Info<< "\nSolving transient heat conduction\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();          

        #include "laser.H"          // laser input 
//        Info<< "Q dimensions: " << Q.dimensions() << endl;

        #include "updateMatProperties.H" 


        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)      // ✅ ρCp·dT/dt
            - fvm::laplacian(k, T)       // ✅ ∇·(k∇T)
            == fvModels.source(rhoCp, T)
        );
        TEqn.source() += Q;         // incoprate the laser energy

        fvConstraints.constrain(TEqn);
        TEqn.solve();
        fvConstraints.constrain(T);


        runTime.write();
        


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //




