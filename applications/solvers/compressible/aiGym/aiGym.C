/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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

Application
    rhoPimpleFoam

Description
    Transient solver for turbulent flow of compressible fluids for HVAC and
    similar applications, with optional mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include <cmath>
#include <zmq.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermophysicalTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct pcm_out {
    double pressure = 0.d;
    int movement_x = 0;
    int movement_y = 0;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & pressure;
        ar & movement_x;
        ar & movement_y;
    }
};
struct pcm_in {
    double pressure = 0.d;
    int requested_movement_x = 0;
    int requested_movement_y = 0;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & pressure;
        ar & requested_movement_x;
        ar & requested_movement_y;
    }
};

int main(int argc, char *argv[])
{
    Info<< "Opening socket to AI controller.\n";
    zmq::context_t context{1};
    zmq::socket_t socket{context, zmq::socket_type::req};
    socket.connect("tcp://localhost:5555");
    Info<< "Socket open on this end.\n";

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop with deltaT = " << runTime.deltaTValue() << "\n" << endl;

    int width = 300;
    int height = 300;
    int x = 100;
    int y = 20;
    int prev_x = x;
    int prev_y = y;
    while (pimple.run(runTime))
    {
        int cell = width * y + x;

        {
            std::ostringstream oss;
            boost::archive::text_oarchive oa(oss);
            pcm_out out;
            out.pressure = rho[cell];
            out.movement_x = x - prev_x;
            out.movement_y = y - prev_y;
            oa << out;
            socket.send(zmq::buffer(oss.str()), zmq::send_flags::none);
        }

        prev_x = x;
        prev_y = y;

        {
            zmq::message_t reply{};
            socket.recv(reply, zmq::recv_flags::none);
            std::istringstream iss(reply.to_string());
            boost::archive::text_iarchive ia(iss);
            pcm_in in;
            ia >> in;
            rho[cell] = in.pressure;
            x += in.requested_movement_x;
            if (x < 0) x = 0;
            if (x >= width) x = width - 1;
            y += in.requested_movement_y;
            if (y < 0) y = 0;
            if (y >= height) y = height - 1;
        }

        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", rho*U);
                }

                fvModels.preUpdateMesh();

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if
            (
                !mesh.steady()
             && !pimple.simpleRho()
             && pimple.firstPimpleIter()
            )
            {
                #include "rhoEqn.H"
            }

            fvModels.correct();

            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        if (!mesh.steady())
        {
            rho = thermo.rho();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
