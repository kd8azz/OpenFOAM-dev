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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <map>
#include <sstream>
#include <zmq.hpp>

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

// 2D for now.
struct Coords {
    int32_t x = 0;
    int32_t y = 0;
    Coords& operator+=(Coords other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    void Clamp(int32_t min_x, int32_t max_x, int32_t min_y, int32_t max_y) {
        if (x < min_x) x = min_x;
        if (x >= max_x) x = max_x - 1;
        if (y < min_y) y = min_y;
        if (y >= max_y) y = max_y - 1;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & x;
        ar & y;
    }
};
struct Actor {
    int32_t id;
    Coords position;
    Coords velocity;
    // aiGym-side only.
    double last_observed_relative_pressure = 0.d;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & id;
        ar & position;
        ar & velocity;
    }
};
// Both the input and output, for a given Actor's experience.
// Input: a request to do the following.
// Output: what actually happened.  E.g. if the actor wants to
//         accelerate into a wall, they're going to hit it instead.
struct ActorStateChange {
    int32_t id = 0;
    // +/- from average pressure
    double relative_pressure = 0.d;
    Coords velocity_change;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & id;
        ar & relative_pressure;
        ar & velocity_change;
    }
};
struct ActionRequest {
    double timestep_duration = 0.d;
    std::vector<ActorStateChange> changes;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & timestep_duration;
        ar & changes;
    }
};
struct ActionResponse {
    std::vector<int32_t> erase_actors;
    std::vector<Actor> new_actors;
    std::vector<ActorStateChange> changes;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & new_actors;
        ar & erase_actors;
        ar & changes;
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

    double average_pressure = 0.0d;
    for (int x = 0; x < width; ++x) {
        for(int y = 0; y < height; ++y) {
            int cell = width * y + x;
            average_pressure += rho[cell];
        }
    }
    average_pressure /= width * height;

    std::map<int32_t, Actor> actors;
    while (pimple.run(runTime))
    {
        {
            Info<< "<Gym Step>\n";
            {
                ActionRequest request;
                request.changes.reserve(actors.size());
                Info<< "Processing " << actors.size() << " actors.\n";
                for (auto& pair : actors) {
                    Actor& actor = pair.second;
                    actor.position += actor.velocity;
                    actor.position.Clamp(0, width, 0, height);
                    int cell = width * actor.position.y + actor.position.x;

                    ActorStateChange state_change;
                    state_change.id = actor.id;
                    state_change.relative_pressure = rho[cell] - average_pressure;
                    actor.last_observed_relative_pressure = rho[cell] - average_pressure;
                    request.changes.push_back(state_change);
                }
                Info<< "Done processing actors.\n";

                Info<< "Encoding.\n";
                std::ostringstream oss;
                boost::archive::text_oarchive oa(oss);
                oa << request;
                Info<< "Sending.\n";
                socket.send(zmq::buffer(oss.str()), zmq::send_flags::none);
                Info<< "Sent.\n";
            }

            {
                Info<< "Receiving.\n";
                zmq::message_t reply{};
                socket.recv(reply, zmq::recv_flags::none);
                Info<< "Decoding.\n";
                std::istringstream iss(reply.to_string());
                boost::archive::text_iarchive ia(iss);
                ActionResponse response;
                ia >> response;
                Info<< "Decoded.\n";

                Info<< "Processing "
                    << response.erase_actors.size() << " erasures, "
                    << response.new_actors.size() << " creations, and "
                    << response.changes.size() << " changes.\n";
                for (const int32_t erase : response.erase_actors) {
                    actors.erase(erase);
                }
                for (Actor& actor: response.new_actors) {
                    int cell = width * actor.position.y + actor.position.x;
                    actors[actor.id] = actor;
                }
                for (const ActorStateChange& change : response.changes) {
                    Actor& actor = actors[change.id];
                    int cell = width * actor.position.y + actor.position.x;
                    rho[cell] += change.relative_pressure - actor.last_observed_relative_pressure;
                    actor.velocity += change.velocity_change;
                }
            }
            Info<< "</Gym Step>.\n";
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
