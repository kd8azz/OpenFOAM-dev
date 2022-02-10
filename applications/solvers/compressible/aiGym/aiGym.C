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
    aiGym

Description
    AI Gym based on rhoPimpleFoam.  Agents exist as point entities and
    interact with their environment by modifying air pressure directly.
    This is intended as a cheap and plausible analogue for audio sensing
    and production.

\*---------------------------------------------------------------------------*/

#include <cmath>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
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
    double x = 0;
    double y = 0;
    Coords(double x_, double y_): x(x_), y(y_) {}
    Coords() {}
    Coords& operator+=(Coords other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    std::pair<bool, bool> Clamp(int32_t min_x, int32_t max_x, int32_t min_y, int32_t max_y) {
        std::pair<bool, bool> clamped;
        if (x < min_x) {x = min_x; clamped.first = true; }
        if (x > max_x - 1) { x = max_x - 1; clamped.first = true; }
        if (y < min_y) { y = min_y; clamped.second = true; }
        if (y > max_y - 1) { y = max_y - 1; clamped.second = true; }
        return clamped;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & x;
        ar & y;
    }
};
template <typename Number>
Coords operator*(const Coords& coords, Number number) {
    Coords result = coords;
    result.x *= number;
    result.y *= number;
    return result;
}
struct Actor {
    int32_t id;
    Coords position;
    double theta = 0.d;
    double speed = 0.d;
    double friction = 0.d;

    // aiGym-side only.
    double previous_pressure = 0.d;
    double previous_theta = 0.d;
    double previous_speed = 0.d;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & id;
        ar & position;
        ar & theta;
        ar & speed;
        ar & friction;
    }
};
// Both the input and output, for a given Actor's experience.
// Input: a request to do the following.
// Output: what actually happened.  E.g. if the actor wants to
//         accelerate into a wall, they're going to hit it instead.
struct ActorStateChange {
    int32_t id = 0;
    // +/- from the average pressure at the start of the simulation.
    double pressure_change = 0.d;
    double theta_change = 0.f;
    double speed_change = 0.f;

    // output only.
    Coords position;
    Coords velocity;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & id;
        ar & pressure_change;
        ar & theta_change;
        ar & speed_change;
        ar & position;
        ar & velocity;
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
    Info<< "argc " << argc << "\n";
    assert(argc == 2);
    std::string port(argv[1]);
    {
        int p = std::stoi(port);
        assert(p > 5000);
        assert(p < 8000);
    }
    std::string url = "tcp://localhost:";
    url.append(port);
    argc = 1;

    Info<< "Opening socket to AI controller.\n";
    zmq::context_t context{1};
    zmq::socket_t socket{context, zmq::socket_type::req};
    socket.connect(url);
    Info<< "Socket open on this end. (" << url << ")\n";

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
    auto getRho = [&](double x, double y) -> double {
      int cell = width * static_cast<int>(y) + static_cast<int>(x);
      return rho[cell];
    };
    auto setRho  = [&](double x, double y, double value) {
      int cell = width * static_cast<int>(y) + static_cast<int>(x);
      rho[cell] = value;
    };

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
            ActionRequest request;
            request.changes.reserve(actors.size());
            for (auto& pair : actors) {
                Actor& actor = pair.second;
                Coords direction = Coords(/*x=*/std::cos(actor.theta), /*y=*/std::sin(actor.theta));
                actor.position += direction * actor.speed * runTime.deltaTValue();
                std::pair<bool, bool> clamped = actor.position.Clamp(0, width, 0, height);
                if (clamped.first) {
                    if (direction.y > 0) actor.theta = M_PI / 2;
                    else actor.theta = 3 * M_PI / 2;
                }
                if (clamped.second) {
                    if (direction.x > 0) actor.theta = 0;
                    else actor.theta = M_PI;
                }
                if (clamped.first && clamped.second) actor.speed = 0.d;

                ActorStateChange state_change;
                state_change.id = actor.id;
                state_change.pressure_change =
                    getRho(actor.position.x, actor.position.y) - actor.previous_pressure;
                state_change.position = actor.position;
                state_change.theta_change = actor.theta - actor.previous_theta;
                state_change.speed_change = actor.speed - actor.previous_speed;
                request.changes.push_back(state_change);

                actor.previous_pressure =
                    getRho(actor.position.x, actor.position.y);
                actor.previous_theta = actor.theta;
                actor.previous_speed = actor.speed;
            }

            std::ostringstream oss;
            boost::archive::binary_oarchive oa(oss);
            oa << request;
            socket.send(zmq::buffer(oss.str()), zmq::send_flags::none);
        }

        {
            zmq::message_t reply{};
            socket.recv(reply, zmq::recv_flags::none);
            std::istringstream iss(reply.to_string());
            boost::archive::binary_iarchive ia(iss);
            ActionResponse response;
            ia >> response;

            for (const int32_t erase : response.erase_actors) {
                actors.erase(erase);
            }
            for (Actor& actor: response.new_actors) {
                actors[actor.id] = actor;
            }
            for (const ActorStateChange& change : response.changes) {
                Actor& actor = actors[change.id];

                double value = getRho(actor.position.x, actor.position.y);
                value += change.pressure_change;
                setRho(actor.position.x, actor.position.y, value);

                actor.theta += change.theta_change;
                actor.speed += change.speed_change;
            }
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
