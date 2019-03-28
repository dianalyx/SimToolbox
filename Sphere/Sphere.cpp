#include "Sphere.hpp"
#include "Util/Base64.hpp"

/*****************************************************
 *  Sphero-cylinder
 ******************************************************/

Sphere::Sphere(const int &gid_, const double &radius_, const double &radiusCollision_, const double pos_[3],
               const double orientation_[4]) {
    gid = gid_;
    radius = radius_;
    radiusCollision = radiusCollision_;
    if (pos_ == nullptr) {
        Emap3(pos).setZero();
    } else {
        for (int i = 0; i < 3; i++) {
            pos[i] = pos_[i];
        }
    }
    if (orientation_ == nullptr) {
        Emapq(orientation).setIdentity();
    } else {
        for (int i = 0; i < 4; i++) {
            orientation[i] = orientation_[i];
        }
    }

    clear();
    return;
}

void Sphere::clear() {
    Emap3(vel).setZero();
    Emap3(omega).setZero();
    Emap3(velCol).setZero();
    Emap3(omegaCol).setZero();
    Emap3(velBrown).setZero();
    Emap3(omegaBrown).setZero();
    Emap3(velNonB).setZero();
    Emap3(omegaNonB).setZero();
    sepmin = std::numeric_limits<double>::max();
}

void Sphere::dumpSphere() const {
    printf("gid %d, R %g, RCol %g, L %g, LCol %g, pos %g, %g, %g\n", gid, radius, radiusCollision, pos[0], pos[1],
           pos[2]);
    printf("vel %g, %g, %g; omega %g, %g, %g\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
    printf("orient %g, %g, %g, %g\n", orientation[0], orientation[1], orientation[2], orientation[3]);
}

void Sphere::writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs) {
    std::vector<std::string> pieceNames;

    std::vector<IOHelper::FieldVTU> pointDataFields;
    pointDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "endLabel");

    std::vector<IOHelper::FieldVTU> cellDataFields;
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "radius");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "radiusCollision");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocity");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omega");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocityBrown");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omegaBrown");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocityCollision");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omegaCollision");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocityNonB");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omegaNonB");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "xnorm");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "znorm");

    for (int i = 0; i < nProcs; i++) {
        pieceNames.emplace_back(std::string("Sphere_") + std::string("r") + std::to_string(i) + "_" + postfix + ".vtp");
    }

    IOHelper::writePVTPFile(prefix + "Sphere_" + postfix + ".pvtp", pointDataFields, cellDataFields, pieceNames);
}

void Sphere::stepEuler(double dt) {
    Emap3(pos) += Emap3(vel) * dt;
    Equatn currOrient = Emapq(orientation);
    EquatnHelper::rotateEquatn(currOrient, Emap3(omega), dt);
    Emapq(orientation).x() = currOrient.x();
    Emapq(orientation).y() = currOrient.y();
    Emapq(orientation).z() = currOrient.z();
    Emapq(orientation).w() = currOrient.w();
}

void Sphere::writeAscii(FILE *fptr) const {
    fprintf(fptr, "S %d %g %g %g %g\n", gid, radius, pos[0], pos[1], pos[2]);
}
