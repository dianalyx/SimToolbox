/**
 * @file Sphere.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief Sphere type
 * @version 1.0
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef SPHERE_HPP_
#define SPHERE_HPP_

#include "FDPS/particle_simulator.hpp"
#include "Util/EigenDef.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/GeoCommon.h"
#include "Util/IOHelper.hpp"

#include <cstdio>
#include <cstdlib>
#include <type_traits>
#include <vector>

/**
 * @brief Sphere class
 *
 */
class Sphere {
  public:
    int gid = GEO_INVALID_INDEX;         ///< unique global id
    int globalIndex = GEO_INVALID_INDEX; ///< unique global index sequentially ordered
    int rank = -1;                       ///< mpi rank
    double radius;                       ///< radius
    double radiusCollision;              ///< radius for collision resolution
    double radiusSearch;                 ///< radiusSearch for short range interactions
    double sepmin;                       ///< minimal separation with its neighbors within radiusSearch

    double pos[3];         ///< position
    double orientation[4]; ///< orientation quaternion. direction norm vector = orientation * (0,0,1)
    double vel[3];         ///< velocity
    double omega[3];       ///< angular velocity

    // these are not packed and transferred
    double velCol[3];     ///< collision velocity
    double omegaCol[3];   ///< collision angular velocity
    double velBrown[3];   ///< Brownian velocity
    double omegaBrown[3]; ///< Brownian angular velocity
    double velNonB[3];    ///< all non-Brownian deterministic velocity beforce collision resolution
    double omegaNonB[3];  ///< all non-Brownian deterministic angular velocity before collision resolution

    /**
     * @brief Construct a new Sphere object
     *
     */
    Sphere() = default;

    /**
     * @brief Destroy the Sphere object
     *
     */
    ~Sphere() = default;

    /**
     * @brief Construct a new Sphere object
     *
     * @param gid_
     * @param radius_
     * @param radiusCollision_
     * @param pos_ if not specified position is set as [0,0,0]
     * @param orientation_ if not specied orientation is set as identity
     */
    Sphere(const int &gid_, const double &radius_, const double &radiusCollision_, const double pos_[3] = nullptr,
           const double orientation_[4] = nullptr);

    /**
     * @brief Copy constructor
     *
     */
    Sphere(const Sphere &) = default;
    Sphere(Sphere &&) = default;
    Sphere &operator=(const Sphere &) = default;
    Sphere &operator=(Sphere &&) = default;

    /**
     * @brief display the data fields for this Sphere
     *
     */
    void dumpSphere() const;

    /**
     * @brief set the velocity data fields to zero
     *
     */
    void clear();

    /**
     * @brief update the position and orientation with internal velocity data fields and given dt
     *
     * @param dt
     */
    void stepEuler(double dt);

    /**
     * @brief return position
     *
     * necessary interface for InteractionManager.hpp
     * @return const double*
     */
    const double *Coord() const { return pos; }
    /**
     * @brief return search radius
     *
     * necessary interface for InteractionManager.hpp
     * @return double
     */
    double Rad() const { return radiusCollision * 4; }

    /**
     * @brief Get the Gid
     *
     * necessary interface for FDPS FullParticle class
     * @return int
     */
    int getGid() const { return gid; }

    /**
     * @brief Get position as a PS::F64vec3 object
     *
     * necessary interface for FDPS FullParticle class
     * @return PS::F64vec3
     */
    PS::F64vec3 getPos() const { return PS::F64vec3(pos[0], pos[1], pos[2]); }

    /**
     * @brief Set position with given PS::F64vec3 object
     *
     * necessary interface for FDPS FullParticle class
     * @param newPos
     */
    void setPos(const PS::F64vec3 &newPos) {
        pos[0] = newPos.x;
        pos[1] = newPos.y;
        pos[2] = newPos.z;
    }

    /**
     * @brief write to a file*
     *
     * FDPS IO interface
     * @param fptr
     */
    void writeAscii(FILE *fptr) const;

    /**
     * @brief write VTK XML PVTP Header file from rank 0
     *
     * @param prefix
     * @param postfix
     * @param nProcs
     */
    static void writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs);

    /**
     * @brief write VTK XML binary base64 VTP data file from every MPI rank
     *
     * Procedure for dumping Spheres in the system:
     * Each Sphere writes a polyline with two (connected) points.
     * Points are labeled with float -1 and 1
     * Sphere data fields are written as cell data
     * Rank 0 writes the parallel header , then each rank write its own serial vtp/vtu file
     *
     * @tparam Container container for local Spheres which supports [] operator
     * @param Sphere
     * @param SphereNumber
     * @param prefix
     * @param postfix
     * @param rank
     */
    template <class Container>
    static void writeVTP(const Container &Sphere, const int SphereNumber, const std::string &prefix,
                         const std::string &postfix, int rank) {
        // for each Sphere:

        // write VTP for basic data
        //  use float to save some space
        // point and point data
        std::vector<double> pos(6 * SphereNumber); // position always in Float64
        std::vector<float> label(2 * SphereNumber);

        // point connectivity of line
        std::vector<int32_t> connectivity(2 * SphereNumber);
        std::vector<int32_t> offset(SphereNumber);

        // Sphere data
        std::vector<int> gid(SphereNumber);
        std::vector<float> radius(SphereNumber);
        std::vector<float> radiusCollision(SphereNumber);
        std::vector<float> vel(3 * SphereNumber);
        std::vector<float> omega(3 * SphereNumber);
        std::vector<float> velCol(3 * SphereNumber);
        std::vector<float> omegaCol(3 * SphereNumber);
        std::vector<float> velNonB(3 * SphereNumber);
        std::vector<float> omegaNonB(3 * SphereNumber);
        std::vector<float> velBrown(3 * SphereNumber);
        std::vector<float> omegaBrown(3 * SphereNumber);
        std::vector<float> xnorm(3 * SphereNumber);
        std::vector<float> znorm(3 * SphereNumber);

#pragma omp parallel for
        for (int i = 0; i < SphereNumber; i++) {
            const auto &sy = Sphere[i];
            // point and point data
            Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
            Evec3 end0 = ECmap3(sy.pos) - direction * (sy.length * 0.5);
            Evec3 end1 = ECmap3(sy.pos) + direction * (sy.length * 0.5);
            pos[6 * i] = end0[0];
            pos[6 * i + 1] = end0[1];
            pos[6 * i + 2] = end0[2];
            pos[6 * i + 3] = end1[0];
            pos[6 * i + 4] = end1[1];
            pos[6 * i + 5] = end1[2];
            label[2 * i] = -1;
            label[2 * i + 1] = 1;

            // connectivity
            connectivity[2 * i] = 2 * i;         // index of point 0 in line
            connectivity[2 * i + 1] = 2 * i + 1; // index of point 1 in line
            offset[i] = 2 * i + 2;               // offset is the end of each line. in fortran indexing

            // Sphere data
            gid[i] = sy.gid;
            radius[i] = sy.radius;
            radiusCollision[i] = sy.radiusCollision;

            Evec3 nx = ECmapq(sy.orientation) * Evec3(1, 0, 0);
            Evec3 nz = ECmapq(sy.orientation) * Evec3(0, 0, 1);
            for (int j = 0; j < 3; j++) {
                vel[3 * i + j] = sy.vel[j];
                omega[3 * i + j] = sy.omega[j];
                velBrown[3 * i + j] = sy.velBrown[j];
                omegaBrown[3 * i + j] = sy.omegaBrown[j];
                velCol[3 * i + j] = sy.velCol[j];
                omegaCol[3 * i + j] = sy.omegaCol[j];
                velNonB[3 * i + j] = sy.velNonB[j];
                omegaNonB[3 * i + j] = sy.omegaNonB[j];

                xnorm[3 * i + j] = nx[j];
                znorm[3 * i + j] = nz[j];
            }
        }

        std::ofstream file(prefix + std::string("Sphere_") + "r" + std::to_string(rank) + std::string("_") + postfix +
                               std::string(".vtp"),
                           std::ios::out);

        IOHelper::writeHeadVTP(file);

        file << "<Piece NumberOfPoints=\"" << SphereNumber * 2 << "\" NumberOfLines=\"" << SphereNumber << "\">\n";
        // Points
        file << "<Points>\n";
        IOHelper::writeDataArrayBase64(pos, "position", 3, file);
        file << "</Points>\n";
        // cell definition
        file << "<Lines>\n";
        IOHelper::writeDataArrayBase64(connectivity, "connectivity", 1, file);
        IOHelper::writeDataArrayBase64(offset, "offsets", 1, file);
        file << "</Lines>\n";
        // point data
        file << "<PointData Scalars=\"scalars\">\n";
        IOHelper::writeDataArrayBase64(label, "endLabel", 1, file);
        file << "</PointData>\n";
        // cell data
        file << "<CellData Scalars=\"scalars\">\n";
        IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
        IOHelper::writeDataArrayBase64(radius, "radius", 1, file);
        IOHelper::writeDataArrayBase64(radiusCollision, "radiusCollision", 1, file);
        IOHelper::writeDataArrayBase64(vel, "velocity", 3, file);
        IOHelper::writeDataArrayBase64(omega, "omega", 3, file);
        IOHelper::writeDataArrayBase64(velBrown, "velocityBrown", 3, file);
        IOHelper::writeDataArrayBase64(omegaBrown, "omegaBrown", 3, file);
        IOHelper::writeDataArrayBase64(velCol, "velocityCollision", 3, file);
        IOHelper::writeDataArrayBase64(omegaCol, "omegaCollision", 3, file);
        IOHelper::writeDataArrayBase64(velNonB, "velocityNonB", 3, file);
        IOHelper::writeDataArrayBase64(omegaNonB, "omegaNonB", 3, file);
        IOHelper::writeDataArrayBase64(xnorm, "xnorm", 3, file);
        IOHelper::writeDataArrayBase64(znorm, "znorm", 3, file);
        file << "</CellData>\n";
        file << "</Piece>\n";

        IOHelper::writeTailVTP(file);
        file.close();
    }
};

/**
 * @brief FDPS writeAscii file header
 */
class SphereAsciiHeader {
  public:
    int nparticle;
    double time;
    void writeAscii(FILE *fp) const { fprintf(fp, "%d \n %lf\n", nparticle, time); }
};

static_assert(std::is_trivially_copyable<Sphere>::value, "");
static_assert(std::is_default_constructible<Sphere>::value, "");

#endif