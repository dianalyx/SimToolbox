/**
 * @file SphereConfig.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief read configuration parameters from a yaml file
 * @version 1.0
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef SPHERECONFIG_HPP_
#define SPHERECONFIG_HPP_

#include <iostream>

/**
 * @brief read configuration parameters from a yaml file
 *
 */
class SphereConfig {
  public:
    // parallel setting
    int ompThreads = -1;
    unsigned int rngSeed = 0;

    // domain setting
    double simBoxHigh[3] = {1.0, 1.0, 1.0};    // simulation box size
    double simBoxLow[3] = {0, 0, 0};           // simulation box size
    bool simBoxPBC[3] = {false, false, false}; // flag of true/false of periodic in that direction
    bool wallLowZ = false;
    bool wallHighZ = false;

    double initBoxHigh[3] = {1.0, 1.0, 1.0}; // initial size
    double initBoxLow[3] = {0, 0, 0};        // initial size
    double initOrient[3] = {2.0, 2.0, 2.0};  // initial orientation for each sylinder. >1 <-1 means random
    bool initCircularX = false;              // set the initial cross-section as a circle

    // physical constant
    double viscosity = 1.0; // pN/(um^2 s)
    double KBT = 0.00411;   // pN.um

    int sphereNumber = 100;
    double sphereRadius = 2.0;
    double sphereRadiusSigma = 0;
    double sphereRadiusColRatio = 1.0;

    // time stepping
    double dt = 0.01;
    double timeTotal = 0.01;
    double timeSnap = 0.01;

    // collision solver
    double colResTol = 0.00001;
    int colMaxIte = 8000;
    bool colNewtonRefine = false;

    SphereConfig() = default;
    SphereConfig(std::string filename);
    ~SphereConfig() = default;

    void dump() const;
};

#endif