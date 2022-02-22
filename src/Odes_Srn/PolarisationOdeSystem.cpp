/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#include <utility>      // std::pair
#include <iostream>
#include <math.h>
#include "PolarisationOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

PolarisationOdeSystem::PolarisationOdeSystem(double celli_theta,
                                             double celli_phi,
                                             double neighbouring_theta,
                                             double neighbouring_phi) //,
                                             // std::vector<double> stateVariables)
    : AbstractOdeSystem(4),
    mCelli_theta(celli_theta),
    mCelli_phi(celli_phi),
    mNeighbouring_theta(neighbouring_theta),
    mNeighbouring_phi(neighbouring_phi)
    {

    mpSystemInfo.reset(new CellwiseOdeSystemInformation<PolarisationOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Theta concentration for this cell
     * 1 - Phi concentration for this cell
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    // SetDefaultInitialCondition(0, 1.0); // soon overwritten
    // SetDefaultInitialCondition(1, 1.0); // soon overwritten

    // cell1 theta = 0.209;
    // cell1 phi = 0.295;
    // GetNeighbouringTheta();
    // GetNeighbouringPhi();

    SetDefaultInitialCondition(0, celli_theta); // soon overwritten
    SetDefaultInitialCondition(1, celli_phi); // soon overwritten
    SetDefaultInitialCondition(2, neighbouring_theta); // soon overwritten
    SetDefaultInitialCondition(3, neighbouring_phi); // soon overwritten
    // SetDefaultInitialCondition(4, prod); // soon overwrittendouble prod,

    // this->mParameters.push_back(0.0);

    // // if (stateVariables != std::vector<double>())
    // // {
    //     stateVariables.push_back(celli_theta);
    //     stateVariables.push_back(celli_phi);
    //     stateVariables.push_back(neighbouring_theta);
    //     stateVariables.push_back(neighbouring_phi);
    //
    //     SetStateVariables(stateVariables);
    // // }
}

PolarisationOdeSystem::~PolarisationOdeSystem()
{
}

void PolarisationOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
//(double time, std::pair<double, double> pol_a, std::pair<double, double> r_ij, const std::vector<double>& rY, std::vector<double>& rDY)
{
    // double notch = rY[0];
    // double delta = rY[1];
    // double mean_delta = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Delta");"
    //
    // // The next two lines define the ODE system by Collier et al. (1996)
    // rDY[0] = mean_delta*mean_delta/(0.01 + mean_delta*mean_delta) - notch;  // d[Notch]/dt
    // rDY[1] = 1.0/(1.0 + 100.0*notch*notch) - delta;                   // d[Delta]/dt

    double theta = rY[0];
    double phi = rY[1];
    double neighbouringTheta = rY[2];
    double neighbouringPhi = rY[3];

    // double mean_polarisation = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Polatisation");"
    // double neighbouring_theta = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Theta");"
    // double neighbouring_phi = this->mParameters[1]; // Shorthand for "this->mParameter("Mean Phi");"
    double prod = sin(theta) * sin(neighbouringTheta) * cos(phi - neighbouringPhi) +
           cos(theta) * cos(neighbouringTheta);
    // The next two lines define the ODE system by Germann et al. (2019)
    rDY[0] = (cos(theta) * sin(neighbouringTheta) * cos(phi - neighbouringPhi) - sin(theta) * cos(neighbouringTheta))*(-prod);

    double sin_Xi_theta = sin(theta);
    if (fabs(sin_Xi_theta) > 1e-10)
    {
        rDY[1] = (-sin(neighbouringTheta) * sin(phi - neighbouringPhi) / sin_Xi_theta)*(-prod);
    }

    rDY[2] = 0.0;
    rDY[3] = 0.0;

}

template<>
void CellwiseOdeSystemInformation<PolarisationOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Theta");
    this->mVariableUnits.push_back("non-dim");
    // this->mInitialConditions.push_back(mCelli_theta);
    this->mInitialConditions.push_back(0.0); // will be filled in later
    // this->mInitialConditions.push_back(0.209); // will be filled in later

    this->mVariableNames.push_back("Phi");
    this->mVariableUnits.push_back("non-dim");
    // this->mInitialConditions.push_back(mCelli_phi);
    this->mInitialConditions.push_back(0.0); // will be filled in later
    // this->mInitialConditions.push_back(0.295); // will be filled in later

    // If this is ever not the first parameter change the line
    // double mean_Polarisation = this->mParameters[0]; in EvaluateYDerivatives().

    // this->mParameterNames.push_back("Mean Polarisation");
    // this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("Neighbour Theta");
    this->mParameterUnits.push_back("non-dim");
    // this->mInitialConditions.push_back(mNeighbouring_theta);
    this->mInitialConditions.push_back(0.0); // will be filled in later
    // this->mInitialConditions.push_back(neighbouring_theta);

    this->mParameterNames.push_back("Neighbour Phi");
    this->mParameterUnits.push_back("non-dim");
    // this->mInitialConditions.push_back(mNeighbouring_phi);
    this->mInitialConditions.push_back(0.0); // will be filled in later
    // this->mInitialConditions.push_back(neighbouring_phi);

    this->mInitialised = true;
}

double PolarisationOdeSystem::GetCelliTheta() const
{
    return mCelli_theta;
}
double PolarisationOdeSystem::GetCelliPhi() const
{
    return mCelli_phi;
}

double PolarisationOdeSystem::GetNeighbouringTheta() const
{
    return mNeighbouring_theta;
}
double PolarisationOdeSystem::GetNeighbouringPhi() const
{
    return mNeighbouring_phi;
}



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PolarisationOdeSystem)
