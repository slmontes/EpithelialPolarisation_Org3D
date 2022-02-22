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

#include "PolarisationSrnModel.hpp"

PolarisationSrnModel::PolarisationSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(4, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<PolarisationSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        // mpOdeSolver = CellCycleModelOdeSolver<PolarisationSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        EulerIvpOdeSolver mpOdeSolver
        mpOdeSolver->Initialise();
        SetDt(0.01);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

PolarisationSrnModel::PolarisationSrnModel(const PolarisationSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    assert(rModel.GetOdeSystem());
    // SetOdeSystem(new PolarisationOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
    // SetOdeSystem(new PolarisationOdeSystem(rModel.GetOdeSystem()->GetInitialConditions());

}

AbstractSrnModel* PolarisationSrnModel::CreateSrnModel()
{
    return new PolarisationSrnModel(*this);
}

void PolarisationSrnModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdatePolarisation();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void PolarisationSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new PolarisationOdeSystem);

    double celli_theta = mpCell->GetCellData()->GetItem("Theta");
    double celli_phi = mpCell->GetCellData()->GetItem("Phi");
    double neighbouring_theta = 0.0;
    double neighbouring_phi = 0.0;

    mpOdeSystem = new PolarisationOdeSystem(celli_theta,celli_phi,neighbouring_theta, neighbouring_phi);
    // mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

}

// void PolarisationSrnModel::AdjustOdeParameters(double currentTime)
// {
//   // if (mpOdeSystem != nullptr)
//   // {
//   //   mpCell->GetCellData()->SetItem("Theta", 0.793);//("Theta", 0.209);
//   //   mpCell->GetCellData()->SetItem("Phi", 0.073);//("Phi", 0.295);
//   // }
//   // else
//   // {
//   // Pass this time step's the cell's current values for Theta and Phi into the solver as a constant over this timestep.
//     mpOdeSystem->rGetStateVariables()[0] = mpCell->GetCellData()->GetItem("Theta");
//     mpOdeSystem->rGetStateVariables()[1] = mpCell->GetCellData()->GetItem("Phi");
//
//   // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
//   // mpOdeSystem->rGetStateVariables()[2] = this->GetNeighbouringTheta();
//   // mpOdeSystem->rGetStateVariables()[3] = this->GetNeighbouringPhi();
//   // }
// }

// void SingleOdeWntCellCycleModel::AdjustOdeParameters(double currentTime)
// {
//     // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
//     mpOdeSystem->rGetStateVariables()[2] = this->GetWntLevel();
//
//     // Use the cell's current mutation status as another input
//     static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
// }

void PolarisationSrnModel::UpdatePolarisation()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i=0; i<4; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }

    double i_theta = mpOdeSystem->rGetStateVariables()[0];
    double i_phi = mpOdeSystem->rGetStateVariables()[1];

    double j_theta = mpOdeSystem->rGetStateVariables()[2];
    double j_phi = mpOdeSystem->rGetStateVariables()[3];

    // Pass this value to the CellData Theta and Phi.
    // mpCell->GetCellData()->SetItem("Theta", new_theta);
    // mpCell->GetCellData()->SetItem("Phi", new_phi);

    // Set those values for next ODE run
    mpOdeSystem->SetParameter("Theta", i_theta);
    mpOdeSystem->SetParameter("Phi", i_phi);

    mpOdeSystem->SetParameter("Neighbour Theta", j_theta);
    mpOdeSystem->SetParameter("Neighbour Phi", j_phi);

}

double PolarisationSrnModel::GetTheta()
{
    assert(mpOdeSystem != nullptr);
    double theta = mpOdeSystem->rGetStateVariables()[0];
    return theta;
}

double PolarisationSrnModel::GetPhi()
{
    assert(mpOdeSystem != nullptr);
    double phi = mpOdeSystem->rGetStateVariables()[1];
    return phi;
}



// double PolarisationSrnModel::GetMeanNeighbouringPolarisation()
// {
//     assert(mpOdeSystem != nullptr);
//     double mean_neighbouring_Polarisation = mpOdeSystem->GetParameter("Mean Polarisation");
//     return mean_neighbouring_Polarisation;
// }
// double PolarisationSrnModel::GetNeighbouringTheta()
// {
//     assert(mpOdeSystem != nullptr);
//     double neighbouring_theta = mpOdeSystem->GetParameter("Neighbour Theta");
//     return neighbouring_theta;
// }
//
// double PolarisationSrnModel::GetNeighbouringPhi()
// {
//     assert(mpOdeSystem != nullptr);
//     double neighbouring_phi = mpOdeSystem->GetParameter("Neighbour Phi");
//     return neighbouring_phi;
// }

void PolarisationSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PolarisationSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(PolarisationSrnModel)
