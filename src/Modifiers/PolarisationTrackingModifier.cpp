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

#include "PolarisationTrackingModifier.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "PolarisationOdeSystem.hpp"
// #include "PolarisationSrnModel.hpp"

template<unsigned DIM>
PolarisationTrackingModifier<DIM>::PolarisationTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PolarisationTrackingModifier<DIM>::~PolarisationTrackingModifier()
{
}

template<unsigned DIM>
void PolarisationTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarisationTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarisationTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    ofstream myfile;
    myfile.open ("polarity_test.txt");

    // First recover each cell's phi and Theta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
         {
             EulerIvpOdeSolver euler_solver;

             CellPtr celli = *cell_iter;

             //Get current polarisation angles for cell i
             double celli_theta = celli->GetCellData()->GetItem("Theta");
             double celli_phi = celli->GetCellData()->GetItem("Phi");

             // Get the set of neighbouring location indices
              std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(celli);

                 // Compute this cell's average neighbouring Polarisation concentration and store in CellData
                 if (!neighbour_indices.empty())
                 {
                     for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                          iter != neighbour_indices.end();
                          ++iter)
                     {
                       CellPtr cellj = rCellPopulation.GetCellUsingLocationIndex(*iter);

                       double celli_theta = celli->GetCellData()->GetItem("Theta");
                       double celli_phi = celli->GetCellData()->GetItem("Phi");

                       double cellj_theta = cellj->GetCellData()->GetItem("Theta");
                       double cellj_phi = cellj->GetCellData()->GetItem("Phi");

                       PolarisationOdeSystem ode_system(celli_theta,celli_phi,cellj_theta,cellj_phi);
                       std::vector<double> initial_conditions = ode_system.GetInitialConditions();

                       // double time_step = 0.25;
                       double time_step = SimulationTime::Instance()->GetTimeStep();
                       // double time = Timer::GetElapsedTime();
                       double time = SimulationTime::Instance()->GetTime();

                       std::vector<double> current_y_values=initial_conditions;
                       std::vector<double> next_y_values(initial_conditions.size());

                       next_y_values= current_y_values;

                       OdeSolution solutions;

                       // Ensure the internal data structures are the right size
                       solutions = euler_solver.Solve(&ode_system, current_y_values, time, time+time_step, time_step, 1e-5);
                        // solution = mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);

                       int num_timesteps = solutions.GetNumberOfTimeSteps();

                       double new_Theta = solutions.rGetSolutions()[num_timesteps-1][0];
                       double new_Phi = solutions.rGetSolutions()[num_timesteps-1][1];

                       celli->GetCellData()->SetItem("Theta", new_Theta);
                       celli->GetCellData()->SetItem("Phi", new_Phi);

                         myfile << "Writing this to a file.\n";
                         myfile << "Number of time steps: " << num_timesteps << std::endl;
                         myfile << "Init Cell i Theta: " << celli_theta << std::endl;
                         myfile << "Init Cell i Phi: " << celli_phi << std::endl;
                         myfile << "Init Cell j Theta: " << cellj_theta << std::endl;
                         myfile << "Init Cell j Phi: " << cellj_phi << std::endl;
                         myfile << "new_Theta: " << celli->GetCellData()->GetItem("Theta") << std::endl;
                         myfile << "new_Phi: " << celli->GetCellData()->GetItem("Phi") << std::endl;


                     }
                   }

                  myfile.close();

         }

    // {
    //     PolarisationSrnModel* p_model = static_cast<PolarisationSrnModel*>(cell_iter->GetSrnModel());
    //     double this_theta = p_model->GetTheta();
    //     double this_phi = p_model->GetPhi();
    //
    //     // Note that the state variables must be in the same order as listed in PolarisationOdeSystem
    //     cell_iter->GetCellData()->SetItem("Theta", this_theta);
    //     cell_iter->GetCellData()->SetItem("Phi", this_phi);
    //
    //     std::cout << "Theta: " << this_theta << std::endl;
    //     std::cout << "Phi: " << this_phi << std::endl;
    //
    // }

    // // Next iterate over the population to compute and store each cell's neighbouring Polarisation concentration in CellData
    // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     // Get the set of neighbouring location indices
    //     std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
    //
    //     // Compute this cell's average neighbouring Polarisation concentration and store in CellData
    //     if (!neighbour_indices.empty())
    //     {
    //         for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
    //              iter != neighbour_indices.end();
    //              ++iter)
    //         {
    //             // CellPtr p_celli = rCellPopulation.GetCellUsingLocationIndex(*cell_iter);
    //             // CellPtr p_cellj = rCellPopulation.GetCellUsingLocationIndex(*iter);
    //             //
    //             // double celli_theta = p_celli->GetCellData()->GetItem("Theta");
    //             // double celli_phi = p_celli->GetCellData()->GetItem("Phi");
    //             //
    //             // double cellj_theta = p_cellj->GetCellData()->GetItem("Theta");
    //             // double cellj_phi = p_cellj->GetCellData()->GetItem("Phi");
    //             //mean_Polarisation += std::make_pair(this_theta/neighbour_indices.size(), this_phi/neighbour_indices.size());
    //
    //             PolarisationSrnModel* p_model = static_cast<PolarisationSrnModel*>(cell_iter->GetSrnModel());
    //
    //             // accPol_celli_theta += this_theta/neighbour_indices.size();
    //             // accPol_celli_phi += this_phi/neighbour_indices.size();
    //         }
    //         // cell_iter->GetCellData()->SetItem("mean Polarisation", mean_Polarisation);
    //         // cell_iter->GetCellData()->SetItem("Theta", accPol_celli_theta);
    //         // cell_iter->GetCellData()->SetItem("Phi", accPol_celli_phi);
    //     }
    //     // else
    //     // {
    //     //     // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
    //     //     // cell_iter->GetCellData()->SetItem("mean Polarisation", 0.0);
    //     //     cell_iter->GetCellData()->SetItem("mean theta", 0.0);
    //     //     cell_iter->GetCellData()->SetItem("mean phi", 0.0);
    //     // }
    // }
}

template<unsigned DIM>
void PolarisationTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);

}

// Explicit instantiation
template class PolarisationTrackingModifier<1>;
template class PolarisationTrackingModifier<2>;
template class PolarisationTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarisationTrackingModifier)
