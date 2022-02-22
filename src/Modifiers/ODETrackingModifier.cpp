/*
 * PolarisationTrackingModifier.cpp
 *
 *  Created on: 10 Mar 2021
 *      Author: Sandra M
 */

#include "ODETrackingModifier.hpp"

// #include "HippoDynamicWntCellCycleModel.hpp"

template<unsigned DIM>
ODETrackingModifier<DIM>::ODETrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{}

template<unsigned DIM>
ODETrackingModifier<DIM>::~ODETrackingModifier()
{}

template<unsigned DIM>
void ODETrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODETrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODETrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Loop over all of the cells within the population and fetch the relevant Nanog concentrations from the NanogInheritedCellCycleModel
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	// // Create a static cast to the cell cycle model in order to fetch the nanog concentration to be allocated to celldata
    	// HippoDynamicWntCellCycleModel* p_model = static_cast<HippoDynamicWntCellCycleModel*>(cell_iter->GetCellCycleModel());
      //   //double this_e2f1 = p_model->GetE2F1Level();
      //   double volume = p_model->GetCelliTheta();
      //   // Note that the state variables must be in the same order as listed in PolarisationOdeSystem
      //   cell_iter->GetCellData()->SetItem("volume_in_cell", volume);
    }


}

template<unsigned DIM>
void ODETrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ODETrackingModifier<1>;
template class ODETrackingModifier<2>;
template class ODETrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ODETrackingModifier)
