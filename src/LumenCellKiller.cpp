/*
* Created on: Feb 24 2021
* Author: Sandra Montes
* Cell killer created for epithelial layer model. Removes any epithelial cells
* that are located in the lumen of the simulation.
*
* Modified on: Jul 09 2021
* Last modified:
* 			Author: Sandra Montes
*/

#include "LumenCellKiller.hpp"
#include "NodeBasedCellPopulation.hpp"


LumenCellKiller::LumenCellKiller(AbstractCellPopulation<3>* pCellPopulation)
        : AbstractCellKiller<3>(pCellPopulation)
        {
          // Sets up output file
          //	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
          //	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
        }
//           mNumberOfNeighbours(NumberOfNeighbours)
// {
//     if ((mNumberOfNeighbours<4) || (mNumberOfNeighbours>10))
//     {
//         EXCEPTION("Number of neighbours of the cell must be between five and ten");
//     }
// }


// double LumenCellKiller::GetNumberOfNeighbours() const
// {
//     return mNumberOfNeighbours;
// }

LumenCellKiller::~LumenCellKiller()
{
//    mAnoikisOutputFile->close();
}


void LumenCellKiller::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    /*
     * We retrieve a vector with the indices of all neighbouring indices according
     * to X distance (Node based cell population interaction radius), calculate
     * the size of the vector and label the cell for appoptosis if the number of
     * neighbours is greater than the defined amount of mNumberOfNeighbours
     */

     NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*> (this->mpCellPopulation);

 		//Update cell population
 		p_tissue->Update();

 		double radius = 1; //1.5

    // Get node index corresponding to cell
    unsigned cell_node_index = mpCellPopulation->GetLocationIndexUsingCell(pCell);

 		std::set<unsigned> neighbouring_node_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(cell_node_index, radius);

    if (!pCell->HasApoptosisBegun() &&
        neighbouring_node_indices.size() < 3)
    {
        pCell->StartApoptosis();
    }
}


void LumenCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        CheckAndLabelSingleCellForApoptosis(*cell_iter);
    }
}


void LumenCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // *rParamsFile << "\t\t\t<NumberOfNeighbours>" << mNumberOfNeighbours << "</NumberOfNeighbours>\n";
    // *rParamsFile << "\t\t\t<RealNumberOfNeighbours>" << real_num_neighbours << "</RealNumberOfNeighbours>\n";

    // Call method on direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(LumenCellKiller)
