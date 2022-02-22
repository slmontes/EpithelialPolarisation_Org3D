/*
* LAST MODIFIED: 24/02/2021
* Cell killer created for epithelial layer model. Removes any epithelial cells
* that are located in the lumen of the simulation.
*
* Created on: Feb 24 2021
* Last modified:
* 			Author: Sandra Montes
*/

#ifndef LUMENCELLKILLER_HPP_
#define LUMENCELLKILLER_HPP_

#include "CheckpointArchiveTypes.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellKiller.hpp"


/**
 * A cell killer that kills cells based on the number of surronding neighbours.
 *
 * Whenever a cell have a greater number of neighbours than the proposed number
 * CheckAndLabelCellsForApoptosis() is called.
 *
 * Note this cell killer will aim to kill all the "centre" cell to create a
 * monolayer with a lumen cleared of cells. This code was based on the
 * RandomCellKiller template.
 *
 * Previously: We assumed a constant time step and that there are an integer
 * number (n = 1/dt) of time steps per hour.
 * We also assumed that this method is called every time step
 * and that the probabilities of dying at different times are independent.
 */

class LumenCellKiller : public AbstractCellKiller<3>
{
private:

    /**
     * Set nuimber of neighbours
      */
     // double mNumberOfNeighbours;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<3> >(*this);
        // archive & mNumberOfNeighbours;
    }

public:

    // Destructor
    ~LumenCellKiller();

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();
    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param NumberOfNeighbours number of neighbours the cell has to achive to
     * be labelled for apoptosis.
     */
    LumenCellKiller(AbstractCellPopulation<3>* pCellPopulation);

    /**
     * @return mNumberOfNeighbours.
     */
    // double GetNumberOfNeighbours() const;

    /**
     * Overridden method to test a given cell for apoptosis.
     *
     * @param pCell the cell to test for apoptosis
     */
    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);

    /**
     * Loop over cells and start apoptosis acording to label
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(LumenCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a LumenCellKiller.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const LumenCellKiller * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    // double numNeighbours = t->GetNumberOfNeighbours();
    // ar << numNeighbours;
}

/**
 * De-serialize constructor parameters and initialise a LumenCellKiller.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, LumenCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_cell_population;
    ar >> p_cell_population;
    // double numNeighbours;
    // ar >> numNeighbours;

    // Invoke inplace constructor to initialise instance
    // ::new(t)LumenCellKiller(p_cell_population, numNeighbours);
    ::new(t)LumenCellKiller(p_cell_population);
}
}
} // namespace ...

#endif /*LUMENCELLKILLER_HPP_*/
