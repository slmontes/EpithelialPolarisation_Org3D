/*
MODIFIED BY Sandra M 01/08/19
*/

#ifndef EPITHELIALLAYERPOLARISATIONFORCE_WMATRIGEL_HPP_
#define EPITHELIALLAYERPOLARISATIONFORCE_WMATRIGEL_HPP_

#include "AbstractForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <math.h>


//-----------------------------------------------------------------------------
template<unsigned SPACE_DIM>
class EpithelialLayerPolarisationForce_wMatrigel : public AbstractForce<SPACE_DIM>
{

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<SPACE_DIM> >(*this);
    }


public:

    /**
     * Constructor.
     */
    EpithelialLayerPolarisationForce_wMatrigel();

    /**
     * Destructor.
     */
    virtual ~EpithelialLayerPolarisationForce_wMatrigel();

    /*
     * c_vector mutiplied by a double
     */
    c_vector<double, SPACE_DIM> cvector_multiplied(c_vector<double, SPACE_DIM> vec, double number);

    /*
     *Get 3D vector from polarity angles
     */
     c_vector<double, SPACE_DIM> pol_to_cvector(double cell_theta, double cell_phi);

    /*
     *Get polarity angles from 3D vector and it's norm
     */
     std::pair<double, double> pt_to_pol(c_vector<double, SPACE_DIM> r, double dist);

    /*
     *Get polarity angles from a vector
     */
     std::pair<double, double> pt_to_pol(c_vector<double, SPACE_DIM> r);

    /*
     *Get dot product from polarity angles
     */
     double pol_dot_product(double celli_theta, double celli_phi, double cellj_theta, double cellj_phi);

    /*
     *Aligning force from the potential U = - Σ(p_i . p_j), such that all
     *polarities point in the same direction.
     */
     std::pair<double, double> unidirectional_polarization_force(double celli_theta,
                                                                 double celli_phi,
                                                                 double cellj_theta,
                                                                 double cellj_phi);

    /*
     *Aligning force from the potential U_Pol = - Σ(p_i . p_j)^2/2, such
     *that all polarities are oriented the same way.
     */
     std::pair<double, double> bidirectional_polarization_force(double celli_theta,
                                                                double celli_phi,
                                                                double cellj_theta,
                                                                double cellj_phi);

    /*
     * Solver for polarisation angles
     */
     std::pair<double, double> EulerSolver_PolarisationAngles(double celli_theta,
                                                              double celli_phi,
                                                              double neighbour_theta,
                                                              double neighbour_phi);

    /*
     * Resistance to bending from the potential U_Epi = Σ(p_i . r_ij/r)^2/2.
     */
     c_vector<double, SPACE_DIM> bending_force(CellPtr p_celli,
                                               double celli_theta,
                                               double celli_phi,
                                               double cellj_theta,
                                               double cellj_phi,
                                               c_vector<double, SPACE_DIM> unit_vector,
                                               double dij);

    /*
     * Calculates the force between two nodes.
     */
    // c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(c_vector<double, SPACE_DIM> r_node_i_location,
    //                                                        c_vector<double, SPACE_DIM> unit_vector);

    void AddForceContribution(AbstractCellPopulation<SPACE_DIM>& rCellPopulation);

    // /**
    //  * Overridden OutputForceParameters() method.
    //  *
    //  * @param rParamsFile the file stream to which the parameters are output
    //  */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EpithelialLayerPolarisationForce_wMatrigel)

#endif /*EPITHELIALLAYERPOLARISATIONFORCE_WMATRIGEL_HPP_*/
