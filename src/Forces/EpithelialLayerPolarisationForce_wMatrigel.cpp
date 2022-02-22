/*
MODIFIED BY Sandra M 11/03/21
*/
#include <utility>      // std::pair
#include <iostream>
#include <math.h>
#include "EpithelialLayerPolarisationForce_wMatrigel.hpp"
#include "IsNan.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "UblasCustomFunctions.hpp"
#include "PolarisationOdeSystem.hpp"
#include "Timer.hpp"
#include "EulerIvpOdeSolver.hpp"


//-----------------------------------------------------------------------------
//------------------------------------------------------------------------------
template<unsigned SPACE_DIM>
EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::EpithelialLayerPolarisationForce_wMatrigel()
   : AbstractForce<SPACE_DIM>()
{
}

template<unsigned SPACE_DIM>
EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::~EpithelialLayerPolarisationForce_wMatrigel()
{
}

/*
 * c_vector mutiplied by a double
 */
template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::cvector_multiplied(c_vector<double, SPACE_DIM> vec, double number)
{
  c_vector<double, 3> multiplied_vec;
  multiplied_vec[0] = vec[0] * number;
  multiplied_vec[1] = vec[1] * number;
  multiplied_vec[2] = vec[2] * number;
  return vec;
}

/*
 * Obtains a 3Dvector from the values of theta (first) and phi (second)
 */
template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::pol_to_cvector(double cell_theta, double cell_phi)
{
  c_vector<double, 3> vec;
  vec[0] = sin(cell_theta) * cos(cell_phi);
  vec[1] = sin(cell_theta) * sin(cell_phi);
  vec[2] = cos(cell_theta);
  return vec;
}

/*
 * Obtains the values of theta and phi from a 3D vector
 */
template<unsigned SPACE_DIM>
std::pair<double, double> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::pt_to_pol(c_vector<double, SPACE_DIM> r, double dist)
{
  double theta = acos(r[2] / dist);  //r[2] is the same as the component z of the vector r (r.z)
  double phi = atan2(r[1], r[0]);   //r[1] is the same as the component y of the vector r (r.y) //2d comand angle to work with the singularity...see the 3d one*** (value rater than infinity)
                                  //r[0] is the same as the component x of the vector r (r.x)

  std::pair<double, double> pol = std::make_pair(theta, phi);
  //used std::pair<double, double> instead of Polarity (struct):
  //pol.first = theta
  //pol.second =  phi
  return pol;
}

/*
 * Obtains the values of phi and theta from the norm of r (called dist)
 * calculated by the square root of the sum of squares of three coordinates
 * of the argument.
 */
 template<unsigned SPACE_DIM>
 std::pair<double, double> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::pt_to_pol(c_vector<double, SPACE_DIM> r)
 {
   auto dist = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
   return pt_to_pol(r, dist);
 }

 /*
  *Obtains the dot product of 2 elements from their polarity angles
  */
  template<unsigned SPACE_DIM>
  double EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::pol_dot_product(double celli_theta,
                                                                        double celli_phi,
                                                                        double cellj_theta,
                                                                        double cellj_phi)
  {
    return sin(celli_theta) * sin(cellj_theta) * cos(celli_phi - cellj_phi) +
           cos(celli_theta) * cos(cellj_theta);
  }

  /*
   *Aligning force from the potential U = - Σ(p_i . p_j), such that all
   *polarities point in the same direction.
   */
   template<unsigned SPACE_DIM>
   std::pair<double, double> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::unidirectional_polarization_force(double celli_theta,
                                                                                                              double celli_phi,
                                                                                                              double cellj_theta,
                                                                                                              double cellj_phi)
   {
    //  dF.theta = cosf(Xi.theta) * sinf(p.theta) * cosf(Xi.phi - p.phi) -
    //            sinf(Xi.theta) * cosf(p.theta);
    // auto sin_Xi_theta = sinf(Xi.theta);
    // if (fabs(sin_Xi_theta) > 1e-10)
    // dF.phi = -sinf(p.theta) * sinf(Xi.phi - p.phi) / sin_Xi_theta;

     /*Lets test this theory*/
     // double celli_theta = p_celli->GetCellData()->GetItem("Theta");
     // double celli_phi = p_celli->GetCellData()->GetItem("Phi");
     //
     // double cellj_theta = p_cellj->GetCellData()->GetItem("Theta");
     // double cellj_phi = p_cellj->GetCellData()->GetItem("Phi");

     double theta = cos(celli_theta) * sin(cellj_theta) * cos(celli_phi - cellj_phi) -
                    sin(celli_theta) * cos(cellj_theta);

     double sin_Xi_theta = sin(celli_theta);
     // if (fabs(sin_Xi_theta) > 1e-10)
     // {
     double phi = -sin(cellj_theta) * sin(celli_phi - cellj_phi) / sin_Xi_theta;
     // }

     std::pair<double, double> temp_pol = std::make_pair(theta, phi);


     return temp_pol;
   }

  /*
   *Aligning force from the potential U_Pol = - Σ(p_i . p_j)^2/2, such
   *that all polarities are oriented the same way.
   */
   template<unsigned SPACE_DIM>
   std::pair<double, double> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::bidirectional_polarization_force(double celli_theta,
                                                                                                             double celli_phi,
                                                                                                             double cellj_theta,
                                                                                                             double cellj_phi)
   {
     // auto prod = pol_dot_product(Xi, p);
     // return prod * unidirectional_polarization_force(Xi, p);

    /*Lets test this theory*/
    // double celli_theta = p_celli->GetCellData()->GetItem("Theta");
    // double celli_phi = p_celli->GetCellData()->GetItem("Phi");
    //
    // double cellj_theta = p_cellj->GetCellData()->GetItem("Theta");
    // double cellj_phi = p_cellj->GetCellData()->GetItem("Phi");

    double prod = pol_dot_product(celli_theta, celli_phi, cellj_theta, cellj_phi);

    std::pair<double, double> unidir_force = unidirectional_polarization_force(celli_theta, celli_phi, cellj_theta, cellj_phi);

    std::pair<double, double> pol_test = std::make_pair(prod*unidir_force.first, prod*unidir_force.second);

    return pol_test;
   }

   /*
   // Euler solver as a functio to call
   // within the force
   */
   template<unsigned SPACE_DIM>
   std::pair<double, double> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::EulerSolver_PolarisationAngles(double celli_theta,
                                                                                                             double celli_phi,
                                                                                                             double neighbour_theta,
                                                                                                             double neighbour_phi)
    {
      EulerIvpOdeSolver euler_solver;
      PolarisationOdeSystem ode_system(celli_theta,celli_phi,neighbour_theta,neighbour_phi);
      std::vector<double> initial_conditions = ode_system.GetInitialConditions();

      // double time_step = 0.25;
      double time_step = SimulationTime::Instance()->GetTimeStep();
      // double time = Timer::GetElapsedTime();
      double time = SimulationTime::Instance()->GetTime();

      std::vector<double> current_y_values=initial_conditions;
      std::vector<double> next_y_values(initial_conditions.size());

      next_y_values = current_y_values;

      OdeSolution solutions;

      // Ensure the internal data structures are the right size
      solutions = euler_solver.Solve(&ode_system, current_y_values, time, time+time_step, time_step, time_step); //1e-5);
       // solution = mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);

      int num_timesteps = solutions.GetNumberOfTimeSteps();

      double new_Theta = solutions.rGetSolutions()[num_timesteps-1][0];
      double new_Phi = solutions.rGetSolutions()[num_timesteps-1][1];

      std::pair<double, double> pol_angles_updated = std::make_pair(new_Theta, new_Phi);

      return pol_angles_updated;

    }




  /*
   * Resistance to bending from the potential U_Epi = Σ(p_i . r_ij/r)^2/2.
   */
   template<unsigned SPACE_DIM>
   c_vector<double, SPACE_DIM> EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::bending_force(CellPtr p_celli,
                                                                                            double celli_theta,
                                                                                            double celli_phi,
                                                                                            double cellj_theta,
                                                                                            double cellj_phi,
                                                                                            c_vector<double, SPACE_DIM> unit_vector,
                                                                                            double dij)
   {

      // auto pi = pol_to_float3(Xi);
     c_vector<double, 3> pi = pol_to_cvector(celli_theta,celli_phi);

     //auto pref_angle = M_PI * (3 - sqrt(5)); //golden_angle
     // auto prodi = (pi.x * r.x + pi.y * r.y + pi.z * r.z) / dist;
     auto prodi = (pi[0] * unit_vector[0] + pi[1] * unit_vector[1] + pi[2] * unit_vector[2]) / dij; //+ cosf(pref_angle);

     // auto r_hat = pt_to_pol(r, dist);
     auto r_hat = pt_to_pol(unit_vector,dij);
     //r_hat.first = theta
     //r_hat.second =  phi

     // dF.x = -prodi / dist * pi.x + powf(prodi, 2) / powf(dist, 2) * r.x;
     // dF.y = -prodi / dist * pi.y + powf(prodi, 2) / powf(dist, 2) * r.y;
     // dF.z = -prodi / dist * pi.z + powf(prodi, 2) / powf(dist, 2) * r.z;
     c_vector<double, 3> dF;
     dF[0] = -prodi / dij * pi[0] + pow(prodi, 2) / pow(dij, 2) * unit_vector[0];
     dF[1] = -prodi / dij * pi[1] + pow(prodi, 2) / pow(dij, 2) * unit_vector[1];
     dF[2] = -prodi / dij * pi[2] + pow(prodi, 2) / pow(dij, 2) * unit_vector[2];

     // // Contribution from (p_j . r_ji/r)^2/2
     // Polarity Xj{Xi.theta - r.theta, Xi.phi - r.phi};
     double Xj_theta = celli_theta-(celli_theta-cellj_theta);
     double Xj_phi = celli_phi-(celli_phi-cellj_phi);

     // auto pj = pol_to_float3(Xj);
     c_vector<double, 3> pj = pol_to_cvector(Xj_theta,Xj_phi);

     // auto prodj = (pj.x * r.x + pj.y * r.y + pj.z * r.z) / dist;
     auto prodj = (pj[0] * unit_vector[0] + pj[1] * unit_vector[1] + pj[2] * unit_vector[2]) / dij; //- cosf(pref_angle);

     // dF.x += -prodj / dist * pj.x + powf(prodj, 2) / powf(dist, 2) * r.x;
     // dF.y += -prodj / dist * pj.y + powf(prodj, 2) / powf(dist, 2) * r.y;
     // dF.z += -prodj / dist * pj.z + powf(prodj, 2) / powf(dist, 2) * r.z;
     dF[0] += -prodj / dij * pj[0] + pow(prodj, 2) / pow(dij, 2) * unit_vector[0];
     dF[1] += -prodj / dij * pj[1] + pow(prodj, 2) / pow(dij, 2) * unit_vector[1];
     dF[2] += -prodj / dij * pj[2] + pow(prodj, 2) / pow(dij, 2) * unit_vector[2];

     return dF;
   }

   /*
    * Add force contribution
    */
    template<unsigned SPACE_DIM>
    void EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::AddForceContribution(AbstractCellPopulation<SPACE_DIM,SPACE_DIM>& rCellPopulation)//(AbstractCellPopulation<SPACE_DIM>& rCellPopulation)
    {
      // This force class is defined for NodeBasedCellPopulations only
      assert(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation) != nullptr);

      NodeBasedCellPopulation<SPACE_DIM>* p_static_cast_cell_population = static_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
      p_static_cast_cell_population->Update();

      for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
            {
              CellPtr celli = *cell_iter;
              // Get the node index corresponding to this cell
              unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(celli);//(*cell_iter);

              Node<SPACE_DIM>* p_node_i = rCellPopulation.GetNode(node_index);

              if (!p_node_i->IsParticle()) // If node A is not a particle, we can do all of the cell-related things
              {

              //Get current polarisation angles for cell i
              double celli_theta = celli->GetCellData()->GetItem("Theta");
              double celli_phi = celli->GetCellData()->GetItem("Phi");

              // Get the location of this node
              const c_vector<double, SPACE_DIM>& r_node_i_location = p_node_i->rGetLocation();

              // Get the radius of this cell
              double radius_of_cell_i = p_node_i->GetRadius();

              // Get the set of node indices corresponding to this cell's neighbours
              std::set<unsigned> neighbouring_node_indices = p_static_cast_cell_population->GetNeighbouringNodeIndices(node_index);
              // std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(node_index);

              // Loop over this set
              for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
                   iter != neighbouring_node_indices.end();
                   ++iter)
                   {
                     CellPtr cellj = rCellPopulation.GetCellUsingLocationIndex(*iter);

                     Node<SPACE_DIM>* p_node_j = rCellPopulation.GetNode(*iter);

                     if (!p_node_j->IsParticle()) // If node A is not a particle, we can do all of the cell-related things
                     {

                     //Get current polarisation angles for cell j
                     double cellj_theta = cellj->GetCellData()->GetItem("Theta");
                     double cellj_phi = cellj->GetCellData()->GetItem("Phi");

                     // Get the location of this node
                     const c_vector<double, SPACE_DIM>& r_node_j_location = p_node_j->rGetLocation();

                     // Get the radius of this cell
                     double radius_of_cell_j = p_node_j->GetRadius();

                     // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
                       c_vector<double, SPACE_DIM> unit_vector = r_node_i_location - r_node_j_location;

                     // Calculate the distance between the two nodes
                     double dij = norm_2(unit_vector);

                     c_vector<double, 3> force;

                     if (dij < radius_of_cell_i + radius_of_cell_j) //Which is the same as dij < 1
                                                                    //(dij <= r_max) if (dist > r_max) return dF;
                     {
                       double pol_dir_ij = pol_dot_product(celli_theta, celli_phi, cellj_theta, cellj_phi);

                       if (pol_dir_ij>0)
                       {

                        auto F = fmax(0.7 - dij, 0) * 2 - fmax(dij - 0.8, 0)/2;
                        // auto F = fmax(0.8 - dij, 0) * 1.0 - fmax(dij - 0.8, 0) * 1.0;


                         force[0] = unit_vector[0] * F / dij;
                         force[1] = unit_vector[1] * F / dij;
                         force[2] = unit_vector[2] * F / dij;

                         force += bending_force(celli,celli_theta,celli_phi,cellj_theta,cellj_phi,unit_vector,dij)*0.3;//*0.3;

                       }
                       else
                       {
                         c_vector<double, 3> pol_vec = pol_to_cvector(celli_theta, celli_phi);

                         auto F = 1; //fmax(0.7 - dij, 0)/2 - fmax(dij - 0.8, 0)*2;
                         // auto F = fmax(0.8 - dij, 0) * 1.0 - fmax(dij - 0.8, 0) * 1.0;

                          force[0] = pol_vec[0] * F / dij;
                          force[1] = pol_vec[1] * F / dij;
                          force[2] = pol_vec[2] * F / dij;

                       }

                     }

                     else
                     {
                       force[0] = 0;
                       force[1] = 0;
                       force[2] = 0;
                     }

                        // p_node_i->ClearAppliedForce();
                        p_node_i->AddAppliedForceContribution(force);
                    }
                  }
                }
            }

    }

template<unsigned SPACE_DIM>
void EpithelialLayerPolarisationForce_wMatrigel<SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
    // Call method on direct parent class
    AbstractForce<SPACE_DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class EpithelialLayerPolarisationForce_wMatrigel<1>;
template class EpithelialLayerPolarisationForce_wMatrigel<2>;
template class EpithelialLayerPolarisationForce_wMatrigel<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EpithelialLayerPolarisationForce_wMatrigel)
