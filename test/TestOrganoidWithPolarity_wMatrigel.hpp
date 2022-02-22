/*
 * Author: Sandra Montes
 * Date: 08/07/2021
 *
 */

#ifndef TESTORGANOIDWITHPOLARITY_HPP_
#define TESTORGANOIDWITHPOLARITY_HPP_

/*
 * == Test a node-based simulation using EpithelialLayerPolarisationForce ==
 *
 * EMPTYLINE
 *
 * As in previous tutorials, we begin by including the necessary header files. We have
 * encountered these files already. Recall that often, either {{{CheckpointArchiveTypes.hpp}}}
 * or {{{CellBasedSimulationArchiver.hpp}}} must be included the first Chaste header.
 */
#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>
#include <cmath>
#include <iostream>

#include "CellBasedSimulationArchiver.hpp" //Must be included before any other cell_based or crypt headers

#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method
#include "SmartPointers.hpp" //Enables macros to save typing

/*
 * Cell-based population and simulation headers
*/
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "NodeBasedCellPopulation.hpp" // Centre-based cell population
#include "CellsGenerator.hpp" //Generates a vector of cells for a given mesh, according cell-cycle model and spatial dimension
#include "FixedG1GenerationalCellCycleModel.hpp" //Fixed-G1phase-based cell cycle model
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"  //Simulates the evolution of the population
#include "StochasticTargetProportionBasedCellCycleModelSandra.hpp"

/*
 * Cell types headers
*/
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PanethCellMutationState.hpp" //Mutation class that defines Paneth cells
#include "TACellMutationState.hpp" //Mutation class that defines TA cells

/*
 * Cell writers
*/
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

/*
 * The next header defines the simulation class modifier that includes the functionality
 * for solving each cell's polarisation vector ODE system at each time step,
 * using information about neighbouring cells through the {{{CellData}}} class.
 * This modifier leads to the {{{CellData}}} cell property being updated at each timestep.
 */
#include "PolarisationTrackingModifier_bending.hpp"
#include "PolarisationTrackingModifier_bending_wWrite.hpp"

/*
 * Forces inside the simulation
 */
#include "EpithelialLayerPolarisationForce.hpp" //Force affected by the polarisation vector of each cell.
#include "EpithelialLayerLinearSpringForce.hpp"
#include "EpithelialLayerPolarisationForce_wMatrigel.hpp"

/*
 * Cell killer
 */
#include "LumenCellKiller.hpp"

//Included for tests (may not be needed anymore)
// #include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"  //to use PRINT_VARIABLE
#include "FakePetscSetup.hpp" //Forbids tests running in parallel


/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestOrganoidWithPolarity : public AbstractCellBasedTestSuite  //CHANE BACK TO THIS FOR ALL OTHER TESTS!!!!!!!!
{

public:

     /*
      * == Test a node-based simulation using EpithelialLayerPolarisationForce ==
      */
      void TestNodeBasedOrganoidWithPolarisation() throw(Exception)
      {
          /* We include the next line because the test is not
           *  yet implemented in parallel. */
          EXIT_IF_PARALLEL;

          /*
           * Example Organoid
           */
           std::vector<Node<3>*> nodes;

           double totalPoints = 31;

           // double dist_to_nb = 0.75;
           // // auto r_max = pow((totalPoints) / 0.64, 1. / 3) * dist_to_nb / 2;
           // auto r_max = pow((totalPoints) / 0.9, 1. / 3) * dist_to_nb / 2;

           std::vector<double> linspaced;

           double startj = 1 - 1.0 / totalPoints;
           double endj = 1.0 / totalPoints - 1;
           double delta = (endj - startj) / (totalPoints - 1);

           for(int j=0; j < totalPoints-1; ++j)
             {
               linspaced.push_back(startj + delta * j);
             }
           linspaced.push_back(endj);

             for (auto i=0; i<totalPoints; i++){
               // auto r = r_max * pow(rand() / (RAND_MAX + 1.), 1. / 3);
               // auto theta = acos(2. * rand() / (RAND_MAX + 1.) - 1);
               // auto phi = rand() / (RAND_MAX + 1.) * 2 * M_PI;
               // points.h_X[i].x = r * sin(theta) * cos(phi);
               // points.h_X[i].y = r * sin(theta) * sin(phi);
               // points.h_X[i].z = r * cos(theta);

               // double r_max = 0.5;

                double golden_angle = M_PI * (3 - sqrt(5));
                double theta = golden_angle * i;
                std::vector<double> z = linspaced;
                double r = sqrt(1 - z[i] * z[i]);
                // double theta = acosf(2. * rand() / (RAND_MAX + 1.) - 1);
                // double phi = rand() / (RAND_MAX + 1.) * 2 * M_PI;
                double Xx = r * cos(theta); //* cos(phi);
                double Xy = r * sin(theta); //* sin(phi);
                double Xz = z[i]; //r * cos(theta); //
                unsigned indx = i;
                nodes.push_back(new Node<3>(indx,  false,  Xx, Xy, Xz));

               // fstream myfile;
               // myfile.open ("Sphere_InitialConditions.txt", std::ios_base::app | std::ios_base::in);
               // if (myfile.is_open())
               // myfile << "NODE: " << i << std::endl;
               // myfile << "r: " << r << std::endl;
               // myfile << "theta: " << theta << std::endl;
               // // myfile << "phi: " << phi << std::endl;
               // myfile.close();
             }

           NodesOnlyMesh<3> mesh;
           mesh.ConstructNodesWithoutMesh(nodes, 1.0); // Distance cut-off usually 1.5

           std::vector<CellPtr> cells;
           // MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);   //%%%%%%%THIS IS TO PREVENT CELL DIVISION!!!!

           boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
           boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
           boost::shared_ptr<AbstractCellProperty> p_paneth_state(CellPropertyRegistry::Instance()->Get<PanethCellMutationState>());
           boost::shared_ptr<AbstractCellProperty> p_TA_state(CellPropertyRegistry::Instance()->Get<TACellMutationState>());
           boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

          // CellsGenerator<StochasticTargetProportionBasedCellCycleModelSandra, 3> cells_generator;
          // cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);

          double target_proportion = 0.8; // + 0.05*tpropn;
          double target_proportion_for_mutants = 0.9;
          double cc_scale = 1.0;

          double init_stem_proportion = 0.42;
          double init_paneth_proportion = 0.09;
          double init_TA_proportion = 0.49;

           std::vector<unsigned> num_paneth;
           std::vector<unsigned> num_TA;
           std::vector<unsigned> num_SC;

          /* Now we need to associate each cell with a {{{DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel}}}
           * and initialise it (each instance is different). So, for each cell we
           * initialise its cell cycle model and randomly set its birthtime.
           */
           // NodeBasedCellPopulation<3> cell_population(mesh, cells);

           for (auto i=0; i<mesh.GetNumNodes(); i++)
          {
            //Set cell cycle
             StochasticTargetProportionBasedCellCycleModelSandra* p_cycle_model = new StochasticTargetProportionBasedCellCycleModelSandra();
             p_cycle_model->SetCellCycleLengthScale(cc_scale);

             //To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
             // ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
             double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf();
             p_cycle_model->SetBirthTime(-birth_time);

             CellPtr p_cell(new Cell(p_state, p_cycle_model));
       			 // p_cell->SetCellCycleModel(p_cycle_model);
             p_cell->SetCellProliferativeType(p_stem_type); //Set the cell to be a differentiated cell
             p_cell->InitialiseCellCycleModel();

             cells.push_back(p_cell);
          }

           // CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
           // cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

           NodeBasedCellPopulation<3> cell_population(mesh, cells);

           for (auto i=0; i<nodes.size(); i++){
             CellPtr p_cell_ib = cell_population.GetCellUsingLocationIndex(i);
             // p_cell_ib->SetCellProliferativeType(p_diff_type);

             //Randomly generate number
             double random_number = RandomNumberGenerator::Instance()->ranf();

             if(random_number < init_paneth_proportion) //Assign cells to be Paneth with 1 - target_proportion
              {
                p_cell_ib->SetMutationState(p_paneth_state); //Set the cell to be paneth cell
                num_paneth.push_back(1);
              }
              else if (random_number < (init_paneth_proportion+init_stem_proportion))
              {
                p_cell_ib->SetCellProliferativeType(p_stem_type); //Set the cell to be stem cell
                num_SC.push_back(1);
              }
              else
              {
                p_cell_ib->SetMutationState(p_TA_state); //Set the cell to be TA cell
                num_TA.push_back(1);
              }

            c_vector<double, 3> cell_location = cell_population.GetLocationOfCellCentre(p_cell_ib);

            double dist = norm_2(cell_location);

            double init_Theta = acosf(cell_location[2] / dist);// + rand() / (RAND_MAX + 1.) * 0.5;
            double init_Phi = atan2(cell_location[1], cell_location[0]);// + rand() / (RAND_MAX + 1.) * 0.5;

            p_cell_ib->GetCellData()->SetItem("Theta",init_Theta);
            p_cell_ib->GetCellData()->SetItem("Phi",init_Phi);

            fstream myfile;
            myfile.open ("Sphere_InitialConditions.txt", std::ios_base::app | std::ios_base::in);
            if (myfile.is_open())
            myfile << "NODE: " << i << std::endl;
            myfile << "theta: " << init_Theta << std::endl;
            myfile << "phi: " << init_Phi << std::endl;
            myfile.close();
           }

           std::cout << "num_paneth: " << num_paneth.size() << std::endl;
           std::cout << "num_SC: " << num_SC.size() << std::endl;
           std::cout << "num_TA: " << num_TA.size() << std::endl;

          cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
          cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
          cell_population.AddCellWriter<CellIdWriter>();
          cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
          cell_population.AddCellWriter<CellAgesWriter>();

          OffLatticeSimulation<3> simulator(cell_population);
          simulator.SetOutputDirectory("TestOrganoidWithPolMatrigel_Organoid_PCStiff_4_SCandTAStiff_1_bendingForceTimes0.3_100hr_x10");

          double dt = 0.005; //Set dt
          double sampling_timestep = 0.5/dt;//4800; //Set sampling timestep
      		double end_time = 110;//168; //Set end time (168hrs is 7days) ... This can be changed

          simulator.SetDt(dt); //Set the timestep dt for force volution
          simulator.SetSamplingTimestepMultiple(sampling_timestep); //Set the sampling timestep multiple for animations
          simulator.SetEndTime(end_time); //Set the number of hours to run the simulation to

          /* As we are using a node-based cell population, we use an appropriate force law. */
          MAKE_PTR(EpithelialLayerPolarisationForce_wMatrigel<3>, p_force);
          simulator.AddForce(p_force);

          MAKE_PTR(EpithelialLayerLinearSpringForce<3>, p_spring_force);
          p_spring_force->SetCutOffLength(1.0);

          //Set all the spring stiffness variables
      		double epithelial_epithelial_stiffness = 1.0; //Epithelial-epithelial spring connections
      		double epithelial_nonepithelial_stiffness = 15.0; //Epithelial-non-epithelial spring connections
      		double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections
          double stiffness_ratio_paneth = 4.0;
          double stiffness_ratio_TA = 1.0;
          //Set the spring stiffnesses
          p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
          p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
          p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
          p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio_paneth);
          p_spring_force->SetTACellStiffnessRatio(stiffness_ratio_TA);
          simulator.AddForce(p_spring_force);

          // Add tracking modifier - we don't have it at this point
          MAKE_PTR(PolarisationTrackingModifier_bending<3>, p_modifier);
          simulator.AddSimulationModifier(p_modifier);

          /* Add a cell killer. */
          // MAKE_PTR_ARGS(LumenCellKiller, p_lumen_killer, (&cell_population));
          // simulator.AddCellKiller(p_lumen_killer);

          // MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force2);
          // simulator.AddForce(p_force2);

          simulator.Solve();
      }

};

#endif /* TESTORGANOIDWITHPOLARITY_HPP_ */
