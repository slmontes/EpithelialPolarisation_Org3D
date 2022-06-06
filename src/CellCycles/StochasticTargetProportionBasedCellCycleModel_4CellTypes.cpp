/*

Copyright (c) 2005-2014, University of Oxford.
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
/*
*By: Sandra M
*For test ran during first half of 2020
*Set proportions test
*/

#include "StochasticTargetProportionBasedCellCycleModel_4CellTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "TACellMutationState.hpp"
#include "EnterocyteCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "Exception.hpp"

StochasticTargetProportionBasedCellCycleModel_4CellTypes::StochasticTargetProportionBasedCellCycleModel_4CellTypes()
    : AbstractSimplePhaseBasedCellCycleModel(),
      mTargetProportion(DOUBLE_UNSET),
      mTargetProportionStifferCells(DOUBLE_UNSET),
	    mCC_Scale(DOUBLE_UNSET)
{
}

StochasticTargetProportionBasedCellCycleModel_4CellTypes::~StochasticTargetProportionBasedCellCycleModel_4CellTypes()
{
}

void StochasticTargetProportionBasedCellCycleModel_4CellTypes::SetTargetProportion(double targetProportion)
{
	mTargetProportion = targetProportion;
}

double StochasticTargetProportionBasedCellCycleModel_4CellTypes::GetTargetProportion()
{
	return mTargetProportion;
}

void StochasticTargetProportionBasedCellCycleModel_4CellTypes::SetTargetProportionStifferCells(double targetProportionStifferCells)
{
	mTargetProportionStifferCells = targetProportionStifferCells;
}

double StochasticTargetProportionBasedCellCycleModel_4CellTypes::GetTargetProportionStifferCells()
{
	return mTargetProportionStifferCells;
}

AbstractCellCycleModel* StochasticTargetProportionBasedCellCycleModel_4CellTypes::CreateCellCycleModel()
{
    // Create a new cell-cycle model
	StochasticTargetProportionBasedCellCycleModel_4CellTypes* p_model = new StochasticTargetProportionBasedCellCycleModel_4CellTypes();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */

    p_model->SetBirthTime(mBirthTime);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetTargetProportion(mTargetProportion);
    p_model->SetTargetProportionStifferCells(mTargetProportionStifferCells);
    p_model->SetCellCycleLengthScale(mCC_Scale);

    return p_model;

}

void StochasticTargetProportionBasedCellCycleModel_4CellTypes::SetG1Duration()
{
    /*
     * Original code from Almet et al
     */
    //  RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    //
    // if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    // {
    //     mG1Duration = GetStemCellG1Duration() + 2*p_gen->ranf(); // U[0,2]
    // }
    // //If we have a transit cell, i.e. soft or hard cell
    // else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    // {
    //     mG1Duration = GetTransitCellG1Duration() + 2*p_gen->ranf(); // U[0,2]
    // }
    // //if we have a non-epithelial cell
    // else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    // {
    //     mG1Duration = DBL_MAX;
    // }
    // else
    // {
    //     NEVER_REACHED;
    // }

    /*
     * Modification made by Sandra M.
     */

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double cc_scaling = GetCellCycleLengthScale();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 2*p_gen->ranf(); // U[0,2]
    }
    //If our transint cell represents a Stem Cell in our simulation
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (!mpCell->GetMutationState()->IsType<TACellMutationState>()) && (!mpCell->GetMutationState()->IsType<PanethCellMutationState>()) && (!mpCell->GetMutationState()->IsType<EnterocyteCellMutationState>()))
    {
      mG1Duration = (GetTransitCellG1Duration() + 2*p_gen->ranf())*(cc_scaling*1); // U[0,2]
    }
    //If our transint cell represents a TA cell
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<TACellMutationState>()))
    {
        mG1Duration = (GetTransitCellG1Duration() + 2*p_gen->ranf())*(cc_scaling*1); // U[0,2]
    }
    //If our transint cell represents a Paneth cell
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<PanethCellMutationState>()))
    {
        // mG1Duration = (GetTransitCellG1Duration() + 2*p_gen->ranf())*(cc_scaling*1); // U[0,2]
        mG1Duration = DBL_MAX;
    }
    //If our transint cell represents a Paneth cell
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<EnterocyteCellMutationState>()))
    {
        // mG1Duration = (GetTransitCellG1Duration() + 2*p_gen->ranf())*(cc_scaling*1); // U[0,2]
        mG1Duration = DBL_MAX;
    }
    //If we have a non-epithelial cell
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void StochasticTargetProportionBasedCellCycleModel_4CellTypes::InitialiseDaughterCell()
{

    /*
     * ORIGINAL FROM ALMET ET AL.
     *
     * We set the daughter cell to have a probability mTargetProportion of being a soft cell or
     * a hard cell
     */
       //
       // //Get the target proportion probability
       // double target_proportion = GetTargetProportion();
       //
       // //Get random number
       // RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
       // double uniform_random_number = p_gen->ranf();
       //
       // //Set daughter cell to be a transit cell
       // boost::shared_ptr<AbstractCellProperty> p_transit_type =
       //     mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
       // mpCell->SetCellProliferativeType(p_transit_type);
       //
       // if (uniform_random_number < target_proportion) // set probability of becoming a hard cell
       // {
       //   boost::shared_ptr<AbstractCellProperty> p_wildtype_state =
       //       mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
       //   mpCell->SetMutationState(p_wildtype_state);
       // }
       // else
       // {
       //   boost::shared_ptr<AbstractCellProperty> p_hard_state =
       //       mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<HardCellMutationState>();
       //   mpCell->SetMutationState(p_hard_state);
       // }
       //
       // AbstractSimpleCellCycleModel::InitialiseDaughterCell();


     // /*
     //  * MODIFIED BY SANDRA M.
     //  *
     //  * We set different proportions for each cell type.
     //  *
     //  */
     //
     //  /*
     //   * WE TRIED TO SET RESTRICTIONS IN THE PRODUCCTION OF CELLS:
     //   * Stem cells can produce: Stem Cells and TA CELLS
     //   * TA (TA) cells can produce: Stem Cells and TA CELLS
     //   * Paneth cells can produce: TA Cells and Paneth cells
     //   * We set different proportions for each cell type.
     //   *
     //   */
     //
     // //Get the target proportion probability
     // double target_proportion = GetTargetProportion();
     //
     // //Get the target proportion probability
     // double target_proportion_stiffer_cells = GetTargetProportionStifferCells();
     //
     // //Get random number
     // RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
     // // RandomNumberGenerator* p_gen2 = RandomNumberGenerator::Instance();
     // double uniform_random_number = p_gen->ranf();
     // double uniform_random_number_forStifferCells = p_gen->ranf();
     //
     //  	//Set daughter cell to be a transit cell
     //  	boost::shared_ptr<AbstractCellProperty> p_transit_type =
     //  			mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
     //  	mpCell->SetCellProliferativeType(p_transit_type);
     //
     //    if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>()))
     //    {
     //
     //    	if (uniform_random_number < target_proportion) // set probability of becoming a stiffer cell --- target_proportion
     //    	   {
     //    		   boost::shared_ptr<AbstractCellProperty> p_wildtype_state =
     //    				  mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
     //           mpCell->SetMutationState(p_wildtype_state);
     //          }
     //    	else
     //    	   {
     //           if (uniform_random_number_forStifferCells <= target_proportion_stiffer_cells)
     //              {
     //                boost::shared_ptr<AbstractCellProperty> p_TA_state =
     //      				      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
     //  		          mpCell->SetMutationState(p_TA_state);
     //              }
     //           else
     //              {
     //                boost::shared_ptr<AbstractCellProperty> p_Paneth_state =
     //      				      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
    	// 	            mpCell->SetMutationState(p_Paneth_state);
     //              }
    	//        }
     //     }
     //     else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<TACellMutationState>()))
     //     {
     //           boost::shared_ptr<AbstractCellProperty> p_TA_state =
     // 				      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
     //           mpCell->SetMutationState(p_TA_state);
     //
     //      }
     //      else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<PanethCellMutationState>()))
     //      {
     //
     //            boost::shared_ptr<AbstractCellProperty> p_Paneth_state =
     //  				      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
  		//           mpCell->SetMutationState(p_Paneth_state);
     //
     //      }

     /*
      * MODIFIED BY SANDRA M.
      *
      * We set different proportions for each cell type.
      *
      */

      /*
       * WE TRIED TO SET RESTRICTIONS IN THE PRODUCCTION OF CELLS:
       * Stem cells can produce: All 3 cell types according to the initial cell
       * proportions 9% Paneth, 42% stem, 49% TA cells. However, we removed the
       * possibility to produce TAs because otherwise it gets overcrowded with it.
       * TA (TA) cells can produce: TA CELLS.
       * Paneth cells can produce its own type but DO NOT reproduce due to its
       * cell cycle length.
       *
       */

     //Get the target proportion probability
     double target_proportion = GetTargetProportion();

     //Get the target proportion probability
     double target_proportion_stiffer_cells = GetTargetProportionStifferCells();

     //Get random number
     RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
     double uniform_random_number = p_gen->ranf();

      	//Set daughter cell to be a transit cell
      	boost::shared_ptr<AbstractCellProperty> p_transit_type =
      			mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
      	mpCell->SetCellProliferativeType(p_transit_type);

        if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (!mpCell->GetMutationState()->IsType<TACellMutationState>()) && (!mpCell->GetMutationState()->IsType<PanethCellMutationState>()) && (!mpCell->GetMutationState()->IsType<EnterocyteCellMutationState>()))
        {

        	if (uniform_random_number < 0.09) // set probability of becoming a Paneth cell --- target_proportion
        	   {
               boost::shared_ptr<AbstractCellProperty> p_Paneth_state =
                   mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
               mpCell->SetMutationState(p_Paneth_state);
             }
          else if (uniform_random_number < 0.98)
            {
              boost::shared_ptr<AbstractCellProperty> p_wildtype_state =
                 mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
              mpCell->SetMutationState(p_wildtype_state);
            }
          else
            {
              boost::shared_ptr<AbstractCellProperty> p_TA_state =
                  mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
              mpCell->SetMutationState(p_TA_state);
            }
         }
         else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<TACellMutationState>()))
         {
           double sim_time = SimulationTime::Instance()->GetTime();
            if (sim_time < 120)
            {
              if (uniform_random_number < 0.30) // set probability of becoming a Paneth cell --- target_proportion
                {
                   boost::shared_ptr<AbstractCellProperty> p_EC_state =
                       mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EnterocyteCellMutationState>();
                   mpCell->SetMutationState(p_EC_state);
                 }
              else
                {
                  boost::shared_ptr<AbstractCellProperty> p_TA_state =
                      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
                  mpCell->SetMutationState(p_TA_state);
                }
              }
            else
            {
              if (uniform_random_number > 0.30) // set probability of becoming a Paneth cell --- target_proportion
                {
                   boost::shared_ptr<AbstractCellProperty> p_EC_state =
                       mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EnterocyteCellMutationState>();
                   mpCell->SetMutationState(p_EC_state);
                 }
              else
                {
                  boost::shared_ptr<AbstractCellProperty> p_TA_state =
                      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
                  mpCell->SetMutationState(p_TA_state);
                }
              }

          }
          else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<PanethCellMutationState>()))
          {

                boost::shared_ptr<AbstractCellProperty> p_Paneth_state =
      				      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
  		          mpCell->SetMutationState(p_Paneth_state);

          }

      // RandomNumberGenerator::Destroy();
      AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();

}

void StochasticTargetProportionBasedCellCycleModel_4CellTypes::SetCellCycleLengthScale(double cc_scale){

	mCC_Scale = cc_scale;

}

double StochasticTargetProportionBasedCellCycleModel_4CellTypes::GetCellCycleLengthScale(){

	return mCC_Scale;

}


void StochasticTargetProportionBasedCellCycleModel_4CellTypes::OutputCellCycleModelParameters(out_stream& rParamsFile)
{

    *rParamsFile << "\t\t\t<TargetProportion>" << mTargetProportion << "</TargetProportion>\n";
    *rParamsFile << "\t\t\t<TargetProportionStifferCells>" << mTargetProportionStifferCells << "</TargetProportionStifferCells>\n";
    *rParamsFile << "\t\t\t<CellCycleLengthScale>" << mCC_Scale << "</CellCycleLengthScale>\n";

    // Nothing to output, so just call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticTargetProportionBasedCellCycleModel_4CellTypes)
