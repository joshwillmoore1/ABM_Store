#ifndef TESTBSCSRNMODELS_HPP_
#define TESTBSCSRNMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractSrnModel.hpp"
#include "NullSrnModel.hpp"
#include "BscSrnModel.hpp"
#include "Goldbeter1991SrnModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestBscSrnModels : public AbstractCellBasedTestSuite
{
public:

    void TestBscSrnCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(BscSrnModel srn_model);

        BscSrnModel* p_srn_model = new BscSrnModel();

        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        p_srn_model->SetInitialConditions(starter_conditions);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->GetCellData()->SetItem("mean delta", 1.0);
        p_cell->GetCellData()->SetItem("mean dkk1", 1.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        // Now updated to initial conditions
        TS_ASSERT_DELTA(p_srn_model->GetNotch(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetDelta(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetHes1(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetBcat(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetTCF(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetDkk1(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetNrg1(), 0.5, 1e-4);

        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        double end_time = 10.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_srn_model->SimulateToCurrentTime();
        }
    }

    void TestBscSrnCreateCopy()
    {
        // Test with DeltaNotchSrnModel
        BscSrnModel* p_model= new BscSrnModel;

        // Set ODE system
        std::vector<double> state_variables;
        state_variables.push_back(2.0);
        state_variables.push_back(3.0);
        state_variables.push_back(2.0);
        state_variables.push_back(3.0);
        state_variables.push_back(2.0);
        state_variables.push_back(3.0);
        state_variables.push_back(2.0);
        p_model->SetOdeSystem(new BscOdeSystem(state_variables));

        p_model->SetInitialConditions(state_variables);

        // Create a copy
        BscSrnModel* p_model2 = static_cast<BscSrnModel*> (p_model->CreateSrnModel());

        // Check correct initializations
        TS_ASSERT_EQUALS(p_model2->GetNotch(), 2.0);
        TS_ASSERT_EQUALS(p_model2->GetDelta(), 3.0);

        // Destroy models
        delete p_model;
        delete p_model2;
    }

    void TestArchiveBscSrnModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "bsc_srn.arch";

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractSrnModel* p_srn_model = new BscSrnModel;

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            // We must create a cell to be able to initialise the cell srn model's ODE system
            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->GetCellData()->SetItem("mean delta", 10.0);
            p_cell->GetCellData()->SetItem("mean dkk1", 10.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_cell->SetBirthTime(0.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Read mean Delta from CellData
            static_cast<BscSrnModel*>(p_srn_model)->UpdateBsc();
            TS_ASSERT_DELTA(static_cast<BscSrnModel*>(p_srn_model)->GetMeanNeighbouringDelta(), 10.0, 1e-12);
            TS_ASSERT_DELTA(static_cast<BscSrnModel*>(p_srn_model)->GetMeanNeighbouringDkk1(), 10.0, 1e-12);

            output_arch << p_srn_model;

            // Note that here, deletion of the cell-cycle model and srn is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractSrnModel* p_srn_model;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_srn_model;

            TS_ASSERT_DELTA(static_cast<BscSrnModel*>(p_srn_model)->GetMeanNeighbouringDelta(), 10.0, 1e-12);
            TS_ASSERT_DELTA(static_cast<BscSrnModel*>(p_srn_model)->GetMeanNeighbouringDkk1(), 10.0, 1e-12);

            delete p_srn_model;
        }
    }



    
};

#endif /* TESTBSCSRNMODELS_HPP_ */
