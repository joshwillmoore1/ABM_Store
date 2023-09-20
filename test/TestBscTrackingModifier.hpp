#ifndef TESTBSCMODIFIER_HPP_
#define TESTBSCMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "BscSrnModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "BscTrackingModifier.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellIdWriter.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBscModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestUpdateAtEndOfTimeStepNodeBased()
    {
        EXIT_IF_PARALLEL;

        // Create a small 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // ASsociate each cell with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Initial condition for delta, notch, mean_delta
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
            p_cc_model->SetDimension(2);

            BscSrnModel* p_srn_model = new BscSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Create and configure cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchNodeBasedUpdateAtEndOfTimeStep");
        simulator.SetEndTime(0.01);

        // No mechanics so cells don't move

        // Create a Delta-Notch tracking modifier and add it to the simulation
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run simulation
        simulator.Solve();

        // Test that the node velocities file exists
        OutputFileHandler output_file_handler("TestBscNodeBasedUpdateAtEndOfTimeStep", false);
        FileFinder generated = output_file_handler.FindFile("results_from_time_0/nodevelocities.dat");
        TS_ASSERT(generated.Exists());
    }

    void TestHeterogeneousBscOnUntetheredTwoCellSystem()
    {
        EXIT_IF_PARALLEL;

        // Two cells close to each other
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.5, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Establish a CCM for the cells and randomise the birth times
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Cell #1:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.9);

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
        p_cc_model->SetDimension(2);

        BscSrnModel* p_srn_model = new BscSrnModel();
        p_srn_model->SetInitialConditions(starter_conditions);
        CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(0.0);
        cells.push_back(p_cell);

        // Cell #2:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions_2;
        starter_conditions_2.push_back(0.19);
        starter_conditions_2.push_back(0.5);
        starter_conditions_2.push_back(0.19);
        starter_conditions_2.push_back(0.5);
        starter_conditions_2.push_back(0.19);
        starter_conditions_2.push_back(0.5);
        starter_conditions_2.push_back(0.19);
        

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model_2 = new UniformCellCycleModel();
        p_cc_model_2->SetDimension(2);

        BscSrnModel* p_srn_model_2 = new BscSrnModel();
        p_srn_model_2->SetInitialConditions(starter_conditions_2);
        CellPtr p_cell_2(new Cell(p_state, p_cc_model_2, p_srn_model_2));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->SetBirthTime(0.0);
        cells.push_back(p_cell_2);

        // Create the cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Set up the simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestBscTwoCell_heterogee");
        simulator.SetEndTime(10.0);

        // No mechanics so cells don't move

        // Add Bsc tracking modifier
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run the simulation
        simulator.Solve();

        // Acquire cell pointers
        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1b = cell_population.GetCellUsingLocationIndex(1);



        // Now move cell so the cells have no neighboursthen they both run to a homogeneous steady state (note here mean delta=0)
        c_vector<double,2> old_point;
        old_point = static_cast<NodesOnlyMesh<2>* >(&(simulator.rGetCellPopulation().rGetMesh()))->GetNode(1)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+2.0;
        new_point.rGetLocation()[1] = old_point[1];
        static_cast<NodesOnlyMesh<2>* >(&(simulator.rGetCellPopulation().rGetMesh()))->SetNode(1, new_point);

        // Run the simulation
        simulator.SetEndTime(20.0);
        simulator.Solve();

        // Avoid memory leaks
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestHomogeneousBscOnUntetheredTwoCellSystem()
    {
        EXIT_IF_PARALLEL;

        // Two cells close to each other
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.5, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Establish a CCM for the cells and randomise the birth times
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Cell #1:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.9);

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
        p_cc_model->SetDimension(2);

        BscSrnModel* p_srn_model = new BscSrnModel();
        p_srn_model->SetInitialConditions(starter_conditions);
        CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(0.0);
        cells.push_back(p_cell);

        // Cell #2:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions_2;
        starter_conditions_2.push_back(0.9);
        starter_conditions_2.push_back(0.5);
        starter_conditions_2.push_back(0.9);
        starter_conditions_2.push_back(0.5);
        starter_conditions_2.push_back(0.9);
        starter_conditions_2.push_back(0.5);
        starter_conditions_2.push_back(0.9);
        

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model_2 = new UniformCellCycleModel();
        p_cc_model_2->SetDimension(2);

        BscSrnModel* p_srn_model_2 = new BscSrnModel();
        p_srn_model_2->SetInitialConditions(starter_conditions_2);
        CellPtr p_cell_2(new Cell(p_state, p_cc_model_2, p_srn_model_2));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->SetBirthTime(0.0);
        cells.push_back(p_cell_2);


        // Create the cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Set up the simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestBscTwoCell_homgee");
        simulator.SetEndTime(10.0);

        // Add bsc tracking modifier
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Define the radius of interaction as we're dealing with a node-based simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(0.75);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Acquire cell pointers
        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1b = cell_population.GetCellUsingLocationIndex(1);

        // Avoid memory leaks
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    
};

#endif /*TESTBSCMODIFIER_HPP_*/
