#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "BscSrnModel.hpp"
#include "BscTrackingModifier.hpp"
#include "OrganoidSpringForce.hpp"
#include "LumenKiller.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "OrganoidG1FixedCellCycleModel.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellData.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "WildTypeCellMutationState.hpp"


//files for periodic geometry
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "Cylindrical2dMesh.hpp"
#include "VoronoiDataWriter.hpp"

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestRadiiConnectSrn: public AbstractCellBasedTestSuite
{
public:
    void TestHexLatticeConnect(){

        //construst a small lattice and check the number of neighbours dependent on cut off length

        static double cut_off_length = 3.5;

        HoneycombMeshGenerator generator(3,3);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, cut_off_length);




        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(3);


        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(3);

            BscSrnModel* p_srn_model = new BscSrnModel();
            std::vector<double> initial_conditions;
            initial_conditions.push_back(0.1);
            initial_conditions.push_back(0.1);
            initial_conditions.push_back(5);
            initial_conditions.push_back(9);
            initial_conditions.push_back(0.5);
            initial_conditions.push_back(0.5);
            initial_conditions.push_back(0.5);
            initial_conditions.push_back(5);
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
            p_cell->SetSrnModel(p_srn_model);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        for (unsigned i = 0;i<mesh.GetNumNodes(); i++){
            Node<2> *node = cell_population.GetNode(i);
            node->SetRadius(1.15);
        }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Hex Connectivity");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(20);

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(0.001); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.3); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(cut_off_length); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);


        simulator.Solve();


    }






};