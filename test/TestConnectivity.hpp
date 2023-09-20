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

class TestConnectivityPolareised: public AbstractCellBasedTestSuite
{
public:
  

//////////////////////////////   NEW TEST /////////////////////////////////////

    void Test2dBilayerEquidistantLayers()
    {
        //geometry setup in polar coordinates
        static double rad = 0.6;
        static double Bas_rad = 4;
        static double epsilon = 1*rad;
        static double lum_rad = Bas_rad-epsilon;
        static double basal_node_dist = 1*rad;
        static double dtheta = 2*asin(basal_node_dist/(2*Bas_rad));
        static double dphi = 2*asin(basal_node_dist/(2*(Bas_rad-epsilon)));

        static double dtheta_turns = round(2*M_PI/dtheta);
        static double dphi_turns = round(2*M_PI/dphi);

        static double cut_off_length = 5; //this value was taken from the comparing paper
 

         std::vector<Node<2>*> nodes;

        for (unsigned z = 0; z<dphi_turns; z++){

            Node<2>* p_node(new Node<2>(nodes.size(), false, lum_rad*cos((2*M_PI)*(z/dphi_turns) ), lum_rad*sin((2*M_PI)*(z/dphi_turns) )));
            nodes.push_back(p_node);

        }

        for (unsigned z = 0; z<dtheta_turns; z++){

            Node<2>* p_node(new Node<2>(nodes.size(), true, Bas_rad*cos((2*M_PI)*(z/dtheta_turns)), Bas_rad*sin((2*M_PI)*(z/dtheta_turns))));
            nodes.push_back(p_node);

        }

        NodesOnlyMesh<2> mesh; //make the mesh object to pass the nodes into
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length); //pass the nodes into the mess and defining the node connectivity radius

        //initialise cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_Lpos_type);
            double birth_time = RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        for (unsigned i = 0;i<mesh.GetNumNodes(); i++){
            Node<2> *node = cell_population.GetNode(i);
            node->SetRadius(rad);
        }


        //next we removed all unwanted nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
        c_vector<double,2> centre = cell_population.GetCentroidOfCellPopulation();

        double r = norm_2(location-centre);

            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_rad)
                {   
                        cell_iter->SetCellProliferativeType(p_ME_type);
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
                        cell_iter->SetSrnModel(p_srn_model);
                          

                }
             if (r <= lum_rad+1e-5)
                {
                     
                        cell_iter->SetCellProliferativeType(p_Lpos_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(5);
                        initial_conditions.push_back(5);
                        initial_conditions.push_back(0.1);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                       
                     
                     
                 }  

       }

        std::cout<< cell_population.GetNumNodes() << endl;

        //create the simulation 
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Periodic 3/4 connect Bi - polar IC");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(50);
        
        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        //this has been made homogeneous to investiage the SRN model - cell don't move to isolate connectivity
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); 
        p_differential_adhesion_force->SetMeinekeSpringStiffness(0.0000);
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); 
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); 
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.4); 
        p_differential_adhesion_force->SetCutOffLength(cut_off_length); 
        simulator.AddForce(p_differential_adhesion_force);

        simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }


    }




    void TestSingleCell(){

         // Two cells close to each other

        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(1, false, 0.0, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 5);


        // Establish a CCM for the cells and randomise the birth times
        std::vector<CellPtr> cells;
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {   

            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);

            BscSrnModel* p_srn_model = new BscSrnModel();
            std::vector<double> initial_conditions;
            initial_conditions.push_back(1);
            initial_conditions.push_back(1);
            initial_conditions.push_back(1);
            initial_conditions.push_back(9);
            initial_conditions.push_back(0.5);
            initial_conditions.push_back(0.5);
            initial_conditions.push_back(0.5);
            initial_conditions.push_back(5);
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_BSC_type);
            p_cell->SetSrnModel(p_srn_model);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
            
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Single ODE test");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(0.001);
        
       
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population)); //depends on cell type - set to only 3 neighbours??
        simulator.AddCellKiller(p_killer);
        

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(20); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.4); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(1); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);


       simulator.Solve();





    }

};