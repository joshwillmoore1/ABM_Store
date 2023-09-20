#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "MamDeltaNotchSrnModel.hpp"
#include "MamDeltaNotchTrackingModifier.hpp"
#include "MamNDMAdaptiveTrackingModifier.hpp"
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
#include "RandomMotionForce.hpp"
#include "CellLabel.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "BiBuckleRandomMotionForce.hpp"


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


class TestBuckleNDM: public AbstractCellBasedTestSuite
{
public:


void Test2dBilayerBucklePressLumenForPolarisation()
    {


       // RandomNumberGenerator::Instance()->Reseed(5);

        double cell_radius = 0.5; // this implies the spring rest length is 1 
        static double ORGANOID_DOMAIN = 14; //31
        static double ORGANOID_RADIUS = (ORGANOID_DOMAIN-1)/2 -1; //14
 
        std::vector<Node<2>*> nodes;

        //first we make a cubic lattice of nodes

        for (unsigned i=0; i<ORGANOID_DOMAIN; i++)
            {
                for(unsigned j=0; j<ORGANOID_DOMAIN; j++ )
                {
                        Node<2>* p_node(new Node<2>(nodes.size(), false, i, j ));
                        p_node->SetRadius(cell_radius);
                        nodes.push_back(p_node);
                    
                }
            }

        NodesOnlyMesh<2> mesh; 
        double cut_off_length = 2; //set cut off length to be 4 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length);

        //initialise cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_Lpos_type);
            double birth_time = RandomNumberGenerator::Instance()->ranf()*2;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        

      
        //next we removed all unwanted nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
        //c_vector<double,2> centre = cell_population.GetCentroidOfCellPopulation();

        double r = pow(location[0]-(ORGANOID_DOMAIN-1)/2,2) + pow(location[1]-(ORGANOID_DOMAIN-1)/2,2);

            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > pow(ORGANOID_RADIUS-1,2))
                {   
                    cell_iter->SetCellProliferativeType(p_ME_type);
                        MamDeltaNotchSrnModel* p_srn_model = new MamDeltaNotchSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.1);
                        initial_conditions.push_back(0.1);
                        initial_conditions.push_back(4);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                     
                          

                }
             if (r <= pow(ORGANOID_RADIUS-1,2) && r > pow(ORGANOID_RADIUS-2,2))
                {           
                    cell_iter->SetCellProliferativeType(p_Lpos_type);
                        MamDeltaNotchSrnModel* p_srn_model = new MamDeltaNotchSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(4);
                        initial_conditions.push_back(4);
                        initial_conditions.push_back(0.1);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);

                     
                
         
                 }  
                

             if (r <= pow(ORGANOID_RADIUS-2,2) ){
  
                    cell_iter->SetCellProliferativeType(p_diff_type);
                    boost::shared_ptr<AbstractCellProperty> pLabel(CellPropertyRegistry::Instance()->Get<CellLabel>());
                    cell_iter->AddCellProperty(pLabel);

                }

                if (r > pow(ORGANOID_RADIUS,2))
                {   
                        cell_iter->Kill();
                 
                }


       }

       cell_population.RemoveDeadCells();

        std::cout<< cell_population.GetNumNodes() << endl;

        //cell type writer:
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        //create the simulation 
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test Srn with Buckle");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetDt(0.002);
        simulator.SetEndTime(100);
        
        //this modifier is how the SRN network interactions with the cells

        //MAKE_PTR(MamNDMAdaptiveTrackingModifier<2>, p_modifier); // for adaptive signal strength
        MAKE_PTR(MamDeltaNotchTrackingModifier<2>, p_modifier); // for fixed signal strength
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(BiBuckleRandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.0001); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        //this has been made homogeneous to investiage the SRN model - cell don't move to isolate connectivity
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); 
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(cell_radius/2); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(2*cell_radius); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);

       simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }








    }

};