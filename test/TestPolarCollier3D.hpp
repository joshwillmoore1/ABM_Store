#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CollierSrnModel.hpp"
#include "CollierTrackingModifier.hpp"
#include "CollierFixedTrackingModifier.hpp"
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
#include "BiBuckleRandomMotionForce.hpp"
#include <vector>

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"





class TestPolarCollier3D: public AbstractCellBasedTestSuite
{
public:


const double spring_const = 25;
const double Random_pert = 0.00025;

const double basal_Delta_IC = 0.2;
const double luminal_Delta_IC = 0.1;

const double basal_Notch_IC = 0.1;
const double luminal_Notch_IC = 0.2;

void TestspheriodFixed3D(){

        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<3>*> nodes;
        static double Lumenial_rad = 3.5;
        static double Basal_rad = Lumenial_rad +(sqrt(3)*cell_radius);

        const double gr=(sqrt(5.0) + 1.0) / 2.0;  // golden ratio = 1.6180339887498948482
        const double ga=(2.0 - gr) * (2.0*M_PI);  // golden angle = 2.39996322972865332

        const int num_points_lum = 80;
        const int num_points_bas = num_points_lum+60;
        
        //BASAL LAYER

        for (size_t i=1; i <= num_points_bas; ++i) 
        {
        const double lat = asin(-1.0 + 2.0 * double(i) / (num_points_bas+1)) + 0.3;
        const double lon = ga * i;

        const double x = Basal_rad*cos(lon + sqrt(3))*cos(lat);
        const double y = Basal_rad*sin(lon + sqrt(3))*cos(lat);
        const double z = Basal_rad*sin(lat);

        Node<3>* p_node(new Node<3>(nodes.size(), false,x ,y,z  ));
                        p_node->SetRadius(cell_radius);
                        nodes.push_back(p_node);
        }
        
        //LUMINAL LAYER
        for (size_t i=1; i <= num_points_lum; ++i) 
        {
        const double lat = asin(-1.0 + 2.0 * double(i) / (num_points_lum+1));
        const double lon = ga * i;

        const double x = Lumenial_rad*cos(lon)*cos(lat);
        const double y = Lumenial_rad*sin(lon)*cos(lat);
        const double z = Lumenial_rad*sin(lat);

        Node<3>* p_node(new Node<3>(nodes.size(), false,x ,y,z  ));
                        p_node->SetRadius(cell_radius);
                        nodes.push_back(p_node);
        }


  

        NodesOnlyMesh<3> mesh; 
        double cut_off_length = 2; //set cut off length to be 4 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length);

        //initialise cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
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
            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
            cells.push_back(p_cell);
        }

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        for (unsigned i = 0;i<mesh.GetNumNodes(); i++){
            Node<3> *node = cell_population.GetNode(i);
            node->SetRadius(cell_radius);
        }


        //next we removed all unwanted nodes
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   


        // Get distance from centre of cell population
        c_vector<double,3> location = cell_population.GetLocationOfCellCentre(*cell_iter);

        double r = norm_2(location);

            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > Lumenial_rad)
                {   
                    cell_iter->SetCellProliferativeType(p_ME_type);
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(basal_Notch_IC);
                        initial_conditions.push_back(basal_Delta_IC);                 
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);  
               
                }
             if (r <= Lumenial_rad+1e-1)
                {

                    cell_iter->SetCellProliferativeType(p_Lpos_type);
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(luminal_Notch_IC);
                        initial_conditions.push_back(luminal_Delta_IC); 
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                  
                     
                 }  

       }

        std::cout<< cell_population.GetNumNodes() << endl;

        //create the simulation 
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Collier adaptive 3d");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.002);
        simulator.SetEndTime(100);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(CollierTrackingModifier<3>, p_modifier); //fixed weights
        p_modifier->SetPolarisationParameter(0.11*0.8);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<3>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(spring_const); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(cell_radius/2); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(2*cell_radius); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);

        MAKE_PTR(RandomMotionForce<3>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert); 
        simulator.AddForce(p_random_force);

       simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }




}

};