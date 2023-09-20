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

#include "NeuGridBscTrackingModifier.hpp"
#include "MooreGridBscTrackingModifier.hpp"
#include "HexBscTrackingModifier.hpp"
#include "RandomMotionForce.hpp"
#include "BiBuckleRandomMotionForce.hpp"
#include <vector>

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPeriodicBilayer: public AbstractCellBasedTestSuite
{
public:

void Test2dGrid()
{

    //set up the geometry in a stacked grid in cartesian coordinates
        static double Length = 30; 
        static double Height = 2; 
        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<2>*> nodes;
        
        //first we make a cubic lattice of nodes - cells a places in a regular lattice with 1 dist 
                /*
                        o           o  o  o
                        |            \ | / 
                    o---o---o      o---o---o                           
                */  

        for (unsigned h=0; h<Height; h++)
            {
                for(unsigned l=0; l<Length; l++ )
                {
                        Node<2>* p_node(new Node<2>(nodes.size(), false, l, h));
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

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);


        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);              

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
           
            cells.push_back(p_cell);
        }

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        for (unsigned i = 0;i<mesh.GetNumNodes(); i++){
            Node<2> *node = cell_population.GetNode(i);
            node->SetRadius(cell_radius);
        }


        //modify cell type dependent on which layer the cell is in:

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

        if (location[1] == 1)
                {   
                   
                        cell_iter->SetCellProliferativeType(p_ME_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(0.6);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                           

                }
             if (location[1] == 0)
                {
                    cell_iter->SetCellProliferativeType(p_Lpos_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.6);
                        initial_conditions.push_back(0.6);
                        initial_conditions.push_back(0.4);
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
        simulator.SetOutputDirectory("Salt Pep");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(150);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(MooreGridBscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.8); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.4); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(1.5); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);

        MAKE_PTR(BiBuckleRandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.01); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

        //simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

}



void TestHexBilayer(){


    //set up the geometry in a stacked grid in cartesian coordinates
        static double Length = 30; 
        static double Height = 2; 
        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<2>*> nodes;
        
        //first we make a hex lattice of nodes - cells a places in a regular lattice with 1 dist in the x directon
                /*
                      o---o---o
                     /   /   /
                    o---o---o                          
                */  

        for (unsigned h=0; h<Height; h++)
            {
                for(unsigned l=0; l<Length; l++ )
                {
                        Node<2>* p_node(new Node<2>(nodes.size(), false, l + h*(cell_radius), h*(sqrt(3)*cell_radius)));
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

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);


        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);              

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
           
            cells.push_back(p_cell);
        }


        


        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        for (unsigned i = 0;i<mesh.GetNumNodes(); i++){
            Node<2> *node = cell_population.GetNode(i);
            node->SetRadius(cell_radius);
        }


        //modify cell type dependent on which layer the cell is in:

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

        if (location[1] == 0 )
                {   
                   
                

                        cell_iter->SetCellProliferativeType(p_Lpos_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.6);
                        initial_conditions.push_back(0.6);
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                           

                }
             else
                {
                    cell_iter->SetCellProliferativeType(p_ME_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(0.6);
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
        simulator.SetOutputDirectory("Periodic hex test SP");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(150);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(HexBscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.8); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(cell_radius/2); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(1.5); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);

        MAKE_PTR(BiBuckleRandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.01); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

        simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

}


void TestCirleBilayer(){


    //set up the geometry in a circle hex grid in cartesian coordinates
                /*
                      o---o---o
                     /   /   /
                    o---o---o                          
                */  

     
        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<2>*> nodes;
        static double lum_radius = 1.5;
        static double bas_radius = lum_radius +(sqrt(3)*cell_radius);

        static double d_theta_lum = 2*cell_radius/lum_radius;
        static double d_theta_bas = 2*cell_radius/bas_radius;

        static double cell_num_lum = ceil(2*M_PI/d_theta_lum);
        static double cell_num_bas = ceil(2*M_PI/d_theta_bas);

        static double bas_shift = atan(cell_radius/bas_radius);

        for (unsigned i = 0; i< cell_num_lum; i++){
            Node<2>* p_node(new Node<2>(nodes.size(), false, lum_radius*cos(2*M_PI*(i/cell_num_lum)), lum_radius*sin(2*M_PI*(i/cell_num_lum))   ));
            p_node->SetRadius(cell_radius);
            nodes.push_back(p_node);
        }

         for (unsigned i = 0; i< cell_num_bas; i++){
            Node<2>* p_node(new Node<2>(nodes.size(), false, bas_radius*cos(2*M_PI*(i/cell_num_bas) + bas_shift ), bas_radius*sin(2*M_PI*(i/cell_num_bas) + bas_shift)   ));
            p_node->SetRadius(cell_radius);
            nodes.push_back(p_node);
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

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);


        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);              

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
           
            cells.push_back(p_cell);
        }


        


        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        for (unsigned i = 0;i<mesh.GetNumNodes(); i++){
            Node<2> *node = cell_population.GetNode(i);
            node->SetRadius(cell_radius);
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
            if (r > lum_radius)
                {   
                    
                       cell_iter->SetCellProliferativeType(p_ME_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(0.4);
                        initial_conditions.push_back(2);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                        
               
                }
             if (r <= lum_radius+1e-3)
                {
                        cell_iter->SetCellProliferativeType(p_Lpos_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.6);
                        initial_conditions.push_back(0.6);
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
        simulator.SetOutputDirectory("small 2d organoid");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(50);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(cell_radius/2); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(1.5); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);

        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.01); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

       //simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

}




void Testspheriod3D(){

        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<3>*> nodes;
        static double Lumenial_rad = 3;
        static double Basal_rad = Lumenial_rad +(sqrt(3)*cell_radius);

        const double gr=(sqrt(5.0) + 1.0) / 2.0;  // golden ratio = 1.6180339887498948482
        const double ga=(2.0 - gr) * (2.0*M_PI);  // golden angle = 2.39996322972865332

        const int num_points_lum = 60;
        const int num_points_bas = num_points_lum*3;
        
        //BASAL LAYER

        for (size_t i=1; i <= num_points_bas; ++i) 
        {
        const double lat = asin(-1.0 + 2.0 * double(i) / (num_points_bas+1)) + 0.3;
        const double lon = ga * i;

        const double x = Basal_rad*cos(lon)*cos(lat);
        const double y = Basal_rad*sin(lon)*cos(lat);
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
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.8);
                        initial_conditions.push_back(0.8);
                        initial_conditions.push_back(1.2);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);  
               
                }
             if (r <= Lumenial_rad+1e-3)
                {

                    cell_iter->SetCellProliferativeType(p_Lpos_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(1.2);
                        initial_conditions.push_back(1.2);
                        initial_conditions.push_back(0.8);
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
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("three_D_ND_polar");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(50);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(BscTrackingModifier<3>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<3>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.8); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(cell_radius/2); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(0.3); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);

        MAKE_PTR(RandomMotionForce<3>, p_random_force);
        p_random_force->SetMovementParameter(0.005); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

       //simulator.Solve();
       std::cout<< cell_population.GetNumNodes() << endl;

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }




}







};
  