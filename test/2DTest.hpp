#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "LumenKiller.hpp"
#include "CellLabel.hpp"
#include "BscSrnModel.hpp"
#include "BscTrackingModifier.hpp"
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DeltaNotchSrnModel.hpp"
#include "DeltaNotchTrackingModifier.hpp"
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


class TestRunningNodeBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayerForBscModel()
    {
        
        //geometry setup in polar coordinates
        static double rad = 1;
        static double Lumenial_rad = 80;
        static double Basal_rad = Lumenial_rad+2;
        static double basal_node_dist = 1;
        static double arc_pos  = round(2*M_PI/(2*asin(basal_node_dist*0.5/(Basal_rad-1))));
        static double lum_dist = basal_node_dist*(Lumenial_rad)/(Lumenial_rad+1);
        std::cout<< "Lum dis = " << lum_dist << endl;
        std::vector<Node<2>*> nodes;

        for (unsigned r=Lumenial_rad; r<Basal_rad; r++)
            {
                for (unsigned t = 0; t<arc_pos; t++)
                {
                    if  (r == Basal_rad-1)
                    {
                        Node<2>* p_node(new Node<2>(nodes.size(), true, r*cos((2*M_PI)*(t/arc_pos) ), r*sin((2*M_PI)*(t/arc_pos) )));
                        p_node->SetRadius(rad);
                        nodes.push_back(p_node);
                    }
                    else
                    {
                        Node<2>* p_node(new Node<2>(nodes.size(), false, r*cos((2*M_PI)*(t/arc_pos) ), r*sin((2*M_PI)*(t/arc_pos) )));
                        p_node->SetRadius(rad);
                        nodes.push_back(p_node);
                    }
                    
                    
                }

            }
            
        NodesOnlyMesh<2> mesh; 
        double cut_off_length = 5; 
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length);

        //initialise cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
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
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
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
            if (r > Lumenial_rad)
                {   
                    
                       cell_iter->SetCellProliferativeType(p_ME_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.1);
                        initial_conditions.push_back(0.1);
                        initial_conditions.push_back(4);
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
        simulator.SetOutputDirectory("Non-equi 2d bi polar-3connect");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(100);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

       // MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population)); //depends on cell type
        //simulator.AddCellKiller(p_killer);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(0); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.4); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(1); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);


        simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        

        
    }



//THIS IS A SECOND TEST OF A SPHERIOD


    void Test3Dspheroid()
    {
        //geometry setup in spherical coordinates
        static double Lumenial_rad = 4;
        static double Basal_rad = Lumenial_rad+1;
        static double basal_node_dist = 1;
        static double arc_pos  = ceil(2*M_PI/(2*asin(basal_node_dist*0.5/(Basal_rad-1))));
        static double arc_pos_theta  = ceil(arc_pos*0.5);
        static double lum_dist = basal_node_dist*(Lumenial_rad)/(Lumenial_rad+1);
        std::cout<< "theta = " << arc_pos_theta  << endl;
        std::cout<< "Lum dis = " << lum_dist << endl;
        std::vector<Node<3>*> nodes;

        const double gr=(sqrt(5.0) + 1.0) / 2.0;  // golden ratio = 1.6180339887498948482
        const double ga=(2.0 - gr) * (2.0*M_PI);  // golden angle = 2.39996322972865332

        const int num_points_lum = 50;
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
                        p_node->SetRadius(1.5);
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
                        p_node->SetRadius(1.5);
                        nodes.push_back(p_node);
        }

            
        NodesOnlyMesh<3> mesh; //make the mesh object to pass the nodes into
        double cut_off_length = 1.5; //this value was taken from the comparing paper 
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length); //pass the nodes into the mess and defining the node connectivity radius


        //Geometry setup in grid format
        /*HoneycombMeshGenerator generator(2,40);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        */

        //initialise cells
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

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<3> cell_population(mesh, cells);


        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,3> location = cell_population.GetLocationOfCellCentre(*cell_iter);
      

        double r =  norm_2(location);
    

            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > Lumenial_rad)
                {   
                    if (RandomNumberGenerator::Instance()->ranf()>0.8)
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
                     else
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

                }
             if (r <= Lumenial_rad)
                {
                    
                     if (RandomNumberGenerator::Instance()->ranf()>0.5)
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
                     else
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

       }


        std::cout<< cell_population.GetNumNodes() << endl;

        //create the simulation 
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Bilayered Spheriod 3D");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10);
        

        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(BscTrackingModifier<3>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

       // MAKE_PTR_ARGS(LumenKiller<3>, p_killer,(&cell_population)); //depends on cell type
        // simulator.AddCellKiller(p_killer);

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<3>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(20); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.3); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(cut_off_length); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);


        //simulator.Solve();

        std::cout<<"Number of cells in simulation = "<< cell_population.GetNumNodes() << endl;

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        

    
    }


    void Test2dBilayer()
    {


        static double Length = 26; //31
        static double Height = 2; //14
        std::vector<Node<2>*> nodes;
        //first we make a cubic lattice of nodes

        for (unsigned h=0; h<Height; h++)
            {
                for(unsigned l=0; l<Length; l++ )
                {
                        Node<2>* p_node(new Node<2>(nodes.size(), false, l, h));
                        p_node->SetRadius(1);
                        nodes.push_back(p_node);
                }
            }

        NodesOnlyMesh<2> mesh; //make the mesh object to pass the nodes into
        double cut_off_length = 1; //this value was taken from the comparing paper 
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length); //pass the nodes into the mess and defining the node connectivity radius


        std::vector<CellPtr> cells;

        //pointers for the cell types and states
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);
       
      

        //For each node - set: cycle type, initial conditons for D/N ODEs and state 
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
            
        }

         //Create the cell population from the mesh and cells that were defined above

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMeinekeDivisionSeparation(0.4);
        

        //next we removed all unwanted nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
        
        

        if (location[1] == 0)
                {   
                    if (RandomNumberGenerator::Instance()->ranf()>0.8)
                     {
                       cell_iter->SetCellProliferativeType(p_Lpos_type);
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
                        cell_iter->SetSrnModel(p_srn_model);
                        
                     }
                     else
                     {
                        cell_iter->SetCellProliferativeType(p_Lpos_type);
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
                        cell_iter->SetSrnModel(p_srn_model);
                     }       

                }
             if (location[1] == 1)
                {
                    
                     if (RandomNumberGenerator::Instance()->ranf()>0.5)
                     {
                        cell_iter->SetCellProliferativeType(p_ME_type);
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
                        cell_iter->SetSrnModel(p_srn_model);
                        
                     }
                     else
                     {
                        cell_iter->SetCellProliferativeType(p_ME_type);
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
                        cell_iter->SetSrnModel(p_srn_model);
                       
                     }
                     
                 }  
           

       }
        cell_population.RemoveDeadCells();

        std::cout<< cell_population.GetNumNodes() << endl;

       // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();

     
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("2D bi Connect");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(50);
         std::cout<<"This has been updated"<<endl;


        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population)); //depends on cell type
        simulator.AddCellKiller(p_killer);
        

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.4); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(cut_off_length); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);


       //simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    
    }


     void TestBoundaryDelta()
    {


        static double Length = 26; //31
        static double Height = 2; //14
        std::vector<Node<2>*> nodes;
        //first we make a cubic lattice of nodes

        for (unsigned h=0; h<Height; h++)
            {
                for(unsigned l=0; l<Length; l++ )
                {
                        Node<2>* p_node(new Node<2>(nodes.size(), false, l, h));
                        p_node->SetRadius(1);
                        nodes.push_back(p_node);
                }
            }

        NodesOnlyMesh<2> mesh; //make the mesh object to pass the nodes into
        double cut_off_length = 1; //this value was taken from the comparing paper 
        mesh.ConstructNodesWithoutMesh(nodes, cut_off_length); //pass the nodes into the mess and defining the node connectivity radius


        std::vector<CellPtr> cells;

        //pointers for the cell types and states
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);
       
      

        //For each node - set: cycle type, initial conditons for D/N ODEs and state 
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_ME_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
            
        }

         //Create the cell population from the mesh and cells that were defined above

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMeinekeDivisionSeparation(0.4);
        

        //next we removed all unwanted nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
        
        

        if (location[1] == 0)
                {   
                   
                        cell_iter->SetCellProliferativeType(p_ME_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                           

                }
             if (location[1] == 1)
                {
                    cell_iter->SetCellProliferativeType(p_ME_type);
                        BscSrnModel* p_srn_model = new BscSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                     
                     
                 }  

                 
           

       }
        cell_population.RemoveDeadCells();

        std::cout<< cell_population.GetNumNodes() << endl;

       // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();

     
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("boundary delta");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(1);
        
       


        MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population)); //depends on cell type
        simulator.AddCellKiller(p_killer);
        

        // Create a force law and pass it to the simulation - differential adhesion - values have been taken from Obsborne et al 2017
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50); //my values (not justified) 20
        p_differential_adhesion_force->SetHomotypicTypeSpringConstantMultiplier(1); //my values (not justified) 5
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1); //my values (not justified) 1
        p_differential_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.4); //my values (not justified) 0.4
        p_differential_adhesion_force->SetCutOffLength(cut_off_length); //my values (not justified) 1
        simulator.AddForce(p_differential_adhesion_force);


       //simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }


    }

};
