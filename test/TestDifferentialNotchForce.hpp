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
#include "CellLabel.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

//notch dependent adhesion
#include "DifferentialNotchAdhesionForce.hpp"

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestDiffAdhesionDependentNotch: public AbstractCellBasedTestSuite
{
public:

const double spring_const = 25;
const double Random_pert = 0.0003;

const double basal_Delta_IC = 0.2;
const double luminal_Delta_IC = 0.1;

const double basal_Notch_IC = 0.1;
const double luminal_Notch_IC = 0.2;

const int End_time = 300;
const double time_step = 0.005;
const double sample_step = 100;

void Test2dBilayerAdaptiveAdhesion()
    {
        //set up the geometry in a circle hex grid in Radial coordinates
                /*
                      o---o---o
                     /   /   /
                    o---o---o                          
                */  

     
        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<2>*> nodes;
        static double lum_radius = 2*sqrt(2);
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
            p_cell->SetCellProliferativeType(p_Lpos_type);
            double birth_time = RandomNumberGenerator::Instance()->ranf()*2.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

      
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
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(basal_Notch_IC);
                        initial_conditions.push_back(basal_Delta_IC);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);

                          

                }
             if (r <= lum_radius+1e-5)
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

        //create the simulation 
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("With_adhesion");
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetDt(time_step);
        simulator.SetEndTime(End_time);
        
        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11*0.95);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

     
        //new force law that uses notch values to encourage homophilic adhesion for notch activation
        MAKE_PTR(DifferentialNotchAdhesionForce<2>, p_differential_adhesion_force); 
        p_differential_adhesion_force->SetMeinekeSpringStiffness(spring_const); 
        p_differential_adhesion_force->SetMaximumSpringConstantMultiplier(2); //undetermined
        p_differential_adhesion_force->SetAdhesionDecayConstant(5); //undetermined
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


    void Test2dBilayerAdaptiveNoAdhesion()
    {
        //set up the geometry in a circle hex grid in Radial coordinates
                /*
                      o---o---o
                     /   /   /
                    o---o---o                          
                */  

     
        double cell_radius = 0.5; // this implies the spring rest length is 1 
        std::vector<Node<2>*> nodes;
        static double lum_radius = 2*sqrt(2);
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
            p_cell->SetCellProliferativeType(p_Lpos_type);
            double birth_time = RandomNumberGenerator::Instance()->ranf()*2.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        

        //define the cell population - putting the cells into the mesh
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

      
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
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(basal_Notch_IC);
                        initial_conditions.push_back(basal_Delta_IC);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                          

                }
             if (r <= lum_radius+1e-5)
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
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("No_adhesion");
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetDt(time_step);
        simulator.SetEndTime(End_time);
        
        //this modifier is how the SRN network interactions with the cells
        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11*0.95);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

     
        MAKE_PTR(OrganoidSpringForce<2>, p_differential_adhesion_force); // the organoid force just considers cell types not labels
        p_differential_adhesion_force->SetMeinekeSpringStiffness(spring_const); //my values (not justified) 20
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