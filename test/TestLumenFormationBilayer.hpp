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
#include "CellLabel.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "RandomMotionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellsGenerator.hpp"

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestLumenFormation: public AbstractCellBasedTestSuite
{
public:

const double spring_const = 25;
const double Random_pert = 0.005;

const double basal_Delta_IC = 0.2;
const double luminal_Delta_IC = 0.1;

const double basal_Notch_IC = 0.1;
const double luminal_Notch_IC = 0.2;

const int End_time = 50;
const double time_step = 0.001;
const double sample_step = 50;

const double lum_radius = 2;
const int cells_across = 10;
const int cells_up = 10;

const double boundaryEnergy = 1;
const double LumLumAdhesionConst = 1;
const double BasBasAdhesionConst = 1;
const double LumBasAdhesionConst = 1;

 void TestVertexWithDiffForce()
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);

        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<LumenERPositiveCellProliferativeType>());
        CellsGenerator<OrganoidG1FixedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);

        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        //pointer for cell labels
         boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

         for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

        CellPtr Centre_cell_ptr = cell_population.GetCellUsingLocationIndex((cells_up*cells_across/2) - cells_across/2);
        c_vector<double,2> centre = cell_population.GetLocationOfCellCentre(Centre_cell_ptr);


        double r = norm_2(location-centre);

            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
                {   

                        cell_iter->SetCellProliferativeType(p_BSC_type);
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(basal_Notch_IC);
                        initial_conditions.push_back(basal_Delta_IC);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);

                        double birth_time = RandomNumberGenerator::Instance()->ranf()*-12.0;
                        cell_iter->SetBirthTime(birth_time);

                }
             if (r <= lum_radius+1e-5)
                {
                        cell_iter->SetCellProliferativeType(p_Lneg_type);
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(luminal_Notch_IC);
                        initial_conditions.push_back(luminal_Delta_IC);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                        double birth_time = RandomNumberGenerator::Instance()->ranf()*-12.0;
                        cell_iter->SetBirthTime(birth_time);
                        
                        // add label to the luminal cells to create diff adhesion  
                        cell_iter->AddCellProperty(p_label);
   
                     
                 }


                 if (r>lum_radius+1 || r<lum_radius-0.8)
                 {
                     cell_iter -> Kill();
                 }

       }
       cell_population.RemoveDeadCells();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("test_vertex_forces");

        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11*0.95);
        simulator.AddSimulationModifier(p_modifier);

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        // Set up force law and pass it to the simulation -  change these for differential adhesion
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(BasBasAdhesionConst);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(LumBasAdhesionConst);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(LumLumAdhesionConst);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(boundaryEnergy);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(1.0);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

       // MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population));
       // simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);

        //add force to introduce lumen pressure


        // Run simulation
        simulator.Solve();


    }
};