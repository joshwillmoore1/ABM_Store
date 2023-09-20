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

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellsGenerator.hpp"

//notch dependent adhesion
#include "DifferentialNotchAdhesionForce.hpp"

//Normal division vector
#include "NormalCentreDivisionPlaneRule.hpp"
#include "LumenPressureForce.hpp"

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestDiffAdhesionDependentNotch: public AbstractCellBasedTestSuite
{
public:

const double spring_const = 25;
const double Random_pert = 0.003;

const double Notch_hss = 0.1936;
const double Delta_hss = 0.049;

const double basal_Delta_IC = 0.2;
const double luminal_Delta_IC = 0.1;

const double basal_Notch_IC = 0.1;
const double luminal_Notch_IC = 0.2;

const int End_time = 50;
const double time_step = 0.001;
const double sample_step = 50;

const double lum_radius = 3;
const int cells_across = 21;
const int cells_up = 21;

const double preBirthTime = 20.0;

const double boundaryEnergy = 2.5; // for rounder cells at apical-basal surfaces
const double SurfaceEnergy = 0.9;
const double BasalBasalEnergy = 1;
const double LuminalLuminalEnergy = 1;
const double BasalLuminalEnergy = 1;

const double force_decay = 0.3;
const double force_scale = 0.3; //0.25

const double SeedValue = 10.0; //6


 void TestVertexNMDwithAdhesion()
    {
/*
        RandomNumberGenerator::Instance()->Reseed(SeedValue);
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
        
        //Normal to centre of the organoid division rule
        MAKE_PTR(NormalCentreDivisionPlaneRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

c_vector<double,2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
           c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location-centre);

            if (dist_to_centre < Centre_tol){
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;

            }
       }

       c_vector<double,2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

         for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

        double r = norm_2(location-centre_cell);

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
                        cell_iter->GetCellData()->SetItem("target area", 0.7);

                        double birth_time = RandomNumberGenerator::Instance()->ranf()*-preBirthTime;
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

                        double birth_time = RandomNumberGenerator::Instance()->ranf()*-preBirthTime;
                        cell_iter->SetBirthTime(birth_time);
   
                     
                 }


                 if (r>lum_radius+1.2 || r<lum_radius-0.5)
                 {
                     cell_iter -> Kill();
                 }

                 


       }
       cell_population.RemoveDeadCells();


        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VERTEX_NDM_large");

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
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(boundaryEnergy);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(1.0);
        simulator.AddForce(p_force);


        MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);


        // Run simulation
        //simulator.Solve();
*/

    }


 void TestVertexNMDwithAdhesionFixed()
    {   
        RandomNumberGenerator::Instance()->Reseed(SeedValue);

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

        //Normal to centre of the organoid division rule
        MAKE_PTR(NormalCentreDivisionPlaneRule<2>, p_DivRule);
        cell_population.SetVertexBasedDivisionRule(p_DivRule);

        c_vector<double,2> centre = cell_population.GetCentroidOfCellPopulation();

        CellPtr Centre_cellPtr;
        double Centre_tol = 1e8;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
           c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double dist_to_centre = norm_2(location-centre);

            if (dist_to_centre < Centre_tol){
                Centre_cellPtr = *cell_iter;
                Centre_tol = dist_to_centre;

            }
       }
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        

       c_vector<double,2> centre_cell = cell_population.GetLocationOfCellCentre(Centre_cellPtr);

         for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {   
        // Get distance from centre of cell population
        c_vector<double,2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

        double r = norm_2(location-centre_cell);

            // define the initial cell types - dependent on distance from the organoid centre - ALL TYPES ARE NON-PRO TO INVESTIAGE NOTCH DELTA PATTERNS
            if (r > lum_radius)
                {       
                        double pert1 = 0.05*( 2*Notch_hss*RandomNumberGenerator::Instance()->ranf() - 1.2*Notch_hss);
                        double pert2 = 0.15*( 2*Delta_hss*RandomNumberGenerator::Instance()->ranf() - 0.5*Delta_hss);
                        cell_iter->SetCellProliferativeType(p_BSC_type);
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        //initial_conditions.push_back(basal_Notch_IC);
                        //initial_conditions.push_back(basal_Delta_IC);
                        initial_conditions.push_back(Notch_hss+pert1);
                        initial_conditions.push_back(Delta_hss+pert2);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                        cell_iter->GetCellData()->SetItem("target area",0.8);

                        double birth_time = RandomNumberGenerator::Instance()->ranf()*-preBirthTime;
                        cell_iter->SetBirthTime(birth_time);

                }
             if (r <= lum_radius+1e-5)
                {       

                        double pert1 = 0.1*( 2*Notch_hss*RandomNumberGenerator::Instance()->ranf() - 0.5*Notch_hss);
                        double pert2 = 0.1*( 2*Delta_hss*RandomNumberGenerator::Instance()->ranf() - 1.2*Delta_hss);
                        cell_iter->AddCellProperty(p_label);
                        cell_iter->SetCellProliferativeType(p_Lneg_type);
                        CollierSrnModel* p_srn_model = new CollierSrnModel();
                        std::vector<double> initial_conditions;
                        initial_conditions.push_back(Notch_hss+pert1);
                        initial_conditions.push_back(Delta_hss+pert2);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);

                        double birth_time = RandomNumberGenerator::Instance()->ranf()*-preBirthTime;
                        cell_iter->SetBirthTime(birth_time);
   
                     
                 }


                 if (r>lum_radius+1.2 || r<lum_radius-0.5)
                 {
                     cell_iter -> Kill();
                 }

       }
       cell_population.RemoveDeadCells();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VERTEX_NDM_ADAPT_Bristol_1");
        
        
        MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        p_modifier->SetPolarisationParameter(0.11*0.95);
        simulator.AddSimulationModifier(p_modifier);

        //MAKE_PTR(CollierTrackingModifier<2>, p_modifier);
        //p_modifier->SetSameCellW1(0.11*0.95);
        //simulator.AddSimulationModifier(p_modifier);

        // Set time step and end time for simulation
        simulator.SetDt(time_step);
        simulator.SetSamplingTimestepMultiple(sample_step);
        simulator.SetEndTime(End_time);

        // Set up force law and pass it to the simulation -  change these for differential adhesion
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(SurfaceEnergy);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(BasalBasalEnergy);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(BasalLuminalEnergy);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(LuminalLuminalEnergy);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(boundaryEnergy);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(boundaryEnergy);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(LumenKiller<2>, p_killer,(&cell_population));
        simulator.AddCellKiller(p_killer);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(Random_pert);
        simulator.AddForce(p_random_force);
        
        //Luminal pressure 
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetMovementParameter(force_decay);
        p_lumen_force->SetScaleForceParameter(force_scale);
        simulator.AddForce(p_lumen_force);

    


        // Run simulation
        simulator.Solve();



    }


};