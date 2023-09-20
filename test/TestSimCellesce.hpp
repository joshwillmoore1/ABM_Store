#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"

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

#include "FakePetscSetup.hpp"

class TestSimulationsForCellesce : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer()
    {
        static double radius = 6;
        static double lum_rad = radius - 2;

        HoneycombVertexMeshGenerator generator(2*radius + 1, 2*radius + 1);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::cout<< "tot num of cells " <<p_mesh->GetNumElements()<< endl;

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(BasalStemCellProliferativeType, p_BSC_type);
        MAKE_PTR(MyoEpiCellProliferativeType, p_ME_type);
        MAKE_PTR(LumenERPositiveCellProliferativeType, p_Lpos_type);
        MAKE_PTR(LumenERNegativeCellProliferativeType, p_Lneg_type);
        

        OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
        p_cc_model->SetDimension(2);

        cout << "here" << endl;

        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            cout << "here" << endl;
            OrganoidG1FixedCellCycleModel* p_cc_model = new OrganoidG1FixedCellCycleModel();
            p_cc_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_cc_model));  
            p_cell->SetCellProliferativeType(p_Lpos_type);
            double birth_time = RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            p_cell->GetCellData()->SetItem("target area", 1.0);
            cells.push_back(p_cell);
        }


         cout << "here pop" << endl;
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);


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
                        initial_conditions.push_back(1-0.1);
                        initial_conditions.push_back(1-0.1);
                        initial_conditions.push_back(1+0.1);
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
                        initial_conditions.push_back(1+0.1);
                        initial_conditions.push_back(1+0.1);
                        initial_conditions.push_back(1-0.1);
                        initial_conditions.push_back(9);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(0.5);
                        initial_conditions.push_back(5);
                        p_srn_model->SetInitialConditions(initial_conditions);
                        cell_iter->SetSrnModel(p_srn_model);
                       
                     
                     
                 } 

                 if (r < lum_rad-2){
                     cell_iter->Kill();

                 }
                
                if (r > lum_rad+1){
                    cell_iter->Kill();

                 }



       }
        cell_population.RemoveDeadCells();
        std::cout<< cell_population.GetNumAllCells() << endl;


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Sim_for_cellesce");
        simulator.SetDt(1.0/500.0);
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetEndTime(10);


        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

       MAKE_PTR(BscTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();

    }
};