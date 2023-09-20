#ifndef TESTLUMENKILLER_HPP_
#define TESTLUMENKILLER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"
#include "CellsGenerator.hpp"
#include "TargetedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "ApoptoticCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "IsolatedLabelledCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "TrianglesMeshReader.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "SmartPointers.hpp"
#include "LumenKiller.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellDeltaNotchWriter.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

#include "DeltaNotchSrnModel.hpp"
#include "DeltaNotchTrackingModifier.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellKiller.
 */

 class TestLumenKillers : public AbstractCellBasedTestSuite
{
public:
void TestLumenKillerMthod()
    {
        static double ORGANOID_DOMAIN = 2; //31
        static double ORGANOID_RADIUS = (ORGANOID_DOMAIN-1)/2 -1; //14

        std::vector<Node<3>*> nodes;

        //first we make a cubic lattice of nodes

        for (unsigned i=0; i<ORGANOID_DOMAIN; i++)
            {
                for(unsigned j=0; j<ORGANOID_DOMAIN; j++ )
                {
                    for (unsigned k = 0; k<ORGANOID_DOMAIN; k++)
                    {
                        nodes.push_back(new Node<3>(nodes.size(), false, i, j, k));
                    }
                }
            }

           
        NodesOnlyMesh<3> mesh; //make the mesh object to pass the nodes into

        mesh.ConstructNodesWithoutMesh(nodes, 1);

        

        std::vector<CellPtr> cells;

        //pointers for the cell types and states
       //pointers for the cell types and states
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(CellLabel, p_label);

        //For each node - set: cycle type, initial conditons for D/N ODEs and state 
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(3);

            std::vector<double> initial_conditions;
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf()); //random initial conditons 
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel(); //make a pointer to the ODE model
            p_srn_model->SetInitialConditions(initial_conditions); // input the initial conditons

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model)); // make a cell pointer with the s
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
            
        }

         //Create the cell population from the mesh and cells that were defined above

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        //next we removed all unwanted nodes
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
       {    

        // Get distance from centre of cell population
        c_vector<double,3> location = cell_population.GetLocationOfCellCentre(*cell_iter);
       
        double r = pow(location[0]-(ORGANOID_DOMAIN-1)/2,2) + pow(location[1]-(ORGANOID_DOMAIN-1)/2,2) + pow(location[2]-(ORGANOID_DOMAIN-1)/2,2);

            // define the initial cell types - dependent on distance from the organoid centre
            if (r > pow(ORGANOID_RADIUS-1,2))
                {   
                    if (RandomNumberGenerator::Instance()->ranf()>0.8)
                     {
                        cell_iter->SetCellProliferativeType(p_stem_type);
                       
                     }
                     else
                     {
                        cell_iter->SetCellProliferativeType(p_transit_type);
                       
                        
                      }       

                }
             if (r <= pow(ORGANOID_RADIUS-1,2))
                {
                 cell_iter->SetCellProliferativeType(p_diff_type);
                 cell_iter->AddCellProperty(p_label);
               
                 
                 }
            
             if (r > pow(ORGANOID_RADIUS,2))
                {
                  cell_iter->Kill();
                    
                }
             if (r < pow(ORGANOID_RADIUS-2,2))
                 {
                  //cell_iter->Kill();

                 }    

       }
        cell_population.RemoveDeadCells();

        std::cout<< cell_population.GetNumNodes() << endl;

       // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
     
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("LumenKillerTest");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(50);
        simulator.SetDt(0.005);

        MAKE_PTR_ARGS(LumenKiller<3>, p_killer,(&cell_population));
        simulator.AddCellKiller(p_killer);
        
        
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force); //**Changed**//
        simulator.AddForce(p_force);


        simulator.Solve();

        //clean up for memory management 
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }


};

#endif /*TESTLUMENKILLER_HPP_*/
