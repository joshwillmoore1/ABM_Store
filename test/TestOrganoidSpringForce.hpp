#ifndef TESTORGANOIDSPRINGFORCE_HPP_
#define TESTORGANOIDSPRINGFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ChemotacticForce.hpp"
#include "RepulsionForce.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "WelikyOsterForce.hpp"
#include "FarhadifarForce.hpp"
#include "DiffusionForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "OrganoidSpringForce.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOrganoidSpringForce : public AbstractCellBasedTestSuite
{
public:

    void TestOrganoidForce()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size(), location_indices);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create force
        DifferentialAdhesionGeneralisedLinearSpringForce<2> force;

        // Test set/get method
        TS_ASSERT_DELTA(force.GetHomotypicLabelledSpringConstantMultiplier(), 1.0, 1e-6);
        TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 1.0, 1e-6);

        force.SetHomotypicLabelledSpringConstantMultiplier(2.0);
        force.SetHeterotypicSpringConstantMultiplier(4.0);

        TS_ASSERT_DELTA(force.GetHomotypicLabelledSpringConstantMultiplier(), 2.0, 1e-6);
        TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 4.0, 1e-6);

        // Initialise a vector of node forces
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point;
        old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];
        p_mesh->SetNode(59, new_point, false);

        double spring_stiffness = force.GetMeinekeSpringStiffness();

        // Test the case where node 59 and its neighbours are unlabelled
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], (-3+4.0/sqrt(7.0))*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

        // Next, test the case where node 59 is labelled but its neighbours are not...
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        cell_population.GetCellUsingLocationIndex(59)->AddCellProperty(p_label);

        force.AddForceContribution(cell_population);

        // ...for which the force magnitude should be increased by 4, our chosen multiplier for heterotypic interactions under attraction
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 4.0*0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+4.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

        // Finally, test the case where node 59 and its neighbours are labelled...
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->AddCellProperty(p_label);
        }

        force.AddForceContribution(cell_population);

        // ...for which the force magnitude should be increased by 2, our chosen multiplier for homotypic labelled interactions, again only for attractive interactions
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 2.0*0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+2.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);
    }

    void TestForceOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with GeneralisedLinearSpringForce
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "GeneralisedLinearSpringForce-2-2");

        out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("linear_results.parameters");
        linear_force.OutputForceParameters(linear_force_parameter_file);
        linear_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("linear_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/linear_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with DifferentialAdhesionGeneralisedLinearSpringForce
        DifferentialAdhesionGeneralisedLinearSpringForce<2> differential_linear_force;
        differential_linear_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(differential_linear_force.GetIdentifier(), "DifferentialAdhesionGeneralisedLinearSpringForce-2-2");

        out_stream differential_linear_force_parameter_file = output_file_handler.OpenOutputFile("differential_linear_results.parameters");
        differential_linear_force.OutputForceParameters(differential_linear_force_parameter_file);
        differential_linear_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("differential_linear_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/differential_linear_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ChemotacticForce
        ChemotacticForce<2> chemotactic_force;
        TS_ASSERT_EQUALS(chemotactic_force.GetIdentifier(), "ChemotacticForce-2");

        out_stream chemotactic_force_parameter_file = output_file_handler.OpenOutputFile("chemotactic_results.parameters");
        chemotactic_force.OutputForceParameters(chemotactic_force_parameter_file);
        chemotactic_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("chemotactic_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/chemotactic_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with RepulsionForce
        RepulsionForce<2> repulsion_force;
        TS_ASSERT_EQUALS(repulsion_force.GetIdentifier(), "RepulsionForce-2");

        out_stream repulsion_force_parameter_file = output_file_handler.OpenOutputFile("repulsion_results.parameters");
        repulsion_force.OutputForceParameters(repulsion_force_parameter_file);
        repulsion_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("repulsion_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/repulsion_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with NagaiHondaForce
        NagaiHondaForce<2> nagai_force;
        TS_ASSERT_EQUALS(nagai_force.GetIdentifier(), "NagaiHondaForce-2");

        out_stream nagai_force_parameter_file = output_file_handler.OpenOutputFile("nagai_results.parameters");
        nagai_force.OutputForceParameters(nagai_force_parameter_file);
        nagai_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("nagai_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/nagai_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with NagaiHondaDifferentialAdhesionForce
        NagaiHondaDifferentialAdhesionForce<2> nagai_da_force;
        TS_ASSERT_EQUALS(nagai_da_force.GetIdentifier(), "NagaiHondaDifferentialAdhesionForce-2");

        out_stream nagai_da_force_parameter_file = output_file_handler.OpenOutputFile("nagai_da_results.parameters");
        nagai_da_force.OutputForceParameters(nagai_force_parameter_file);
        nagai_da_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("nagai_da_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/nagai_da_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with WelikyOsterForce
        WelikyOsterForce<2> weliky_force;
        TS_ASSERT_EQUALS(weliky_force.GetIdentifier(), "WelikyOsterForce-2");

        out_stream weliky_force_parameter_file = output_file_handler.OpenOutputFile("weliky_results.parameters");
        weliky_force.OutputForceParameters(weliky_force_parameter_file);
        weliky_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("weliky_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/weliky_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        FarhadifarForce<2> force;
        TS_ASSERT_EQUALS(force.GetIdentifier(), "FarhadifarForce-2");

        out_stream farhadifar_force_parameter_file = output_file_handler.OpenOutputFile("farhadifar_results.parameters");
        force.OutputForceParameters(farhadifar_force_parameter_file);
        farhadifar_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("farhadifar_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/farhadifar_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with DiffusionForce
        DiffusionForce<2> diffusion_force;
        TS_ASSERT_EQUALS(diffusion_force.GetIdentifier(), "DiffusionForce-2");

        out_stream diffusion_force_parameter_file = output_file_handler.OpenOutputFile("diffusion_results.parameters");
        diffusion_force.OutputForceParameters(diffusion_force_parameter_file);
        diffusion_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("diffusion_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/diffusion_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
    }


#endif /*TESTORGANOIDSPRINGFORCE_HPP_*/
