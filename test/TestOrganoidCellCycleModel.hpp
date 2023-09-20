#ifndef TESTORGANOIDCELLCYCLEMODELS_HPP_
#define TESTORGANOIDCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "ApoptoticCellProperty.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellCycleTimesGenerator.hpp"
#include "CellLabel.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "GammaG1CellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "OrganoidG1FixedCellCycleModel.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestOrganoidCellCycleModels : public AbstractCellBasedTestSuite
{
public:


void TestOrganoidG1FixedCellCycleModel()
    {
        //write some tests
    }
};

#endif /*TESTORGANOIDCELLCYCLEMODELS_HPP_*/
