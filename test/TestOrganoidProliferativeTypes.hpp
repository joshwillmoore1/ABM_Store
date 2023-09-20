#ifndef TESTORGANOIDPROLIFERATIVETYPES_HPP_
#define TESTORGANOIDPROLIFERATIVETYPES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellProliferativeType.hpp"
#include "CellPropertyRegistry.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "BscSrnModel.hpp"
#include "BscTrackingModifier.hpp"

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestOrganoidProliferativeTypes : public AbstractCellBasedTestSuite
{
public:

//make the test for the basal stem
    void TestArchiveBasalStemCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BasalStemCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            BasalStemCellProliferativeType* p_type = new BasalStemCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 10u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            BasalStemCellProliferativeType* p_real_state = dynamic_cast<BasalStemCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 10u);

            // Tidy up
            delete p_type;
        }
    }

     void TestArchiveMyoEpiCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "MyoEpiCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            MyoEpiCellProliferativeType* p_type = new MyoEpiCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 40u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            MyoEpiCellProliferativeType* p_real_state = dynamic_cast<MyoEpiCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 40u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveLumenERPositiveCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "LumenERPositiveCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            LumenERPositiveCellProliferativeType* p_type = new LumenERPositiveCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 20u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            LumenERPositiveCellProliferativeType* p_real_state = dynamic_cast<LumenERPositiveCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 20u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveLumenERNegativeCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "LumenERNegativeCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            LumenERNegativeCellProliferativeType* p_type = new LumenERNegativeCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 30u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            LumenERNegativeCellProliferativeType* p_real_state = dynamic_cast<LumenERNegativeCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 30u);

            // Tidy up
            delete p_type;
        }
    }
};

#endif /*TESTORGANOIDPROLIFERATIVETYPES_HPP_*/
