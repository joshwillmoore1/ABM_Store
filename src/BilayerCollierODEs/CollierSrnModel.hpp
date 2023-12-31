#ifndef COLLIERSRNMODEL_HPP_
#define COLLIERSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CollierOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

/**
 * A subclass of AbstractOdeSrnModel that includes a Basal stem cell ODE system in the sub-cellular reaction network.
 *
 * \todo #2752 document this class more thoroughly here
 */
class CollierSrnModel : public AbstractOdeSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:
    /**
     * Protected copy-constructor for use by CreateSrnModel().  The only way for external code to create a copy of a SRN model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent SRN model will have had ResetForDivision() called just before CreateSrnModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel  the SRN model to copy.
     */
    CollierSrnModel(const CollierSrnModel& rModel);

public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    CollierSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of
     * this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    AbstractSrnModel* CreateSrnModel();

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new MamDeltaNotch ODE system.
     */
    void Initialise(); // override

    /**
     * Overridden SimulateToTime() method for custom behaviour.
     *
     * \todo #2752 say what it does in this class
     */
    void SimulateToCurrentTime();

    /**
     * Update the current levels of Delta and Notch in the cell.
     *
     * N.B. Despite the name, this doesn't update the levels of delta or notch, or compute mean levels.
     * It just copies the current mean delta from the CellData
     * (set by DeltaNotchTrackingModifier) to the DeltaNotchOdeSystem.
     *
     * \todo #2752 Improve the name of this method!
     */
    void UpdateCollier();

    /**
     * @return the current NICD level in this cell.
     */
    double GetNotch();

    /**
     * @return the current Delta level in this cell.
     */
    double GetDelta();

    /**
     * @return the current level of mean Delta in the neighbouring cells.
     *
     * N.B. This doesn't calculate anything, it just returns the parameter
     * from the DeltaNotchOdeSystem.
     */
    double GetMeanNeighbouringDelta();

     
    /**
     * Output SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSrnModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CollierSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CollierSrnModel)

#endif /* COLLIERSRNMODEL_HPP_ */
