#include "EnterocyteCellMutationState.hpp"

EnterocyteCellMutationState::EnterocyteCellMutationState()
    : AbstractCellMutationState(6)
{}


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(EnterocyteCellMutationState)
