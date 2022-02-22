#include "PanethCellMutationState.hpp"

PanethCellMutationState::PanethCellMutationState()
    : AbstractCellMutationState(4)
{}


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PanethCellMutationState)
