#include "TACellMutationState.hpp"

TACellMutationState::TACellMutationState()
    : AbstractCellMutationState(3)
{}


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TACellMutationState)
