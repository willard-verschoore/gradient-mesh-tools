#pragma once

#include "slot-map.hpp"

/**
 * This file contains type aliases for the storage containers used by the
 * gradient mesh data structure. This allows us some flexibility to change them
 * later on without having to alter significant parts of the codebase.
 */

template <typename T>
using Storage = SlotMap<T>;
template <typename T>
using Id = SlotMapHandle<T>;
