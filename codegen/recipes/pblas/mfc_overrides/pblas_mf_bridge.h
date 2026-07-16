/* pblas_mf_bridge.h -- thin wrapper that pulls in the unified
 * multifloats_bridge.h. Kept as a separate file so the PBLAS recipe
 * and compile_pblas.sh do not need to change their -include flag. */
#ifndef PBLAS_MF_BRIDGE_H
#define PBLAS_MF_BRIDGE_H
#include "multifloats_bridge.h"
#endif
