#pragma once
#include "ResponseDebugLogger.hpp"

#define DEBUG_LOG_VALUE(world, logger, key, expression)                  \
  do {                                                                         \
    TimedValueLogger __timed_logger__(world, key, logger);                     \
    auto __timed_val__ = (expression);                                         \
    __timed_logger__.log_value(__timed_val__);                                 \
  } while (0)

#define DEBUG_TIMED_BLOCK(world, logger, key, block)                           \
  do {                                                                         \
    TimedValueLogger __timed_logger__(world, key, logger);                     \
    block __timed_logger__.log();                                              \
  } while (0)
