#include <stdio.h>
#include <stdlib.h>
#include "dogleg.h"

int main(void)
{
  dogleg_parameters2_t p = {.debug_vnlog = true};
  if(p.dogleg_debug != DOGLEG_DEBUG_VNLOG)
  {
    printf("ERROR: DOGLEG_DEBUG_VNLOG bit alignment does NOT match libdogleg v0.16\n");
    return 1;
  }
  printf("OK: DOGLEG_DEBUG_VNLOG bit alignment DOES match libdogleg v0.16\n");

  return 0;
}
