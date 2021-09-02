#include "inter_folder/interface.h"

#include "jurassic.h"
#include <stdio.h>


void formod_GPU(ctl_t const *ctl, atm_t *atm, obs_t *obs);

void formod_multiple_packages(ctl_t *ctl, atm_t *atm, obs_t *packages, int n) {
  printf("number of packages.. %d\n", n);
  for(int i = 0; i < n; i++) {
    printf("call formod_GPU..\n");
    formod_GPU(ctl, atm, &packages[i]);
  }
}
