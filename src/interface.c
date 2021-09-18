#include "inter_folder/interface.h"

#include "jurassic.h"
#include <stdio.h>


//void jur_formod(ctl_t const *ctl, atm_t *atm, obs_t *obs);

void formod_multiple_packages(ctl_t *ctl, atm_t *atm,  aero_t *aero, int n, obs_t *obs_packages, los_t **los_packages) { 
  printf("DEBUG number of packages.. %d\n", n);
  printf("DEBUG call jur_formod..\n");
	jur_formod(ctl, atm, obs_packages, aero, los_packages, n);
}
