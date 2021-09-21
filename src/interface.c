#include "inter_folder/interface.h"

#include "jurassic.h"
#include <stdio.h>


//void jur_formod(ctl_t const *ctl, atm_t *atm, obs_t *obs);

void formod_multiple_packages(ctl_t *ctl, atm_t *atm,  aero_t *aero, int n, obs_t *obs_packages, los_t ***los_packages) { 
  printf("DEBUG #%d number of packages.. %d\n", ctl->MPIglobrank, n);
  printf("DEBUG #%d call jur_formod..\n", ctl->MPIglobrank);
	jur_formod(ctl, atm, obs_packages, aero, los_packages, n);
}
