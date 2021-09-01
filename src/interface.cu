#include "inter_folder/interface.h"

#include "jurassic.h"
#include <stdio.h>

__host__ void print_obs_size(obs_t *p) {
  printf("jurassic-gpu OBS size: %d\n", p->nr);
}
