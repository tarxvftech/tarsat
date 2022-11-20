#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "satellite.h"
#include "maidenhead.h"

int main(int argc, char **argv){
   size_t num_satellites = 0;
   topo_pos_t obs = {0};
   latlon me = {0}; //need this because topo_pos_t is floats but maidenhead uses doubles
   me = maidenhead_to_latlon("FN41");
   obs.lat = me.lat;
   obs.lon = me.lon;

   sat_mem_t satellites[10] = {0};
   parse_elements( 
         "1 22825U 93061C   22321.58780667  .00000066  00000-0  43288-4 0  9996", 
         "2 22825  98.8716 355.7510 0007856 331.8797  28.1960 14.30210513520047", 
         &(satellites[0].tle));
   strncpy(satellites[0].name, "AO27", 15);
   satellites[0].name[15] = 0;
   num_satellites += 1;
   return 0;
}
