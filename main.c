#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "misc.h"
#include "satellite.h"
#include "maidenhead.h"
#include "linux.h"
#include "sattests.h"


void print_pos(double jd, topo_pos_t obs, sat_pos_t p){
    printf("NORAD:\t%d\tErrors:\t%d\n", p.satid, p.err);
    printf("Az:   \t%.2f\tEl:\t%.2f\tRange:\t%.1f\n", p.az, p.elev, p.dist);
    printf("\t\t\t\tRaDec \t%.4f %.4f\n", p.ra, p.dec);
}
void print_n2yo(tle_t tle){
    printf("https://n2yo.com/?s=%d\n", tle.norad_number);
}

char algo[33] = {0};

int main(int argc, char **argv){
    topo_pos_t obs = {0};
    latlon me = {0}; //need this because topo_pos_t is floats but maidenhead uses doubles
    me = maidenhead_to_latlon("FN41");
    obs.lat = me.lat;
    obs.lon = me.lon;

    sat_mem_t satellites[10] = {0};
    makelut();
    test_1();
    printf("LL: %.8f\t%.8f\n", obs.lat, obs.lon);
    return 0;
}
