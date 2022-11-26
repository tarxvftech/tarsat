#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "satellite.h"
#include "maidenhead.h"


#define Pf(t) printf("%s = %f\n", #t, t)

void print_pos(double jd, topo_pos_t obs, sat_pos_t p){
    printf("NORAD:\t%d\tErrors:\t%d\n", p.satid, p.err);
    printf("Az:   \t%.2f\tEl:\t%.2f\tRange:\t%.1f\n", p.az, p.elev, p.dist);
    printf("\t\t\t\tRaDec \t%.4f %.4f\n", p.ra, p.dec);
}
void print_n2yo(tle_t tle){
    printf("https://n2yo.com/?s=%d\n", tle.norad_number);
}

double unix2jd(time_t unix_epoch){
    double epoch0_in_jd = 2440587.5;
    double epoch_f = unix_epoch;
    double days =  epoch_f / 86400.0 ;
    double jd = days + epoch0_in_jd;
    return jd;
    //might have issues with leap seconds? idk
}
double getjd(){
    time_t nowepoch = time(NULL);
    printf("Now: %d\n", nowepoch);
    return unix2jd(nowepoch);
}
void get_passes(){
    //np.rise = bisectSearchJD( startjd, endjd, tle, observer, BISECT_RISING); 
    //np.max = bisectSearchJD( np.rise.jd, endjd, tle, observer, BISECT_MAX); 
    //np.set = bisectSearchJD( np.max.jd,  endjd, tle, observer, BISECT_FALLING); 
}

char algo[33] = {0};

int main(int argc, char **argv){
    size_t num_satellites = 0;
    topo_pos_t obs = {0};
    latlon me = {0}; //need this because topo_pos_t is floats but maidenhead uses doubles
    me = maidenhead_to_latlon("FN41");
    obs.lat = me.lat;
    obs.lon = me.lon;

    sat_mem_t satellites[10] = {0};
#include "generated/sat_tles.c"
    printf("LL: %.8f\t%.8f\n", obs.lat, obs.lon);
    double start_jd = getjd();
    Pf(start_jd);
    float days = 0.57;
    for( int i = 0; i < 1; i++ ){
    /*for( int i = 0; i < num_satellites; i++ ){*/
        double jd_inc = 10.0/86400;
        int M = 0;
        strncpy(algo, "groundtruth", 25);
        for( double jd = start_jd; jd < start_jd+days; jd+=jd_inc) {
            /*printf("JD: %f\n", jd);*/
            /*for( int i = 0; i < 1; i++ ){*/
            /*int x = printf("\n%s\t", satellites[i].name );*/
            /*if( x < 8 ){*/
            /*printf("\t");*/
            /*}*/
            sat_pos_t current = calcSat( &satellites[i].tle, jd, obs);
            M+=1;
            /*print_pos(jd, obs, current );*/
            /*printf("\t\t\t\t"); print_n2yo(satellites[i].tle);*/
        }
        int N = 0;
        int mins_per_day = 1440;
        double pi = 3.141592653589793238462643383279502884197;
        double mean_motion_rpd = satellites[i].tle.xno * mins_per_day / (2*pi);
        Pf(mean_motion_rpd);
        double one_orbit_jd = 1/mean_motion_rpd;
        Pf(one_orbit_jd);
        int upward = false;
        strncpy(algo, "mmdev", 25);
        for( double jd = start_jd; jd < start_jd+days; jd+=jd_inc) {
            sat_pos_t pos1 = calcSat( &satellites[i].tle, jd, obs);
            sat_pos_t pos2 = calcSat( &satellites[i].tle, jd+1.0/86400, obs);
            N+=2;
            sat_pos_t current = pos1;
            double slope = pos2.elev - pos1.elev;
            slope *= 1000;
            /*Pf(slope);*/

            /*direction = 1;*/
            /*direction = -1;*/

            //bisect search to find a local maxima
            //jump ahead by one orbit
            //bisect search to find a local maxima
            //jump ahead by one orbit
            //repeat until orbit max is above horizon
            //find and store the rise, max, and fall times
            //go to top of loop
            if( slope >= 0 ){
                jd_inc = 10.0/86400;
                upward = true;
                /*sat_nextpass(&satellites[i].tle, start_jd, obs);*/
            } else if( slope < 0 && upward == true){
                //we hit the peak, now let's skip to the next peak if we can
                jd_inc = one_orbit_jd;
                upward = false;
            } else if( slope < 0 && upward == false){
                jd_inc = -60.0/86400;
            } else {
                jd_inc = 3000.0/86400; //guess from lookign at the graphs, should really be based on mean motion i guess
            }

        }
        printf("%s: M %d, N %d, N/M %f\n", satellites[i].name, M, N, ((float)N)/M);
        strncpy(algo, "old", 25);
        Pf(start_jd);
        sat_nextpass(&satellites[i].tle, start_jd, obs);
    }

        return 0;
    }
