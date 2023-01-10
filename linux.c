#include <time.h>
#include "misc.h"
#include "satellite.h"

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
sof loadtle(char * filename, int catalog_number, tle_t * tle){
    //requires format like amsat TLE format e.g.
    //friendlyname1
    //1 25544U ...
    //2 25544  ...
    //friendlyname2 ...
    FILE* fd = fopen(filename, "r");

    char friendlyname[32] = {0};
    char line1[96] = {0};
    char line2[96] = {0};

    if(!fd) { return -1; }
    while (1){
        char * r = NULL;
        r = fgets(friendlyname, 32, fd);  if( r == NULL ){ return 0; }
        r = fgets(line1, 96, fd);         if( r == NULL ){ return 0; }
        r = fgets(line2, 96, fd);         if( r == NULL ){ return 0; }

        parse_elements( line1, line2, tle);
        if( tle->norad_number == catalog_number ){
            return true;
        } 
    }
    return false;
}
