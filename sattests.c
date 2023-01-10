
#include <assert.h>
#include <string.h>
#include "satellite.h"
#include "linux.h"
#include "sattests.h"

int _sat_algoN_read(int satid, char * algo);
void test_1(){
    jd_ts startjd = 2459943.332685; //NOT FOUND
    Pf(startjd);
    tle_t tle1 = {0};
    tle_t tle2 = {0};
    int r = loadtle("data/tles/historical/22321.txt",25544, &tle1); //ISS
    assert( r == 1 );
    r = loadtle("data/tles/historical/22363.txt",25544, &tle2); //ISS
    assert( r == 1 );
    topo_pos_t obs = {0};
    /*sat_pass_t pass1 = nextpass(obs, &tle1, startjd);*/
    /*sat_pass_t pass2 = nextpass(obs, &tle2, startjd);*/

    //https://stackoverflow.com/questions/1275484/good-plotting-library-for-c
    //https://stackoverflow.com/questions/63987163/how-do-you-draw-a-circle-using-gnuplot
    /*plot_pass(obs,&tle1,pass1);*/
    /*plot_pass(obs,&tle2,pass2);*/
    double frequency = 10e9;
    /*points1 = doppler_curve(pass1, frequency, 1000);*/
    //should use doppler_point somewhere in there? I guess?
    //doppler_curve should be unscaled by frequency
    /*points2 = doppler_curve(pass2, 1000);*/
}
void test_2(){
    jd_ts startjd = 2459943.332685;
    tle_t tle1 = {0};
    tle_t tle2 = {0};
    int r = loadtle("data/tles/historical/22321.txt",25544, &tle1); //ISS
    if( r != 1 ){ return; } //TODO ERROR
    r = loadtle("data/tles/historical/22363.txt",25544, &tle2); //ISS
    if( r != 1 ){ return; } //TODO ERROR
    topo_pos_t obs = {0};
    /*plot_ra_dec_alt();*/
    //3d plot in ECI xyz?
    //plot just the TLE data changing over time! That ought to be interesting enough
    //Plot the slope of the elevation curve too and see if there are any hints there
    //on how to do it faster
}
void test_3(){
}


soe passes_equalish( sat_pass_t p1, sat_pass_t p2, jd_ts precision_time, double precision_degrees){
}
errorcount compare_N_passes_between(tle_t tle, topo_pos_t obs, jd_ts jd, int N,
        nextpassfn np1, char * np1name,
        nextpassfn np2, char * np2name
        ){ 
    //take two implementations of nextpass and their names and compare their outputs
    int errors = 0;
    strncpy(algo, np1name, 25);
    sat_pass_t * passes1 = malloc( sizeof(sat_pass_t) * N);
    if( passes1 == NULL ){
        errors += 1;
        return errors;
    }
    errors += get_next_N_passes( tle, obs, jd, N, passes1, np1);
    assert( errors == 0 );

    strncpy(algo, np2name, 25);
    sat_pass_t * passes2 = malloc( sizeof(sat_pass_t) * N);
    if( passes2 == NULL ){
        errors += 1;
        return errors;
    }
    errors += get_next_N_passes( tle, obs, jd, N, passes2, np2);
    assert( errors == 0 ); 
    //these assertions are because we're returning
    //an errorcount, but the actual utility of this function is to compare
    //algorithms while they're under development - so we have to enforce
    //that there are zero errors prior to comparing the passes return.

    //given two nextpass algorithms for the same satellite, satellite data, observer position, and start time, 
    //we expect every pass to be found, to be found in the same order, etc.
    jd_ts precision_t = 1.0/86400; //+- 1second.
    jd_ts precision_d = .5; //+- 0.5 degree
    for( int i = 0; i < N; i++ ){
        errors += passes_equalish(passes1[i], passes2[i], precision_t, precision_d);
        //fprintf to stderr or a file here for debugging, maybe
        //or write a report with graphs out, would be nicer
    }
    return errors;
}
