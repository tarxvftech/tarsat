
#ifndef SATELLITE_H
#define SATELLITE_H

#define PI  3.141592653589793238462643383279
#define MINS_PER_DAY 1440
#include <norad.h>
#include <stdbool.h>
#include <stdint.h>
#include "misc.h"
#define radial_vel_offset_sec 5  //was time_offset
                                 //double time_offset = 5;
                                 //
typedef double jd_ts; //julian day timestamp
typedef int errorcode_satpass;
typedef int errorcode_satpos;

typedef struct {
    float ra;
    float dec;
} equa_pos_t;

typedef struct {
    char name[16];
    float ra;
    float dec;
    float mag;
} star_t;


typedef struct {
    float lat;   //degrees
    float lon;   //degrees
    float alt;   //meters MSL
} topo_pos_t;

typedef struct {
    jd_ts jd;  //julian day/time stamp
               //could also commit fully to j2k throughout, if i were so inclined...
    float  az;   //degrees
    float  elev; //degrees
    float  dist; //meters
} sat_azel_t;

typedef struct {
    //doesn't include observer because it can only be generated with data from observer, so calcSat caller keeps track of it
    jd_ts jd;  //time of position

    double az;   //azimuth from observer, degrees (0/360 north, 180 south, etc
    double elev; //elevation, degrees (0 is horizon, 90 is straight up normal to the ground)
    double dist; //distance, meters (range)

    double ra;  //right ascension, degrees
    double dec; //declination, degrees

    int satid; //catalog number of the satellite
    int err;     // ==0 if all is okay, so make sure to check this
} sat_pos_t;

typedef struct {
    sat_pos_t rise; //satellite position when it comes over horizon (elev - to +)
    sat_pos_t max;  //when it's at max elevation
    sat_pos_t set;  //when it goes over horizon (elev + to -)
    int err;     // ==0 if all is okay, so make sure to check this
} sat_pass_t;

typedef struct {
    char       name[16];
    tle_t      tle;
    sat_pos_t  current;  
    sat_pass_t nextpass; 
} sat_mem_t;

extern sat_mem_t satellites[];
extern int num_satellites;
extern int satellites_initialized;
extern const star_t stars[];
extern int num_stars;


//But first, we have to talk about units.
//_rad and _deg suffixes indicate non-standard units where necessary, or where it might just be ambiguous
//
//RA (Right Ascension) is assumed to be in hours.
//      You can convert RA_hrs to RA_deg with the HRS2DEG() macro (which is just *15)
//DEC (Declination) is assumed to be in degrees
//latitude and longitude are assumed to be in degrees
//

topo_pos_t getObserverPosition();
double jd_to_j2k(double jd);
double local_sidereal_degrees(double j2k, double longitude);
double hour_angle_degrees(double lst_deg, double ra_hrs);
void ra_dec_to_az_alt(double jd,
        double latitude, double longitude,
        double ra_hrs, double dec_deg,
        double* az_out_deg, double* alt_out_deg);

sat_pos_t  calcSat(tle_t * tle, topo_pos_t observer_degrees, jd_ts time_jd);


/*
   sat_pass_t sat_nextpass (
//sat in question
tle_t * tle,
//start time
double startjd,
//observer location, in degrees latitude and longitude, and altitude in meters
topo_pos_t observer
); //TODO convert to a satellite pass type
sat_pos_t bisectSearchJD(
double startjd,
double endjd,
tle_t * tle,
topo_pos_t observer,
enum BisectSearchType searchtype
);
*/


//stuff for the "asteroids" "game"
typedef struct {
    float rot; //radians!
    float spdx;
    float spdy;
    float x;
    float y;
} game_obj_2d_t;

void game_move(game_obj_2d_t* o, unsigned long long td);
void game_addvel(game_obj_2d_t* o, float vel, float rot);
void game_obj_init(game_obj_2d_t* o);
void game_obj_screenwrap(game_obj_2d_t* o);
#define DEG(x) ( (x)*180/PI )
#define RAD(x) ( (x)*PI/180 )
#define HRS2DEG(x) ((x)*15)

typedef enum Direction
{
    BOTH  = 0, //endjd = startjd + .5 orbit; startjd -= .5 orbit
    RIGHT = 1, //endjd = startjd + 1 orbit
    LEFT  = 2, //endjd = startjd; startjd = startjd - 1 orbit
} direction_t;

typedef enum SearchType
{
    SEARCH_FAILURE =0,      //for when returning values, indicates not found
    SEARCH_RISING = 1,      //crossing zero / low to high 
    SEARCH_FALLING= 1<<1,   //crossing zero \ high to low 
    SEARCH_MAX    = 1<<2,   //local maximum ^
    SEARCH_MIN    = 1<<3,   //local minimum v
    SEARCH_PASSMAX= 1<<4,   //highest point above the horizon (implies a rising before and a falling after), plus passmaxes are strictly a subset of the max points
} search_t;
typedef struct foundfeature {
    sat_pos_t pos; //pos has a jd in it
    search_t feature;
} foundfeature_t;

foundfeature_t search_hillclimb(
        topo_pos_t obs,
        tle_t* tle,

        jd_ts stepsize,
        jd_ts startjd,
        jd_ts maxjd
        );
foundfeature_t search_simple(
        topo_pos_t obs,
        tle_t* tle,
        search_t searchtype,

        direction_t dir,
        jd_ts stepsize,
        jd_ts startjd,
        jd_ts maxjd
        );

foundfeature_t bisectSearchJD(
        topo_pos_t obs,
        tle_t* tle,
        search_t searchtype,
        jd_ts precision,
        jd_ts startjd,
        jd_ts maxjd
        );
typedef sat_pass_t(*nextpassfn)(tle_t tle, topo_pos_t obs, jd_ts startjd);
errorcount get_next_N_passes(tle_t tle, topo_pos_t obs, jd_ts startjd, int N, sat_pass_t * passes, nextpassfn np);

uint32_t makelut();
uint32_t lutkey(search_t st, int slope, int elev);
uint32_t lutval(int st, int slope, int elev, int dir);

#endif
