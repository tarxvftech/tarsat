#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include <observe.h>

#include "satellite.h"
#include <maidenhead.h>
#include <watdefs.h>
#include <afuncs.h>


double jd_to_j2k(double jd)
{
    return jd - 2451545.0;
}
double wrap( double in, double max)
{
    //TODO replace with fmod i think
    while( in < 0 ) {
        in += max;
    }
    while( in >= max ) {
        in -= max;
    }
    return in;
}
double UT_dechrs_from_jd(double jd){
    double UT = 24 * (jd - floor(jd) + .5);
    return UT;
}
double local_sidereal_degrees(double j2k, double longitude )
{
    double UT = UT_dechrs_from_jd(j2k); //hours, naturally
    /*printf("J2K: %.6f \n", j2k);*/
    /*printf("UT: %.6f \n", UT);*/
    double degrees_rotation_per_day = .985647;
    double gmst_j2k_correction = 100.46;
    double local_sidereal_time = gmst_j2k_correction + degrees_rotation_per_day * j2k + longitude + 15 * UT;
    local_sidereal_time = wrap(local_sidereal_time, 360);
    return local_sidereal_time;
}
double hour_angle_degrees( double lst_deg, double ra_hrs)
{
    double ha = lst_deg - HRS2DEG(ra_hrs);
    ha = wrap( ha, 360 );
    //note that hour angle is in range -PI,PI in liblunar...and we're 0,360
    return ha; //hour angle is in degrees
}
void ra_dec_to_alt_az( const double hr_ang, const double dec, double *alt, double *az, const double lat)
{
   double temp, cos_lat = cos( lat);

   *alt = asine( sin( lat) * sin( dec) + cos_lat * cos( dec) * cos( hr_ang));
   if( cos_lat < .00001)         /* polar case */
      *az = hr_ang;
   else
      {
      temp = (sin( dec) - sin( *alt) * sin( lat)) / (cos( *alt) * cos_lat);
      temp = PI - acose( temp);
      *az = ((sin( hr_ang) < 0.) ? temp : -temp);
      }
}
void ra_dec_to_az_alt(double jd,
                      double latitude, double longitude,
                      double ra_hrs, double dec_deg,
                      double * az_out_deg, double * alt_out_deg)
{
    if( 1 ){
        double j2k = jd_to_j2k(jd);
        double lst = local_sidereal_degrees(j2k, longitude);
        /*printf("lst %.4f\n", lst);*/
        double ha = RAD(hour_angle_degrees(lst, ra_hrs )); //radians
        /*printf("HA mine = %f\n", DEG(ha));*/

        double latrad = RAD(latitude);
        double dec = RAD(dec_deg);

        double alt = asin( sin(latrad) * sin(dec) + cos(latrad)*cos(dec)*cos(ha));
        double A = acos(
                (sin(dec) - sin(alt)*sin(latrad))
                /
                (cos(alt)*cos(latrad))
                );
        double az = sin(ha) < 0? A: 2*PI-A;

        az = wrap(az, 2 * PI);

        *az_out_deg = DEG(az);
        *alt_out_deg = DEG(alt);
        return;
    } else {
        /*DPT radec = {RAD(dec_deg), RAD(HRS2DEG(ra_hrs))}; // from looking at precess.cpp, both expected to be in suitable form for sin,cos etc*/

        //full_ra_dec_to_alt_az
        //DPT ra_dec comes in
        //gets pushed to precess_pt, which spits out a DPT loc_at_epoch
        //DPT ra_dec is renamed DPT ipt (in point)
        //temp_ipt[0] = -ipt->x  <-RA
        //temp_ipt[1] = ipt->y <- DEC
        //precess_ra_dec:
        //      where temp_ipt[0] is "old_ra"
        //      so:
        //loc_at_epoch->x = -temp_opt[0] <- RA
        //loc_at_epoch->y = -temp_opt[1] <- DEC
        //when loc_at_epoch is handed to ra_dec_to_alt_az, ra_dec.y is DEC!
        //so DPT radec is correct {ra,dec}
        //
        DPT radec = {RAD(HRS2DEG(ra_hrs)), RAD(dec_deg)};
        // from looking at precess.cpp, both expected to be in suitable form for sin,cos etc
            //so needs to be radians
        /*DPT latlon = {RAD(latitude), RAD(longitude)};*/
        DPT latlon = {RAD(longitude), RAD(latitude)}; //x,y; and x is _longitude_! BAD NAMING!
        //BEWARE -- lunar library stores data in DPT
        //BUT!
        //DPT alt_az = {double x, y}
        //and they consistently use it y, then x!
        //which means when you use it, latlon is stored x=lon, y=lat
        //alt_az is stored alt=y, az=x
        DPT altaz = {0};
        //by convention, alt should be (-PI,PI), az should be (0,2*PI)
        //except there is some evidence in liblunar that az is likely to be (-PI,PI)
        //hence the +PI in there at the end
        double j2k = jd_to_j2k(jd);
        double ha;
        full_ra_dec_to_alt_az(
                &radec, // DPT *ra_dec,
                &altaz, //DPT *alt_az,
                NULL, //DPT *loc_epoch,
                &latlon, //DPT *latlon,
                jd, //double jd_utc,
                &ha //double *hr_ang
                );
        *alt_out_deg = DEG(altaz.y);
        *az_out_deg  = DEG(altaz.x +PI);
        /*printf("ha = %f\n", DEG(ha+PI));*/
    }
}

extern char algo[];
sat_pos_t  calcSat( tle_t * tle, double time_jd, topo_pos_t observer_degrees)
{
    topo_pos_t obs = { //observer position, but in radians
        RAD(observer_degrees.lat),
        RAD(observer_degrees.lon),
        observer_degrees.alt
    };
    int is_deep = select_ephemeris( tle );
    int ephem = 1;       /* default to SGP4 */
    double sat_params[N_SAT_PARAMS], observer_loc[3];
    double rho_sin_phi;
    double rho_cos_phi;

    double ra;
    double dec;
    double dist_to_satellite;

    double t_since;
    double pos[3];
    //double vel[3]; //pass in place of the NULL after pos to get velocities
    int err_val = 0;

    //remember to put lat and lon in rad by this point
    earth_lat_alt_to_parallax( obs.lat, obs.alt, &rho_cos_phi, &rho_sin_phi);
    observer_cartesian_coords( time_jd, obs.lon, rho_cos_phi, rho_sin_phi, observer_loc);
    if( is_deep && (ephem == 1 || ephem == 2)) {
        ephem += 2;
    }
    if( !is_deep && (ephem == 3 || ephem == 4)) {
        ephem -= 2;
    }

    t_since = (time_jd - tle->epoch) * 1440;
    switch( ephem) {
    case 0:
        SGP_init( sat_params, tle);
        err_val = SGP( t_since, tle, sat_params, pos, NULL);
        break;
    case 1:
        SGP4_init( sat_params, tle);
        err_val = SGP4( t_since, tle, sat_params, pos, NULL);
        break;
    case 2:
        SGP8_init( sat_params, tle);
        err_val = SGP8( t_since, tle, sat_params, pos, NULL);
        break;
    case 3:
        SDP4_init( sat_params, tle);
        err_val = SDP4( t_since, tle, sat_params, pos, NULL);
        break;
    case 4:
        SDP8_init( sat_params, tle);
        err_val = SDP8( t_since, tle, sat_params, pos, NULL);
        break;
    default:
        //printf( "? How did we get here? ephem = %d\n", ephem);
        err_val = 0;
        break;
    }
    if( err_val ) {
        /*printf( "Ephemeris error %d\n", err_val);*/
    }
    get_satellite_ra_dec_delta( observer_loc, pos, &ra, &dec, &dist_to_satellite);
    epoch_of_date_to_j2000( time_jd, &ra, &dec);
    double az = 0;
    double elev = 0;
    ra_dec_to_az_alt(time_jd, DEG(obs.lat), DEG(obs.lon), DEG(ra)/15, DEG(dec), &az, &elev);
    //printf("POS: %.4f,%.4f,%.4f\n", pos[0], pos[1], pos[2] );
    /*printf("VEL: %.4f,%.4f,%.4f\n", vel[0], vel[1], vel[2] );*/
    sat_pos_t ret;
    ret.az = az;
    ret.elev = elev;
    ret.ra = DEG(ra)/15;
    ret.dec = DEG(dec);
    ret.dist = dist_to_satellite;
    ret.jd = time_jd;
    ret.satid = tle->norad_number;
    ret.err = err_val;
    char fn[64] = {0};
    snprintf(fn, 64, "%d_%s.csv", ret.satid, algo);
    FILE * fd = fopen(fn, "a");
    fprintf(fd,"%f,%f,%f,%f,%f,%f\n", time_jd, ret.az, ret.elev, ret.dist, ret.ra, ret.dec);
    fclose(fd);
    snprintf(fn, 64, "%d_%s_N.csv", ret.satid, algo);

    fd = fopen(fn, "r");
    int n = 0;
    if( fd ){
       fscanf(fd, "%d", &n);
       fclose(fd);
    } 

    fd = fopen(fn, "w");
    n += 1;
    if( fd ){
       fprintf(fd, "%d\n", n);
    } 
    fclose(fd);
    return ret;
}




sat_pos_t bisectSearchJD(
        double startjd,
        double endjd,
        tle_t * tle,
        topo_pos_t observer,
        enum BisectSearchType searchtype
        ){
    //find zero crossings +,-, and local maxima --and-minima-- (minima not implemented)
    //makes me wish for currying ... https://stackoverflow.com/a/66865994
    //
    //hard requirement: The desired feature must be between startjd and endjd, and there should be only one or you'll get unexpected answers
    //similarly, understand that these tests are based on slope - and too long a search area will have multiple areas of similar slope 
    //which will break the assumptions for bisect search
    //keep your search area to less than 3/4 (or better, 1/2) "wavelength" for periodic stuff
    //TODO helper fn should handle this correctly 
    //
    //Ideally you'll set startjd and endjd knowing what the orbit period is and having a rough idea of where you can expect a pass to happen 
    // - parts of which you can figure out from the bare TLE data (mean motion) and a sample or two.
    //https://space.stackexchange.com/a/23451
    //TODO: finish this fn and turn nextpass into a helper function to do just that
/*
 
        Curve measured at midpoint of values (*)
        |            Negative elevation, negative slope   |
        |         /^\                                     |
        |       /    \                                    |
        |   - / ----- \ ----------*---------------------  | 
        |   /          \                                  |
        | /              \                                |
        |                  \                              |
        |   negative elevation, positive slope            |
        |                                 /^\             |
        |                               /    \            |
        |    ---------------------*-- / ----- \ --------  |
        |                           /          \          |
        |                         /              \        |
        |                                          \      |
        |           Positive elevation, negative slope    | (and pos elev with positive slope is obvious from this)
        |                     /^\                         |
        |                   /    \                        |
        |                 /       \                       |
        |    ---------- / --------*\ -----------------    |
        |             /              \                    |
        |                              \                  |


        elevation sign, slope sign (direction to feature)
        (= means no measurable slope or zero for elevation)
Max:
        -,- (left)
        -,+ (right)
        +,- (left)
        +,+ (right)
        +,= (return, we found it!) (nah we'll just say left or something or special case it in base)
Rising: (dawning)
        -,- (left)
        -,+ (right)
        +,- (left)
        +,+ (left)
Falling: (setting)
        -,- (left)
        -,+ (right)
        +,- (right)
        +,+ (right)
enum BisectSearchType
{
    BISECT_RISING = 0,
    BISECT_FALLING,
    BISECT_MAX,
    BISECT_MIN
};
searchtype = 00,01,10,(11) -> rising,falling,max, min (min unimplemented)
elevation = elevation > 0?1:0;
slope = slope > 0?1:0;
key = (st&3) << 2 | (elevation&1) << 1 | (slope&1) << 0 
goRight (0 for left, 1 for right)
bool lut[16] = {
        //least significant bit
    0, //0b0000 Rising (0b00), negative elevation, negative slope, go left
    1, //0b0001 negative elevation, positive slope, go right 
    0, //0b0010 positive elevation, negative slope, go left
    0, //0b0011 ++,left
 
    0, //0b0100 Falling (0b01), negative elevation, negative slope, go left
    1, //0b0101 -+,right
    1, //0b0110 +-,right
    1, //0b0111 ++,right

    0, //0b1000 Max (0b10), --,left
    1, //0b1001 -+, right
    0, //0b1010 +-, left
    1, //0b1011 ++, right

    0, //0b1100 min (0b11) unimplemented so not filled out
    0, //0b1101
    0, //0b1110
    0, //0b1111 
       //most significant bit
};
That same table as a 16bit short:
lut = 0xae2
and we can pull out the bit in question by bit shifting appropriately and masking to 1 bit
goRight = (lut >> key) & 1; (if(goRight){...}else{...}

*/
    //printf("bisect search: startjd: %f, endjd:%f\n", startjd, endjd);

    double jd_1s = 1.0 / 86400; //1 second in decimal days
    double precision = 1*jd_1s;

    double width = endjd-startjd; 
    //endjd must be beyond startjd
    if( width < 0 ){
        //TODO error! out of order.
        sat_pos_t x;
        x.err = -1;
        return x;
    }

    double midpoint = startjd + width/2;
    if( width <= precision ){
        //if we're at the feature - and we should be able to sanity check that, then 
        sat_pos_t l = calcSat( tle, midpoint-5*jd_1s, observer);
        sat_pos_t m = calcSat( tle, midpoint, observer);
        sat_pos_t r = calcSat( tle, midpoint+5*jd_1s, observer);
        switch(searchtype){
            case BISECT_RISING:
                if( !( l.elev < 0 && l.elev < m.elev && m.elev < r.elev && 0 < r.elev )){
                    m.err = -1;
                } else {
                    m.err = 0;
                }
                break;
            case BISECT_FALLING:
                if( !( l.elev > 0 && l.elev > m.elev && m.elev > r.elev && 0 > r.elev )){
                    m.err = -1;
                } else {
                    m.err = 0;
                }
                break;
            case BISECT_MAX:
                if( !( l.elev < m.elev && m.elev > r.elev )){
                    m.err = -1;
                } else {
                    m.err = 0;
                }
                break;
        }
        //printf("Returning with %s at elev %.1f\n", m.err == 0? "SUCCESS":"FAILURE", m.elev);
        return m;
    } else {
        sat_pos_t m1 = calcSat( tle, midpoint, observer);
        sat_pos_t m2 = calcSat( tle, midpoint+1*jd_1s, observer); //pop right a little so we can measure the change in elevation
        double mslope_elev = m2.elev - m1.elev; //positive mslope_elev represents increasing values left-to-right
        char elevation = m1.elev >= 0?1:0; //positive is truthy
        char slope = mslope_elev >= 0?1:0; //negative is falsy
        char key = (searchtype&3) << 2 | (elevation&1) << 1 | (slope&1) << 0;
        uint16_t lut = 0xae2; //direction lookup table
        bool goRight = (lut >> key) & 1;
        //printf("\tkey: 0x%x goRight: %d\n", key, goRight);
        
        if( goRight ){
            return bisectSearchJD(midpoint, endjd, tle, observer, searchtype);
        } else { 
            return bisectSearchJD(startjd, midpoint, tle, observer, searchtype);
        }
    }
}
sat_pass_t sat_nextpass(
    //sat in question
    tle_t * tle,
    //start time
    double startjd,
    //observer location, in degrees latitude and longitude, and altitude in meters
    topo_pos_t observer
)
{
    //determine sampling from TLE ?
    //getting a list of passes involves successive calls to this function, incrementing startjd
    sat_pass_t np; //this represent the entirety of the next pass, we'll return this
                   //np.rise,max,set, np.rise.jd,az,elev,dist
                   
    //first, some quick notes: mean motion is the average (hence
    //'mean') number of orbits per day. "Day" is complicated but close enough to what you expect to not matter here.
    //We want to use mean motion to estimate how far to search for satellite orbits.
    //Revolutions per day is usually abbreviated as "N".
    //When expressed in radians per minute, that same speed is written as "Xno" (or 2*PI*N / 1440)
    //Bill Gray's sat_code parses the TLE data (orbits per day) into this Xno notation and the associated radians per minute units.
    //I just want decimal days per orbit, though - so I'm going to undo that a bit
    //to get back to orbits per day, and then invert it to get days per orbit, like so:
    double orbit_period_jd = 1. / (1440*tle->xno/(2*PI)); 

    double endjd = startjd + 1*orbit_period_jd;  
    //orbit period does not guarantee a pass over us within that time! But if there's a pass, we only want one. 
    //To find passes we'll have to offset and search those spaces just like we will this first time
    double jd_1s = 1.0 / 86400; //1 second in decimal days
    int max_orbits_search = 10;
    int i = 0;
    do {
        //printf("searching ... %d\n",i);
        np.rise = bisectSearchJD( startjd, endjd, tle, observer, BISECT_RISING); 
        startjd += 3*orbit_period_jd;
        endjd = startjd + 3*orbit_period_jd;
        i++;
    } while ( np.rise.err != 0 && i < max_orbits_search );
    if(np.rise.err != 0 || i > max_orbits_search){
        //ERR!
        //TODO handle errors when we can't find a pass, maybe just by continuing to search
        //still need to report negative results
        //printf("ERROR, no pass found\n");
        np.err = np.rise.err;
        return np;
    } else {
        //printf("SUCCESS, pass found\n");
        /*printf("JD: %f, azel %.0f,%.0f\n", np.rise.jd, np.rise.az, np.rise.elev);*/
        np.err = 0;
    }
    endjd = np.rise.jd + orbit_period_jd / 4;
    np.max = bisectSearchJD( np.rise.jd, endjd, tle, observer, BISECT_MAX); 
    np.set = bisectSearchJD( np.max.jd,  endjd, tle, observer, BISECT_FALLING); 
                                        
    return np;
}


void init_sat_global(){
    //Hardcoded TLEs must be updated at compile time
    //This also implies updating TLEs requires a firmware update.
    //This is a temporary measure until nvm is ready.
//#include "generated/sat_tles.c"
    satellites_initialized = 0; //only init_sat_global is allowed to increment from -1 to 0

}
sat_mem_t satellites[10] = {0};
int num_satellites = 0;
int satellites_initialized = -1;
const star_t stars[] = {
// { name, right ascension, declination, magnitude }
    {"Polaris",         2+32/60,  89.3,  1.99},
    {"Sirius",          6+45/40, -16.7, -1.46},
    {"Canopus",         6+24/40, -52.7, -0.73},
    {"A. Centauri",    14+40/40, -60.8, -0.29},
    {"Arcturus",       14+16/40,  19.2, -0.05},
    {"Vega",           18+37/40,  38.8,  0.03},
    {"Capella",         5+17/40,  46.0,  0.07},
    {"Rigel",           5+15/40,  -8.2,  0.15},
    {"Procyon",         7+39/40,   5.2,  0.36},
    {"Achernar",        1+38/40, -57.2,  0.45},
    {"Betelgeuse",      5+55/40,   7.4,  0.55},
    {"Hadar",          14+ 4/40, -60.4,  0.61},
    {"Altair",         19+51/40,   8.9,  0.77},
    {"Acrux",          12+27/40, -63.1,  0.79},
    {"Aldebaran",       4+36/40,  16.5,  0.86},
    {"Antares",        16+29/40, -26.4,  0.95},
    {"Spica",          13+25/40, -11.2,  0.97},
    {"Pollux",          7+45/40,  28.0,  1.14},
    {"Fomalhaut",      22+58/40, -29.6,  1.15},
    {"Deneb",          20+41/40,  45.3,  1.24},
    {"Mimosa",         12+48/40, -59.7,  1.26},
    {"Regulus",        10+ 8/40,  12.0,  1.36},
    {"Adhara",          6+59/40, -29.0,  1.50},
    {"Castor",          7+35/40,  31.9,  1.58},
    {"Shaula",         17+34/40, -37.1,  1.62},
    {"Gacrux",         12+31/40, -57.1,  1.63},
    {"Bellatrix",       5+25/60,   6.3,  1.64},
    {"Elnath",          5+26/60,  28.6,  1.66},
    {"Miaplacidus",     9+13/60, -69.7,  1.67},
    {"Alnilam",         5+36/60,  -1.2,  1.69},
    {"Alnair",         22+ 8/60, -47.0,  1.74},
    {"Alnitak",         5+41/60,  -1.9,  1.75},
    {"Alioth",         12+54/60,  56.0,  1.77},
    {"Mirfak",          3+24/60,  49.9,  1.80},
    {"Dubhe",          11+ 4/60,  61.8,  1.80},
    {"Regor",           8+10/60, -47.3,  1.81},
    {"Wezen",           7+ 8/60, -26.4,  1.83},
    {"Kaus Australis", 18+24/60, -34.4,  1.84},
    {"Alkaid",         13+48/60,  49.3,  1.86},
    {"Sargas",         17+37/60, -43.0,  1.86},
    {"Avior",           8+23/60, -59.5,  1.87},
    {"Menkalinan",      6+ 0/60,  44.9,  1.90},
    {"Atria",          16+49/60, -69.0,  1.92},
    {"Alhena",          6+38/60,  16.4,  1.93},
    {"Peacock",        20+26/60, -56.7,  1.93},
    {"Koo She",         8+45/60, -54.7,  1.95},
    {"Mirzam",          6+23/60, -18.0,  1.98},
    {"Alphard",         9+28/60,  -8.7,  1.98},
    {"Polaris",         2+32/60,  89.3,  1.99},
    {"Algieba",        10+20/60,  19.8,  2.00},
    {"Hamal",           2+ 7/60,  23.5,  2.01},
    {"Diphda",          0+44/60, -18.0,  2.04},
    {"Nunki",          18+55/60, -26.3,  2.05},
    {"Menkent",        14+ 7/60, -36.4,  2.06},
    {"Alpheratz",       0+ 8/60,  29.1,  2.07},
    {"Mirach",          1+10/60,  35.6,  2.07},
    {"Saiph",           5+48/60,  -9.7,  2.07},
    {"Kochab",         14+51/60,  74.2,  2.07},
    {"Al Dhanab",      22+43/60, -46.9,  2.07},
    {"Rasalhague",     17+35/60,  12.6,  2.08},
    {"Algol",           3+ 8/60,  41.0,  2.09},
    {"Almach",          2+ 4/60,  42.3,  2.10},
    {"Denebola",       11+49/60,  14.6,  2.14},
    {"Cih",             0+57/60,  60.7,  2.15},
    {"Muhlifain",      12+42/60, -49.0,  2.20},
    {"Naos",            8+ 4/60, -40.0,  2.21},
    {"Aspidiske",       9+17/60, -59.3,  2.21},
    {"Alphecca",       15+35/60,  26.7,  2.22},
    {"Suhail",          9+ 8/60, -43.4,  2.23},
    {"Mizar",          13+24/60,  54.9,  2.23},
    {"Sadr",           20+22/60,  40.3,  2.23},
//#include "data/stars/stars.table" //~3150 more stars
//(not individually named, constellation abbreviations for the names)
//ahahaha yeah no, uses way too much space
//and WAY too long to draw, damn!
};
int num_stars = sizeof( stars ) / sizeof( star_t );
