#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include <observe.h>

#include "satellite.h"
#include <maidenhead.h>
#include <watdefs.h>
#include <afuncs.h>
#include <assert.h>

int _sat_algoN_write(int satid, char * algo, int n);
int _sat_algoN_read(int satid, char * algo);


double jd_to_j2k(double jd)
{
    return jd - 2451545.0;
}
double wrap(double in, double max)
{
    //TODO replace with fmod i think
    while(in < 0) {
        in += max;
    }
    while(in >= max) {
        in -= max;
    }
    return in;
}
double UT_dechrs_from_jd(double jd)
{
    double UT = 24 * (jd - floor(jd) + .5);
    return UT;
}
double local_sidereal_degrees(double j2k, double longitude)
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
double hour_angle_degrees(double lst_deg, double ra_hrs)
{
    double ha = lst_deg - HRS2DEG(ra_hrs);
    ha = wrap(ha, 360);
    //note that hour angle is in range -PI,PI in liblunar...and we're 0,360
    return ha; //hour angle is in degrees
}
void ra_dec_to_alt_az(const double hr_ang, const double dec, double* alt, double* az, const double lat)
{
    double temp, cos_lat = cos(lat);

    *alt = asine(sin(lat) * sin(dec) + cos_lat * cos(dec) * cos(hr_ang));
    if(cos_lat < .00001) {        /* polar case */
        *az = hr_ang;
    } else {
        temp = (sin(dec) - sin(*alt) * sin(lat)) / (cos(*alt) * cos_lat);
        temp = PI - acose(temp);
        *az = ((sin(hr_ang) < 0.) ? temp : -temp);
    }
}
void ra_dec_to_az_alt(double jd,
        double latitude, double longitude,
        double ra_hrs, double dec_deg,
        double* az_out_deg, double* alt_out_deg)
{
    if(1) {
        double j2k = jd_to_j2k(jd);
        double lst = local_sidereal_degrees(j2k, longitude);
        /*printf("lst %.4f\n", lst);*/
        double ha = RAD(hour_angle_degrees(lst, ra_hrs));  //radians
        /*printf("HA mine = %f\n", DEG(ha));*/

        double latrad = RAD(latitude);
        double dec = RAD(dec_deg);

        double alt = asin(sin(latrad) * sin(dec) + cos(latrad)*cos(dec)*cos(ha));
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
sat_pos_t  calcSat(tle_t* tle, topo_pos_t observer_degrees, double time_jd)
{
    topo_pos_t obs = { //observer position, but in radians
        RAD(observer_degrees.lat),
        RAD(observer_degrees.lon),
        observer_degrees.alt
    };
    int is_deep = select_ephemeris(tle);
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
    earth_lat_alt_to_parallax(obs.lat, obs.alt, &rho_cos_phi, &rho_sin_phi);
    observer_cartesian_coords(time_jd, obs.lon, rho_cos_phi, rho_sin_phi, observer_loc);
    if(is_deep && (ephem == 1 || ephem == 2)) {
        ephem += 2;
    }
    if(!is_deep && (ephem == 3 || ephem == 4)) {
        ephem -= 2;
    }

    t_since = (time_jd - tle->epoch) * 1440;
    switch(ephem) {
        case 0:
            SGP_init(sat_params, tle);
            err_val = SGP(t_since, tle, sat_params, pos, NULL);
            break;
        case 1:
            SGP4_init(sat_params, tle);
            err_val = SGP4(t_since, tle, sat_params, pos, NULL);
            break;
        case 2:
            SGP8_init(sat_params, tle);
            err_val = SGP8(t_since, tle, sat_params, pos, NULL);
            break;
        case 3:
            SDP4_init(sat_params, tle);
            err_val = SDP4(t_since, tle, sat_params, pos, NULL);
            break;
        case 4:
            SDP8_init(sat_params, tle);
            err_val = SDP8(t_since, tle, sat_params, pos, NULL);
            break;
        default:
            //printf( "? How did we get here? ephem = %d\n", ephem);
            err_val = 0;
            break;
    }
    if(err_val) {
        /*printf( "Ephemeris error %d\n", err_val);*/
    }
    get_satellite_ra_dec_delta(observer_loc, pos, &ra, &dec, &dist_to_satellite);
    epoch_of_date_to_j2000(time_jd, &ra, &dec);
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
    FILE* fd = fopen(fn, "a");
    fprintf(fd,"%f,%f,%f,%f,%f,%f\n", time_jd, ret.az, ret.elev, ret.dist, ret.ra, ret.dec);
    fclose(fd);

    int n = _sat_algoN_read(ret.satid, algo);
    n += 1;
    _sat_algoN_write(ret.satid, algo, n);

    return ret;
}

int _sat_algoN_write(int satid, char * algo, int n){
    char fn[64] = {0};
    snprintf(fn, 64, "%d_%s_N.csv", satid, algo);
    FILE * fd = fopen(fn, "w");
    if(fd) {
        fprintf(fd, "%d\n", n);
    }
    fclose(fd);
}
int _sat_algoN_read(int satid, char * algo){
    char fn[64] = {0};
    snprintf(fn, 64, "%d_%s_N.csv", satid, algo);
    FILE * fd = fopen(fn, "r");
    int n = 0;
    if(fd) {
        fscanf(fd, "%d", &n);
        fclose(fd);
    }
    return n;
}


foundfeature_t search_hillclimb(
        tle_t* tle,
        topo_pos_t obs,

        jd_ts stepsize,
        jd_ts startjd,
        jd_ts maxjd
        ){
    ;
}
foundfeature_t search_simple(
        tle_t* tle,
        topo_pos_t obs,
        search_t searchtype,
        direction_t dir,

        jd_ts stepsize,
        jd_ts startjd,
        jd_ts maxjd
        )
{
    foundfeature_t found = {0};
    double jd = startjd;
    sat_pos_t m1 = calcSat(tle, obs, jd);
    double lastslope = 0;
    do {
        sat_pos_t m2 = calcSat(tle, obs, jd);
        double elevslope = m2.elev - m1.elev;
        if( searchtype & SEARCH_MIN ){
            if(lastslope < 0 && elevslope > 0) {
                found.pos = m1;
                found.feature = SEARCH_MIN;
                return found;
            } 
        }
        if( searchtype & SEARCH_PASSMAX){
            if(lastslope > 0 && elevslope < 0 && m2.elev > 0) {
                found.pos = m1;
                found.feature = SEARCH_PASSMAX | SEARCH_MAX;
                return found;
            }
        }
        if( searchtype & SEARCH_MAX){
            if(lastslope > 0 && elevslope < 0) {
                found.pos = m1;
                found.feature = SEARCH_MAX;
                return found;
            }
        }
        m1 = m2;
        lastslope = elevslope;
        if(dir == LEFT) {
            jd-=stepsize;
        } else {
            jd+=stepsize;
        }
    } while(jd < maxjd);
    found.pos = m1;
    found.feature = SEARCH_FAILURE;
    return found;
}

uint32_t lutkey(search_t st, int slope, int elev){
    return ( 
            ((st&15)<<2) 
            | 
            ((slope&1) << 1) 
            |
            ((elev&1) << 0) 
           ) ;
}
uint32_t lutval(int st, int slope, int elev, int dir){
    return (dir << 
            lutkey(st,slope,elev)
           );
}
uint32_t makelut(){
    uint32_t lut = 0;
#define NEG (0)
#define POS (1)

    lut |= lutval(SEARCH_RISING, NEG, NEG, LEFT);
    lut |= lutval(SEARCH_RISING, NEG, POS, LEFT);
    lut |= lutval(SEARCH_RISING, POS, NEG, RIGHT);
    lut |= lutval(SEARCH_RISING, POS, POS, LEFT);

    lut |= lutval(SEARCH_FALLING, NEG, NEG, LEFT);
    lut |= lutval(SEARCH_FALLING, NEG, POS, RIGHT);
    lut |= lutval(SEARCH_FALLING, POS, NEG, RIGHT);
    lut |= lutval(SEARCH_FALLING, POS, POS, RIGHT);

    lut |= lutval(SEARCH_MAX, NEG, NEG, LEFT);
    lut |= lutval(SEARCH_MAX, NEG, POS, LEFT);
    lut |= lutval(SEARCH_MAX, POS, NEG, RIGHT);
    lut |= lutval(SEARCH_MAX, POS, POS, RIGHT);

    lut |= lutval(SEARCH_MIN, NEG, NEG, RIGHT);
    lut |= lutval(SEARCH_MIN, NEG, POS, RIGHT);
    lut |= lutval(SEARCH_MIN, POS, NEG, LEFT);
    lut |= lutval(SEARCH_MIN, POS, POS, LEFT); 

    printf("LUT: 0x%x\n", lut);
    return lut;
}

foundfeature_t bisectSearchJD(
        tle_t* tle,
        topo_pos_t observer,
        search_t searchtype,
        jd_ts precision,

        double startjd,
        double maxjd
        )
{

    double jd_1s = 1.0 / 86400; //1 second in decimal days

    double width = maxjd-startjd;
    //maxjd must be beyond startjd
    if(width < 0) {
        //TODO error! out of order.
        foundfeature_t x;
        x.feature = SEARCH_FAILURE;
        x.pos.err = -1;
        return x;
    }
    foundfeature_t ret;

    double halfwidth = width / 2;
    double midpoint = startjd + halfwidth;
    if(width <= precision) {
        sat_pos_t m = calcSat(tle, observer, midpoint);
        if ( (searchtype & SEARCH_MAX ) || (searchtype & SEARCH_MIN ) ){
            //trust that the bisection algorithm worked appropriately, assuming we have a min or max as appropriate within our search range and return
            ret.pos = m;
            return ret;
        }
        sat_pos_t l = calcSat(tle, observer, startjd);
        sat_pos_t r = calcSat(tle, observer, startjd+width);
        //rising, falling, and passmaxes require interaction with the horizon (elev = 0).
        //so we need to sanity check, as our parent/caller may need to continue searching
        if( (searchtype & SEARCH_RISING ) &&
                l.elev < m.elev && m.elev < r.elev && l.elev <= 0 && r.elev >= 0) {
            ret.pos = l;
            return ret;
        } 
        if ( (searchtype & SEARCH_FALLING ) &&
                l.elev > r.elev && l.elev >= 0 && r.elev <= 0) {
            ret.pos = l;
            return ret;
        } 
        if ( (searchtype & SEARCH_PASSMAX ) &&
                l.elev <= m.elev && m.elev >= r.elev && m.elev > 0) {
            ret.pos = m;
            return ret;
        } 
        printf("NOT FOUND\n");
        l.err = 1;
        ret.pos = l;
        return ret;
    } else {
        sat_pos_t m1 = calcSat(tle, observer, midpoint);
        sat_pos_t m2 = calcSat(tle, observer, midpoint+1*jd_1s);  //pop right a little so we can measure the change in elevation
        double mslope_elev = m2.elev - m1.elev; //positive mslope_elev represents increasing values left-to-right
        char elevation = m1.elev >= 0?1:0; //positive is truthy
        char slope = mslope_elev >= 0?1:0; //negative is falsy
        uint32_t key = lutkey(searchtype, slope, elevation);
        uint32_t lut = 0xe0f7b; //direction lookup table, generated by makelut
                              
        bool goRight = (lut >> key) & 1;
        /*printf("\tkey: 0x%x goRight: %d\n", key, goRight);*/

        if(goRight) {
            return bisectSearchJD(
                    tle, 
                    observer, 
                    searchtype,
                    precision, 
                    midpoint, 
                    maxjd
                    );
        } else {
            return bisectSearchJD(
                    tle, 
                    observer, 
                    searchtype,
                    precision, 
                    startjd, 
                    midpoint 
                    );
        }
    }
}



sat_pos_t search_simple_better(
        tle_t* tle,  //sat in question
        topo_pos_t obs, //observer location, in degrees latitude and longitude, and altitude in meters
        search_t featuretype, // SEARCH_RISING etc
        direction_t dir, //search direction

        jd_ts startjd //start time
        )
{
    //simple_search to the last min or max, then jump ahead and hillclimb to maxes until we find a pass
    jd_ts stepsize = 300.0/86400;
    jd_ts maxjd = startjd + 1;
    double mean_motion_rpd = tle->xno * MINS_PER_DAY / (2*PI);
    jd_ts one_orbit_jd = 1/mean_motion_rpd;
    foundfeature_t found = search_simple(
            tle,
            obs, 
            SEARCH_MIN|SEARCH_MAX,
            LEFT,
            stepsize,
            startjd, 
            maxjd
            );
    jd_ts nextjd = 0;
    if( found.feature & SEARCH_MIN ){
        //jump half an orbit 
        nextjd = found.pos.jd + one_orbit_jd/2;
    } else if( found.feature & SEARCH_MAX ){
        //jump 1 orbit 
        nextjd = found.pos.jd + one_orbit_jd;
    }
    stepsize = 30.0/86400;
    //subtract a bit so we don't miss it
    nextjd -= 2*stepsize;
    //nextjd should be just before a max, now
    if( featuretype == SEARCH_PASSMAX ){
        while( 1 ){
            found = search_simple(
                    tle,
                    obs, 
                    SEARCH_MIN|SEARCH_MAX|SEARCH_PASSMAX,
                    RIGHT,
                    stepsize,
                    nextjd, 
                    nextjd + one_orbit_jd
                    );
            nextjd = found.pos.jd + one_orbit_jd;
            stepsize = 30.0/86400;
            nextjd -= 3*stepsize;
            if( found.feature & featuretype ){
                return found.pos;
            }
            if( found.feature & SEARCH_MIN ){
                nextjd += one_orbit_jd/2 - 300.0/86400;
            }
            if( found.feature & SEARCH_MAX && found.pos.elev < -20){
                //don't keep this
                nextjd += 3*one_orbit_jd;
            }
        }
    } else {
        found = search_simple(
                tle,
                obs, 
                featuretype,
                RIGHT,
                stepsize,
                nextjd, 
                nextjd + one_orbit_jd
                );
        return found.pos;
    }
}
sat_pos_t nextpass_bisect_only(
        tle_t* tle,  //sat in question
        topo_pos_t obs, //observer location, in degrees latitude and longitude, and altitude in meters
        search_t featuretype, // SEARCH_RISING etc
        direction_t dir, //search direction

        jd_ts startjd //start time
        )
{
    //bisect to a max, then jump ahead and bisect to a max
}
sat_pass_t nextpass(
        tle_t* tle,  //sat in question
        topo_pos_t obs, //observer location, in degrees latitude and longitude, and altitude in meters
        jd_ts startjd //start time
        ){
    sat_pass_t pass = {0};
    /*pass.rise = get_orbit_local_feature( tle, obs, SEARCH_RISING,  RIGHT, startjd);*/
    /*pass.max  = get_orbit_local_feature( tle, obs, SEARCH_PASSMAX, RIGHT, pass.rise.jd);*/
    /*pass.set  = get_orbit_local_feature( tle, obs, SEARCH_FALLING, RIGHT, pass.max.jd);*/
    return pass;
}

errorcount get_next_N_passes(tle_t tle, topo_pos_t obs, jd_ts startjd, int N, sat_pass_t * passes, nextpassfn np){
    //find and store the next N passes in the sat_pass_t array passes.
    
    //This is meant to be used for writing automated and semi-automated
    //tests and test cases as a way of sanity checking more efficient
    //algorithms.

    jd_ts nextjd = startjd;
    sat_pass_t pass = {0};
    double mean_motion_rpd = tle.xno * MINS_PER_DAY / (2*PI);
    jd_ts one_orbit_jd = 1/mean_motion_rpd;
    jd_ts post_pass_offset = .05*one_orbit_jd; //how long after a pass before looking for the next pass
        //.05 is just to get it past the pass a bit, 
        //and 5% of an orbit period seems like it would be 
        //more consistant across satellites 
        //than some constant in seconds
    int errors = 0;
    for( int i = 0; i < N; i++ ){
        pass = np(tle,obs,nextjd);
        errors += pass.err; //yes, it will log the error but just continue trying anyway. 
                            //I have not defined the errors very well so it's more useful for now to collect them all, if any, 
                            //and I can look at them in gdb when they happen
        passes[i] = pass;
        nextjd = pass.set.jd + post_pass_offset; //actually, won't this fail all future passes if the pass isn't correctly calculated, 
                                                 //since the pass won't be fully filled out?
        assert( pass.set.jd != 0 ); //let's find out
    }
    return errors;
}




