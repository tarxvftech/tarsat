
sat_pos_t bad_bisectSearchJD(
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
sat_pass_t bad_sat_nextpass(
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


