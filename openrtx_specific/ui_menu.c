#include <maidenhead.h>
#include <satellite.h>
#include <math.h>

void _ui_drawMenuSatPredict(ui_state_t* ui_state){
    //actually, it's "asteroids"
    //this is a single game frame
    static unsigned long long last_t = 0;
    static game_obj_2d_t me = {0};
    static game_obj_2d_t asteroids[3] = {0};
    static int frame = 0;
    unsigned long long t = getTick();
    unsigned long long td = t - last_t;
    point_t center = {SCREEN_WIDTH/2, SCREEN_HEIGHT/2};
    /*printf("TD: %lu \n", td);*/
    if( last_t == 0 ){
        //init
        //everything should already be init'd to zero
        /*printf("init\n");*/
        me.x = center.x;
        me.y = center.y;
        me.rot = RAD(-45);
        game_obj_init( &asteroids[0] );
        game_obj_init( &asteroids[1] );
        game_obj_init( &asteroids[2] );

    } else {
        //regular call
        game_move(&me, td);
        if( ui_state->keys & KEY_4 ){
            me.rot -= RAD(2);
        }
        if( ui_state->keys & KEY_6 ){
            me.rot += RAD(2);
        }
        if( ui_state->keys & KEY_5 ){
            game_addvel(&me, 1, me.rot);
        }
        if( ui_state->keys & KEY_ENTER ){
            //fire
        }
        game_obj_screenwrap(&me);
        for( int i = 0; i < 3; i++ ){
            game_move(&asteroids[i], td);
            game_obj_screenwrap(&asteroids[i]);
        }
    }
    point_t mepos = {me.x, me.y};
    fflush(stdout);

    gfx_clearScreen();
    gfx_drawDeltaArrow(mepos, 8, DEG(me.rot)+90, yellow_fab413);
    if( ui_state->keys & KEY_5 ){
        //draw rocket exhaust
    }
    for( int i = 0; i < 3; i++ ){
        point_t apos = {asteroids[i].x, asteroids[i].y};
        gfx_drawCircle(apos, 4, color_white);
    }

    frame++;

    last_t = t;
}
void _ui_drawMenuSatSkyView(ui_state_t * ui_state, int whichsat, int passidx, int drawstars ){
    //whichsat == 0 means all sats (so you can pass in menu_selected directly)
    //whichpass == 0 means first pass
    //(this is partly because pass selection is not menu based, but based around custom inputs)
    topo_pos_t obs = getObserverPosition();
    double jd = curTime_to_julian_day(last_state.time);

    //handle zoom and pan events here
    //+- radius for zoom, +- plot_center for pan
    //0 means reset
    //
    //normalize this away from pixels
    //and store the zoom coords not as a modification of plot center
    //but as an azel to center! that way we can zoom in on whats under our viewport
    //without things sliding away on us
    static point_t viewport_center = {0, 0};
    int screen_scale = SCREEN_HEIGHT/2;
    static int radius = SCREEN_WIDTH/2/1.5;
    static int snap_item = -1;
    point_t plot_center = {SCREEN_WIDTH/2, SCREEN_HEIGHT/2+5};


    //TODO draw scale in degrees/arcseconds, etc
    //90/radius is degrees per pixel
    //debug prints
    //gfx_print(layout.bottom_pos, FONT_SIZE_5PT, TEXT_ALIGN_RIGHT, color_grey, "r=%dpx",radius );
    float deg_per_pix = 90.0 / radius;
    float pix_per_deg = radius/90.0;
    gfx_print(layout.bottom_pos, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_grey, "1px=%f deg", deg_per_pix);
    //gfx_print(layout.bottom_pos, FONT_SIZE_5PT, TEXT_ALIGN_RIGHT, color_grey, "trk%d", snap_item);

    if( ui_state->keys & KEY_0 ){
        radius = SCREEN_WIDTH/2/1.5;
        viewport_center.x = 0;
        viewport_center.y = 0;
    }
    //the more zoomed in we are, the more things i want to print for each satellite or star
    //az, el, speed, next pass, etc
    //as we zoom in a LOT and get more space
    //want to make the font larger too
    //
    if( ui_state->keys & KEY_UP ){
        radius *= 1.2;
        viewport_center.x *= 1.2;
        viewport_center.y *= 1.2;
    }
    if( ui_state->keys & KEY_DOWN ){
        radius /= 1.2;
        viewport_center.x /= 1.2;
        viewport_center.y /= 1.2;
    }
    if( ui_state->keys & KEY_2 ){
        viewport_center.y += fmax(screen_scale * 10/ radius, 1);
    }
    if( ui_state->keys & KEY_8 ){
        viewport_center.y -= fmax(screen_scale * 10/ radius, 1);
    }
    if( ui_state->keys & KEY_4 ){
        viewport_center.x += fmax(screen_scale * 10/ radius, 1);
    }
    if( ui_state->keys & KEY_6 ){
        viewport_center.x -= fmax(screen_scale * 10/ radius, 1);
    }
    if( ui_state->keys & KEY_1 ){
        snap_item++; //continued in whichsat == 0 section below!
    }
    plot_center = offset_point( plot_center, 1, viewport_center);

    gfx_drawPolarAzElPlot( plot_center, radius, color_grey );

    /*int retrograde = degrees(tle->xincl) > 90;*/
    //if( retrograde ) satellite passes east to west
    //else it passes west to east

    if( whichsat == 0 ){
        //draw current positions of all known sats
        for( int i = 0; i < num_satellites; i++){
            point_t temppos = {0,0};
            sat_mem_t selected = satellites[ i ];
            double jd_1s = 1.0 / 86400; //1 second in decimal days
            sat_pos_t sat = selected.current;
            if( sat.elev < -15 ){ //-15 is 15 degrees below the horizon, mind you! We want to show where things are for situational awareness
                continue;
            }
            gfx_drawPolar( plot_center, radius, sat.az, sat.elev, '+', yellow_fab413 );


            //draws a little trail of where it was two times in the past
            //hoping that shows motion well enough, with the fab413 -> white -> grey indicating age of data
            //ideally we'd save a number of track histories (with decreasing granularity as data ages)
            //and just draw those, rather than recalculating old stuff
            //...but that's a TODO later maybe, we'll see how this goes
            /*sat_pos_t sat_p2 = calcSat( selected.tle, jd-60*jd_1s, obs);*/
            /*sat_pos_t sat_p1 = calcSat( selected.tle, jd-30*jd_1s, obs);*/
            /*gfx_drawPolar( plot_center, radius, sat_p2.az, sat_p2.elev, 0, color_grey );*/
            /*gfx_drawPolar( plot_center, radius, sat_p1.az, sat_p1.elev, 0, color_white );*/


            point_t pos = azel_deg_to_xy( sat.az, sat.elev, radius);
            if( snap_item == i ){
                viewport_center.x = -1*pos.x;
                viewport_center.y = -1*pos.y;
                //mult -1 because the viewport gets shifted away in the opposite direction
                //(if you want the (1,1) corner of a piece of paper to
                //be centered on your desk when the center of the desk and
                //center of the paper are in the same place, you have to
                //move the paper in the (-1,-1) direction
            }
            point_t text_offset = {6,6};

            temppos = offset_point( plot_center, 2, pos, text_offset);
            gfx_print(temppos, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_white,
                    "%s", selected.name);
        }
        if( snap_item >= num_satellites ){ //if snap_item is 1 and num_satellites is 1, zero indexing means we're past it and should reset already
            snap_item = -1; //reset to not tracking
        }
        if( snap_item >= 0 ){
            sat_mem_t trk = satellites[ snap_item ];
            gfx_print(layout.bottom_pos, FONT_SIZE_5PT, TEXT_ALIGN_RIGHT, color_white, "TRK: %s", trk.name);
        }
    } else { //if we're drawing a specific satellite...
        (void) passidx;
        sat_mem_t selected = satellites[ whichsat-1 ]; //-1 because whichsat is one indexed (and C is zero indexed)
        
        //draw a pass and current position for only a specific satellite
        //Beneath the pass, allow showing a doppler plot for sanity checking doppler correction later
        //add a setting for enabling doppler correction
        #define PASS_PTS 20
        static sat_azel_t pass_azel[PASS_PTS] = {0};
        static int pass_pos_idx = 0;
        static int lastsat = -1;
        static int lastpass = -1;
        if( lastsat == -1 && lastpass == -1 ){
            lastsat = whichsat;
            lastpass = passidx;
        }
        if( lastsat != whichsat || lastpass != passidx ){
            lastsat = whichsat;
            lastpass = passidx;
            pass_pos_idx = 0;
        }
        double jd_width = selected.nextpass.set.jd - selected.nextpass.rise.jd;
        double jd_step = jd_width / PASS_PTS;
        #define PASS_PTS_PER_CALL 4
        for( int j = 0; j < PASS_PTS_PER_CALL && pass_pos_idx < PASS_PTS; j++ ){
            double jd_at = selected.nextpass.rise.jd + pass_pos_idx * jd_step;
            sat_pos_t prediction = calcSat( &(selected.tle), jd_at , obs);
            pass_azel[pass_pos_idx].jd = prediction.jd;
            pass_azel[pass_pos_idx].az = prediction.az;
            pass_azel[pass_pos_idx].elev = prediction.elev;
            pass_azel[pass_pos_idx].dist = prediction.dist;
            pass_pos_idx++;
        }
        //int num_points_pass = sizeof(pass_azel) / sizeof(sat_azel_t);
        int num_points_pass = pass_pos_idx;


        for( int i = 0; i < num_points_pass; i+=1 ){
            gfx_drawPolar( plot_center, radius, pass_azel[i].az, pass_azel[i].elev,
                    0, //0 here means set a pixel, don't draw a character
                    color_white );
        }
        point_t rise_rel = azel_deg_to_xy( pass_azel[0].az, pass_azel[0].elev, radius);
        point_t left_text_offset = {-26, 3};
        point_t set_rel = azel_deg_to_xy( pass_azel[num_points_pass-1].az, pass_azel[num_points_pass-1].elev, radius);
        point_t right_text_offset = {9, 3};
        point_t rise = offset_point( plot_center, 2, rise_rel, left_text_offset );
        point_t set = offset_point( plot_center, 2, set_rel, right_text_offset );

        //at start and end points, print rise and set time
        //need to convert jd to hours/mins if i want time here
        //gfx_print(rise, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_white, "%02d:%02d", 4, 53);
        //gfx_print(set, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_white, "%02d:%02d", 4, 59);


        point_t line_offset_5pt = {0,8};
        point_t temppos = {0,0};

        temppos = offset_point( rise, 1, line_offset_5pt );
        gfx_print(temppos, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_white,
                "%03.0f AZ", pass_azel[0].az);

        gfx_drawPolar( plot_center, radius, selected.nextpass.max.az, selected.nextpass.max.elev, 
                '^', color_white );
        point_t max_rel = azel_deg_to_xy( selected.nextpass.max.az, selected.nextpass.max.elev, radius);
        point_t max = offset_point( plot_center, 2, max_rel, left_text_offset );
        temppos = offset_point( max, 1, line_offset_5pt );
        gfx_print(temppos, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_white,
                "%02.0f EL", selected.nextpass.max.elev);

        temppos = offset_point( set, 1, line_offset_5pt );
        gfx_print(temppos, FONT_SIZE_5PT, TEXT_ALIGN_LEFT, color_white,
                "%03.0f AZ", pass_azel[num_points_pass-1].az);


        double az = selected.current.az;
        double elev = selected.current.elev;
        gfx_drawPolar( plot_center, radius, az, elev, '+', yellow_fab413 );
    }


    if( drawstars ){

        for( int i = 0; i < num_stars; i+=1 ){
            double az, alt;
            double ra, dec;
            ra  = stars[i].ra*15; //in decimal hours, so *15->deg
            dec = stars[i].dec;
            ra_dec_to_az_alt(jd, obs.lat, obs.lon, ra, dec, &az, &alt);
            az = az;
            alt = alt;
            /*printf("%.1f %.1f %.0f %.0f\n", ra, dec, az, alt);*/
            if( alt < -2.5 ){
                //draw stars that are juuust below the horizon
                continue;
            }
            int brt_scale = 0xff;
            if( stars[i].mag > 5 && radius < screen_scale*3 ){
                //too dim, not zoomed in enough to see it
                continue;
            }
            if( radius > screen_scale * 3){
                brt_scale <<= 2;
            }
            if( radius > screen_scale * 8){
                brt_scale <<= 2;
            }
            int brt = pow(2.51,-1*stars[i].mag)*brt_scale;
            if( brt > 0xff ) brt = 0xff;
            if( brt < 0x10 ) brt = 0x10;

            color_t clr = {0xff, 0xff, 0x00, brt};
            color_t clr2 = {0xff, 0xff, 0x00, 0x30};
            if( i < pow(log2(fmax(radius,58)-55),2) || radius > screen_scale * 12){ //just an eyeball magic number
                //the more zoomed in we are, the more names of stars we print!
                point_t relpos = azel_deg_to_xy( az, alt, radius);
                point_t offset = {3,3};
                fontSize_t ft= FONT_SIZE_5PT;
                point_t pos = offset_point( plot_center, 2, relpos, offset );
                gfx_print(pos, ft, TEXT_ALIGN_LEFT, clr2, stars[i].name);
            }
            gfx_drawPolar( plot_center, radius, az, alt, 0, clr );
        }
    }


}

void _ui_drawMenuSatTrack(ui_state_t * ui_state)
{

    char gridsquare[7] = {0}; //we want to use this as a c-style string so that extra byte stays zero

    gfx_clearScreen();
    _ui_drawMainTop();
    topo_pos_t obs = getObserverPosition();
    double jd = curTime_to_julian_day(last_state.time);

    if( ui_state->menu_selected == 0 ){
        //first menu entry (0) is to show all satellites in sky view
        _ui_drawMenuSatSkyView(ui_state, 0, 0,
                1 //draw stars?
                );
        return;
    }
    _ui_drawMainBackground();
    //_ui_drawBottom();
    int idx = ui_state->menu_selected - 1; //because 0 is "all sats"
    sat_mem_t selected = satellites[ idx ];
    sat_pos_t sat = selected.current;
    //TODO recalculate when we re-enter from a different satellite entry
    double tdiff = 0;
    double nextpass = selected.nextpass.rise.jd;

    tdiff = (nextpass - jd) * 24 * 60 * 60; //seconds
                                            //
    gfx_print(layout.bottom_pos, FONT_SIZE_8PT, TEXT_ALIGN_LEFT, color_white,
            "AOS %ds", (int)tdiff);

    // left side
    // relative coordinates to satellite
    gfx_print(layout.line2_pos, FONT_SIZE_8PT, TEXT_ALIGN_LEFT, color_white,
            "AZ %.1f", sat.az);

    gfx_print(layout.line2_pos, FONT_SIZE_8PT, TEXT_ALIGN_RIGHT, color_white,
            "EL %.1f", sat.elev);

    gfx_print(layout.line1_pos, FONT_SIZE_8PT, TEXT_ALIGN_CENTER, color_white, selected.name);

    //right side
    //doppler correction readout
    /*snprintf(sbuf, 25, "%.1fk DOP", ((float)doppler_offset)/1000);*/
    /*gfx_print(layout.line1_pos, sbuf, FONT_SIZE_8PT, TEXT_ALIGN_RIGHT, color_white);*/
    //draw gridsquare text

    int gslen = strnlen(last_state.settings.gridsquare_override, 7);
    if( gslen > 0 ){
        //adapt maidenhead printed length to whatever was specified
        //we could just use the maidenhead as is, but i like this as a sanity check
        //to make sure the location used is what we expect. 
        //This way the displayed gridsquare is definitely the location used because obs is a single source of truth
        lat_lon_to_maidenhead(obs.lat, obs.lon, gridsquare, gslen/2); 
    } else {
        lat_lon_to_maidenhead(obs.lat, obs.lon, gridsquare, 3); //precision=3 here means 6 characters like FN41uq
    }
    gfx_print(layout.line3_pos, FONT_SIZE_8PT, TEXT_ALIGN_LEFT, color_white, gridsquare);

    //center bottom - show
    //satellite and AOS/LOS countdown

    //draw Az/El
    int radius = SCREEN_WIDTH/2/4;
    point_t plot_center = {SCREEN_WIDTH/2, SCREEN_HEIGHT/2+5};
    gfx_drawPolarAzElPlot( plot_center, radius, color_grey );
    /*for( int i = 0; i < num_points_pass; i+=2 ){*/
        /*gfx_drawPolar( plot_center, radius, pass_azel[i].az,*/
            /*pass_azel[i].elev, 0, color_white );*/
    /*}*/
    gfx_drawPolar( plot_center, radius, sat.az, sat.elev, '+', yellow_fab413 );

    /*
    char * pass_state;
    double mark;
    double pass_start_jd = 0; //pass_azel[0].jd;
    double pass_end_jd = 0; //pass_azel[num_points_pass-1].jd;
    //mark is the time we care about most - e.g. before a pass, it's the time the pass starts
    //or during a pass, it's the LOS time (when the sat will go below the horizon)
    if( jd < pass_start_jd ){
        //before this pass comes over the horizon
        pass_state = "AOS";
        mark = pass_start_jd;
    } else if ( jd > pass_start_jd && jd < pass_end_jd ){
        //during the pass
        pass_state = "LOS";
        mark = pass_end_jd;
    } else {
        //now it's gone over the horizon, so same as the elif above (just will be a
        //negative number to show it was in the past)
        //left here for clarity to show the actual LOS condition
        pass_state = "LOS";
        mark = pass_end_jd;
    }
    float diff = (mark - jd)*86400; //diff is seconds until (+) or since (-) the mark timestamp
    snprintf(sbuf, 25, "%s %s %.0fs", sat_name, pass_state, diff);
    gfx_print(layout.line3_pos, sbuf, FONT_SIZE_8PT, TEXT_ALIGN_CENTER, color_white);
    */
}

