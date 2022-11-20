
double curTime_to_julian_day(curTime_t t)
{
    //t appears to correctly be UTC which we expect, good.
    //this appears to be correct to within at least a second

    //expects t to be after year 2000
    //many thanks to Peter Baum, and his "Date Algorithms" reference.
    uint8_t s = t.second;
    uint8_t m = t.minute;
    uint8_t h = t.hour;
    uint8_t D = t.date; //.day is the _weekday_, date is day-of-month
    uint8_t M = t.month;
    short Y = CENTURY + t.year; //CENTURY is from rtc interface
    /*printf("%04d/%02d/%02d  %02d:%02d:%02d\n", Y,M,D,h,m,s);*/

    int Z = Y + (M - 14) / 12; //relies on int truncation
    const short Fvec[] = {306, 337, 0, 31, 61, 92, 122, 153, 184, 214, 245, 275};
    short F = Fvec[M - 1];

    //! note difference between floor and int truncation
    //floor(-1.5) => -2
    //int(-1.5) => -1
    double jd = D + F + 365 * Z + floor(Z / 4) - floor(Z / 100) + floor(Z / 400) + 1721118.5;
    //that +.5 is because julian .0 is actually _noon_, not midnight
    //so JD #.5 is _halfway through that julian day_, at _midnight_
    //so we add hours, minutes, and seconds into a fractional day and add that to our new "midnight" epoch and...
    jd += ((float)h) / 24 + ((float)m) / 1440 + ((float)s) / 86400;
    //voila!

    /*printf("jd: %.6f\n",jd);*/
    return jd;
}
topo_pos_t getObserverPosition(){
    topo_pos_t obs = {0};
    if( strnlen(last_state.settings.gridsquare_override,7) > 0 ){
       maidenhead_to_lat_lon( last_state.settings.gridsquare_override, 
             &obs.lat, &obs.lon);
       obs.alt = 0; //msl geoid meters //redundant assignment, just being
       //explicit that it's there in case we want to add manual altitude
       //input later (for SOTA, if we're doing more than just satellite work later with this data?)
   #ifdef GPS
    } else {
        //TODO: need a way to show gps enabled/disable, gps fix/nofix
          /*|| !last_state.settings.gps_enabled || last_state.gps_data.fix_quality == 0 */
        //fix_type is 1 sometimes when it shouldn't be, have to use fix_quality
        //gfx_print(layout.line3_pos, "no gps fix", FONT_SIZE_12PT, TEXT_ALIGN_CENTER, color_white);
        obs.lat = last_state.gps_data.latitude;
        obs.lon = last_state.gps_data.longitude;
        obs.alt = last_state.gps_data.altitude; //msl geoid meters
   #endif
    }
    return obs;
}
void game_move(game_obj_2d_t * o, unsigned long long td){
    o->x += (o->spdx)*td/1000;
    o->y += (o->spdy)*td/1000;
}
void game_addvel( game_obj_2d_t * o, float vel, float rot ){
    float vy = vel * sin(rot);
    float vx = vel * cos(rot);
    o->spdx += vx;
    o->spdy += vy;
}
void game_obj_init( game_obj_2d_t * o ){
    o->x = rand() % SCREEN_WIDTH;
    o->y = rand() % SCREEN_HEIGHT;
    o->rot = ((float)(rand() % ((int)(2*PI*100))))/100;
    float vel = (rand() % 10);
    game_addvel( o, vel, o->rot);
}
void game_obj_screenwrap( game_obj_2d_t * o ){
    if( o->x < 0 ){
        o->x = SCREEN_WIDTH + o->x;
    }
    if( o->x >= SCREEN_WIDTH ){
        o->x = o->x - SCREEN_WIDTH;
    }
    if( o->y < 0 ){
        o->y = SCREEN_HEIGHT + o->y;
    }
    if( o->y >= SCREEN_HEIGHT ){
        o->y = o->y - SCREEN_HEIGHT;
    }
}
