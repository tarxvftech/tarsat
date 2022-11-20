
sat_pos_t calcSatNow( tle_t * tle ){
    double jd;
    topo_pos_t obs = getObserverPosition();
    jd = curTime_to_julian_day(last_state.time);
    return calcSat(tle, jd, obs);
}
