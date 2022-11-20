/**
 * \internal Satellite calculations thread
 * to offload that processing from the UI thread
 * Runs through all the satellites and checks if it needs to update the satellite state.
 * For each sat, it calculates the current position, next N passes if needed, etc
 *
 * TODO Can I have this disabled by default?
 *
 */
void *satcalc_task(void *arg)
{
    (void) arg;
    while(1)
    {
        if( satellites_initialized == 0 ){ //ui task ran and inited satellite data
            //double jd = curTime_to_julian_day(last_state.time);  //crash due to bad data from "RTC"
            sleepFor(2u, 0u); //delay 2s to test something TEMPORARY
            //double jd = curTime_to_julian_day(last_state.time);  //nocrash
        }
        if( satellites_initialized >= 0 ){ 
            //printf("start calculating satellites...\r\n");
            topo_pos_t obs = getObserverPosition(); //implicitly depends on last_state settings/gps
            if( obs.lat != 0 ){
                double jd = curTime_to_julian_day(last_state.time); 
                //only search for passes that are coming up soon and then don't update them until they pass completely
                //make a special UI page for calculating passes with feedback and UI around which sats when, and to what criteria (minimum elevation, etc)
                for( int i = 0; i < num_satellites; i++){
                    sat_mem_t * sat = &satellites[i]; 
                    //printf("Sat %d %s\r\n", i, sat->name);
                    //C structures assign by value...so if i want to change the global satellites i need to be careful to actually change the global structure!
                    
                    sat->current = calcSat( &sat->tle, jd, obs ); 

                    double tdiff = 0;
                    double nextpass = sat->nextpass.rise.jd;
                    tdiff = (nextpass - jd) * 24; //hours * 60; //minutes * 60; //seconds
                    //printf("\n\nSat %d %s jd %f nextpass %f tdiff %fhrs\n", i, sat->name, jd, nextpass, tdiff);
                    if( sat->nextpass.err != 0 || tdiff <= 0 ){
                        //printf("recalc sat %d pass %d\n",i,0);
                        sat->nextpass = sat_nextpass( &sat->tle, jd, obs);
                        if( sat->nextpass.err == 0 ){
                            //printf("pass found!!\n");
                        } else {
                            //printf("no pass found in search time\n");
                            sat->nextpass.err = -1;
                        }
                    }
                }
                if( satellites_initialized == 0 ){ //ui task ran and inited satellite data
                    satellites_initialized = 1; //indicate ui task has valid data to play with
                }
            }
            //printf("done calculating satellites...\r\n");
        }

        sleepFor(2u, 0u); //update every 2s
    }
}

