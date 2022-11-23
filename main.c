#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "satellite.h"
#include "maidenhead.h"

void print_pos(double jd, topo_pos_t obs, sat_pos_t p){
   printf("NORAD:\t%d\tErrors:\t%d\n", p.satid, p.err);
   printf("Az:   \t%.2f\tEl:\t%.2f\tRange:\t%.1f\n", p.az, p.elev, p.dist);
   printf("\t\t\t\tRaDec \t%.4f %.4f\n", p.ra, p.dec);
}
void print_n2yo(tle_t tle){
   printf("https://n2yo.com/?s=%d\n", tle.norad_number);
}

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

int main(int argc, char **argv){
   size_t num_satellites = 0;
   topo_pos_t obs = {0};
   latlon me = {0}; //need this because topo_pos_t is floats but maidenhead uses doubles
   me = maidenhead_to_latlon("FN41");
   obs.lat = me.lat;
   obs.lon = me.lon;

   sat_mem_t satellites[10] = {0};
#include "generated/sat_tles.c"
   printf("LL: %.8f\t%.8f\n", obs.lat, obs.lon);
   double start_jd = getjd();
   double jd_inc = 300.0/86400;
   for( int i = 0; i < num_satellites; i++ ){
      char fn[64] = {0};
      snprintf(fn, 64, "%s.csv", satellites[i].name);
      FILE * fd = fopen(fn,"w");
      for( double jd = start_jd; jd < start_jd+10; jd+=jd_inc) {
         /*printf("JD: %f\n", jd);*/
         /*for( int i = 0; i < 1; i++ ){*/
         /*int x = printf("\n%s\t", satellites[i].name );*/
         /*if( x < 8 ){*/
            /*printf("\t");*/
         /*}*/
         sat_pos_t current = calcSat( &satellites[i].tle, jd, obs);
         /*print_pos(jd, obs, current );*/
         /*printf("\t\t\t\t"); print_n2yo(satellites[i].tle);*/
         fprintf(fd,"%f,%f,%f,%f,%f,%f\n", jd, current.az, current.elev, current.dist, current, current.dec);
      }
      fclose(fd);
   }

   return 0;
}
