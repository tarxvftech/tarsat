
sat_pos_t calcSatNow( tle_t * tle ){
    double jd;
    topo_pos_t obs = getObserverPosition();
    jd = curTime_to_julian_day(last_state.time);
    return calcSat(tle, jd, obs);
}
void init_sat_global()
{
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
int num_stars = sizeof(stars) / sizeof(star_t);
