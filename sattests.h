#ifndef MM_SATTESTS_H
#define MM_SATTESTS_H
extern char algo[];
soe passes_equalish( sat_pass_t p1, sat_pass_t p2, jd_ts precision_time, double precision_degrees);
errorcount compare_N_passes_between(tle_t tle, topo_pos_t obs, jd_ts jd, int N,
        nextpassfn np1, char * np1name,
        nextpassfn np2, char * np2name
        );

void test_1();
void test_2();
void test_3();
#endif
