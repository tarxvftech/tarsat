#ifndef MM_LINUX_H
#define MM_LINUX_H
jd_ts unix2jd(time_t unix_epoch);
jd_ts getjd();
sof loadtle(char * filename, int catalog_number, tle_t * tle);

#endif
