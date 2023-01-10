#ifndef MINUNIT_H
#define MINUNIT_H

//#include <assert.h>
//assert(truthy)

//or test with python
//https://news.ycombinator.com/item?id=14293984


//minunit
//https://jera.com/techinfo/jtns/jtn002
#define mu_assert(message, test) do { if (!(test)) return message; } while (0)
#define mu_run_test(test) do { char *message = test(); tests_run++; \
                                if (message) return message; } while (0)
extern int tests_run;

#endif
