#ifndef COMMON_H
#define COMMON_H

#define DEF_PI 3.141592654
#define DEF_2PI 6.28318531
#define DEF_RE 6.3712e+06
#define DEF_TO_RAD (DEF_PI/180.)
#define DEF_Q 1.6021e-19
#define DEF_C 2.99725e+08
#define DEF_HM0_E 9.1093826e-31
#define DEF_HM0_P 1.6725e-27

typedef union {
    struct {
        double x;
        double y;
        double z;
    };
    struct {
        double r;
        double theta;
        double phi;
    };
} Vector;

typedef struct {
    double r0;
    double the0;
    double fi0;
    double the1;
    double fi1;
    int iy;
    int mes;
    int ide;
    int id;
    int ih;
    int min;
    int is;
    union {
        struct {
            double zn;
            int nk1;
            int iopt;
            int ist;
            double rig;
            double del;
            double rk;
        };
        struct {
            double lat_step;
            double lon_step;
        };
    };
} InfileHeader;

#endif // COMMON_H
