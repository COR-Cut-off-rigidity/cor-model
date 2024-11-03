#include <math.h>
#include "full_rc.h"

void full_rc(double X, double Y, double Z, double PHI,
             double SC_SY, double SC_PR, double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC,
             double *BYPRC, double *BZPRC) {

    static const double C_SY[86] = {-957.2534900,-817.5450246,583.2991249,758.8568270,
                             13.17029064,68.94173502,-15.29764089,-53.43151590,27.34311724,
                             149.5252826,-11.00696044,-179.7031814,953.0914774,817.2340042,
                             -581.0791366,-757.5387665,-13.10602697,-68.58155678,15.22447386,
                             53.15535633,-27.07982637,-149.1413391,10.91433279,179.3251739,
                             -6.028703251,1.303196101,-1.345909343,-1.138296330,-0.06642634348,
                             -0.3795246458,.07487833559,.2891156371,-.5506314391,-.4443105812,
                             0.2273682152,0.01086886655,-9.130025352,1.118684840,1.110838825,
                             .1219761512,-.06263009645,-.1896093743,.03434321042,.01523060688,
                             -.4913171541,-.2264814165,-.04791374574,.1981955976,-68.32678140,
                             -48.72036263,14.03247808,16.56233733,2.369921099,6.200577111,
                             -1.415841250,-0.8184867835,-3.401307527,-8.490692287,3.217860767,
                             -9.037752107,66.09298105,48.23198578,-13.67277141,-16.27028909,
                             -2.309299411,-6.016572391,1.381468849,0.7935312553,3.436934845,
                             8.260038635,-3.136213782,8.833214943,8.041075485,8.024818618,
                             35.54861873,12.55415215,1.738167799,3.721685353,23.06768025,
                             6.871230562,6.806229878,21.35990364,1.687412298,3.500885177,
                             0.3498952546,0.6595919814};

    static const double C_PR[86] = {-64820.58481,-63965.62048,66267.93413,135049.7504,
                             -36.56316878,124.6614669,56.75637955,-87.56841077,5848.631425,
                             4981.097722,-6233.712207,-10986.40188,68716.52057,65682.69473,
                             -69673.32198,-138829.3568,43.45817708,-117.9565488,-62.14836263,
                             79.83651604,-6211.451069,-5151.633113,6544.481271,11353.03491,
                             23.72352603,-256.4846331,25.77629189,145.2377187,-4.472639098,
                             -3.554312754,2.936973114,2.682302576,2.728979958,26.43396781,
                             -9.312348296,-29.65427726,-247.5855336,-206.9111326,74.25277664,
                             106.4069993,15.45391072,16.35943569,-5.965177750,-6.079451700,
                             115.6748385,-35.27377307,-32.28763497,-32.53122151,93.74409310,
                             84.25677504,-29.23010465,-43.79485175,-6.434679514,-6.620247951,
                             2.443524317,2.266538956,-43.82903825,6.904117876,12.24289401,
                             17.62014361,152.3078796,124.5505289,-44.58690290,-63.02382410,
                             -8.999368955,-9.693774119,3.510930306,3.770949738,-77.96705716,
                             22.07730961,20.46491655,18.67728847,9.451290614,9.313661792,
                             644.7620970,418.2515954,7.183754387,35.62128817,19.43180682,
                             39.57218411,15.69384715,7.123215241,2.300635346,21.90881131,
                             -.01775839370,.3996346710};

    double HXSRC, HYSRC, HZSRC, HXPRC, HYPRC, HZPRC;
    src_prc(SC_SY, SC_PR, PHI, X, Y, Z, &HXSRC, &HYSRC, &HZSRC, &HXPRC, &HYPRC, &HZPRC);

    double X_SC = SC_SY - 1.;
    double FSX, FSY, FSZ;
    rc_shield(C_SY, X_SC, X, Y, Z, &FSX, &FSY, &FSZ);

    X_SC = SC_PR - 1.;
    double FPX, FPY, FPZ;
    rc_shield(C_PR, X_SC, X, Y, Z, &FPX, &FPY, &FPZ);

    *BXSRC = HXSRC + FSX;
    *BYSRC = HYSRC + FSY;
    *BZSRC = HZSRC + FSZ;

    *BXPRC = HXPRC + FPX;
    *BYPRC = HYPRC + FPY;
    *BZPRC = HZPRC + FPZ;
}

void src_prc(double SC_SY, double SC_PR, double PHI,
             double X, double Y, double Z, double *BXSRC, double *BYSRC, double *BZSRC,
             double *BXPRC, double *BYPRC, double *BZPRC) {

    const double XT = X;
    const double ZT = Z;

    const double XTS = XT / SC_SY;
    const double YTS = Y / SC_SY;
    const double ZTS = ZT / SC_SY;

    const double XTA = XT / SC_PR;
    const double YTA = Y / SC_PR;
    const double ZTA = ZT / SC_PR;

    double BXS = 0.;
    double BYS = 0.;
    double BZS = 0.;
    rc_symm(XTS, YTS, ZTS, &BXS, &BYS, &BZS);

    double BXA_S = 0.;
    double BYA_S = 0.;
    double BZA_S = 0.;
    prc_symm(XTA, YTA, ZTA, &BXA_S, &BYA_S, &BZA_S);

    const double CP = cos(PHI);
    const double SP = sin(PHI);
    const double XR = XTA*CP - YTA*SP;
    const double YR = XTA*SP + YTA*CP;

    double BXA_QR = 0.;
    double BYA_QR = 0.;
    double BZA_Q = 0.;
    prc_quad(XR, YR, ZTA, &BXA_QR, &BYA_QR, &BZA_Q);

    const double BXA_Q = BXA_QR*CP + BYA_QR*SP;
    const double BYA_Q = -BXA_QR*SP + BYA_QR*CP;
    const double BXP = BXA_S + BXA_Q;
    const double BYP = BYA_S + BYA_Q;
    const double BZP = BZA_S + BZA_Q;

    *BXSRC = BXS;
    *BYSRC = BYS;
    *BZSRC = BZS;

    *BXPRC = BXP;
    *BYPRC = BYP;
    *BZPRC = BZP;
}

void rc_symm(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double DS = 1.E-2;
    const double DC = 0.99994999875;
    const double D = 1.E-4;
    const double DRD = 5E3;

    const double RHO2 = X*X + Y*Y;
    const double R2 = RHO2 + Z*Z;
    const double R = sqrt(R2);
    const double RP = R + D;
    const double RM = R - D;
    const double SINT = sqrt(RHO2)/R;
    const double COST = Z/R;

    if(SINT < DS) {
        const double A = ap(R, DS, DC)/DS;
        const double DARDR = (RP*ap(RP, DS, DC) - RM*ap(RM, DS, DC))*DRD;
        const double FXY = Z*(2.*A - DARDR)/(R*R2);
        *BX = FXY*X;
        *BY = FXY*Y;
        *BZ = (2.*A* (COST*COST)+DARDR*(SINT*SINT))/R;
    } else {
        const double THETA = atan2(SINT, COST);
        const double TP = THETA + D;
        const double TM = THETA - D;
        const double SINTP = sin(TP);
        const double SINTM = sin(TM);
        const double COSTP = cos(TP);
        const double COSTM = cos(TM);
        const double BR = (SINTP*ap(R, SINTP, COSTP) - SINTM*ap(R, SINTM, COSTM))/(R*SINT)*DRD;
        const double BT = (RM*ap(RM, SINT, COST) - RP*ap(RP, SINT, COST))/R*DRD;
        const double FXY = (BR + BT*COST/SINT)/R;
        *BX = FXY*X;
        *BY = FXY*Y;
        *BZ = BR*COST - BT*SINT;
    }
}

double ap(double R, double SINT, double COST) {
    const double A1 = -456.5289941;
    const double A2 = 375.9055332;
    const double RRC1 = 4.274684950;
    const double DD1 = 2.439528329;
    const double RRC2 = 3.367557287;
    const double DD2 = 3.146382545;
    const double P1 = -0.2291904607;
    const double R1 = 3.746064740;
    const double DR1 = 1.508802177;
    const double DLA1 = 0.5873525737;
    const double P2 = 0.1556236119;
    const double R2 = 4.993638842;
    const double DR2 = 3.324180497;
    const double DLA2 = 0.4368407663;
    const double P3 = 0.1855957207;
    const double R3 = 2.969226745;
    const double DR3 = 2.243367377;

    int PROX = 0;
    double SINT1 = SINT;
    double COST1 = COST;
    if(SINT1 < 1.E-2) {
        SINT1 = 1.E-2;
        COST1 = .99994999875;
        PROX = 1;
    }

    const double R_R1_DR1 = (R - R1) / DR1;
    const double R_R2_DR2 = (R - R2) / DR2;
    const double R_R3_DR3 = (R - R3) / DR3;
    const double COST1_DLA1 = COST1 / DLA1;
    const double COST1_DLA2 = COST1 / DLA2;
    const double ARG1 = -(R_R1_DR1*R_R1_DR1) - (COST1_DLA1*COST1_DLA1);
    const double ARG2 = -(R_R2_DR2*R_R2_DR2) - (COST1_DLA2*COST1_DLA2);
    const double ARG3 = -(R_R3_DR3*R_R3_DR3);

    double DEXP1 = 0.;
    if(ARG1 >= -500.) {
        DEXP1 = exp(ARG1);
    }

    double DEXP2 = 0.;
    if(ARG2 >= -500.) {
        DEXP2 = exp(ARG2);
    }

    double DEXP3 = 0.;
    if(ARG3 >= -500.) {
        DEXP3 = exp(ARG3);
    }

    const double ALPHA = (SINT1*SINT1) / R;
    const double ALPHA_S = ALPHA*(1. + P1*DEXP1 + P2*DEXP2 + P3*DEXP3);

    const double GAMMA = COST1 / (R*R);
    const double GAMMA_S = GAMMA;
    const double GAMMAS2 = GAMMA_S*GAMMA_S;

    const double ALSQH = (ALPHA_S*ALPHA_S)/2.;
    const double F = 64./27.*GAMMAS2 + (ALSQH*ALSQH);
    const double Q = pow((sqrt(F) + ALSQH),(1./3.));
    double C = Q - 4.*pow(GAMMAS2,(1./3.))/(3.*Q);
    if(C < 0.){
        C = 0.;
    }
    const double G = sqrt(C*C + 4.*pow(GAMMAS2,(1./3.)));
    const double RS = 4./((sqrt(2.*G - C)+sqrt(C))*(G + C));
    const double COSTS = GAMMA_S * (RS*RS);
    const double SINTS = sqrt(1. - COSTS*COSTS);
    const double RHOS = RS*SINTS;
    const double ZS = RS*COSTS;
    const double ZS2 = ZS*ZS;
    const double RRC1_RHOS = RRC1 + RHOS;

    double P = RRC1_RHOS*RRC1_RHOS + ZS2 + DD1*DD1;
    double XK2 = 4.*RRC1*RHOS/P;
    double XK = sqrt(XK2);
    double XKRHO12 = XK*sqrt(RHOS);

    double XK2S = 1. - XK2;
    double DL = log(1./XK2S);
    double ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383 + XK2S*(0.03742563713 + XK2S*0.01451196212))) +
                 DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    double ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
                 DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI1 = ((1. - XK2*0.5)*ELK - ELE)/XKRHO12;
    const double RRC2_RHOS = RRC2 + RHOS;
    P = RRC2_RHOS*RRC2_RHOS + ZS2 + DD2*DD2;
    XK2 = 4.*RRC2*RHOS/P;
    XK = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);
    XK2S = 1. - XK2;
    DL = log(1./XK2S);
    ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383 + XK2S*(0.03742563713 + XK2S*0.01451196212))) +
          DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
          DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI2 = ((1. - XK2*0.5)*ELK-ELE)/XKRHO12;
    double AP = A1*APHI1 + A2*APHI2;

    if(PROX) {
        AP = AP*SINT/SINT1;
    }

    return AP;
}

void prc_symm(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double DS = 1.E-2;
    const double DC = 0.99994999875;
    const double D = 1.E-4;
    const double DRD = 5.E3;

    const double RHO2 = X*X + Y*Y;
    const double R2 = RHO2 + Z*Z;
    const double R = sqrt(R2);
    const double RP = R + D;
    const double RM = R - D;
    const double SINT = sqrt(RHO2)/R;
    const double COST = Z/R;

    double FXY;
    if(SINT < DS) {
        const double A = apprc(R, DS, DC)/DS;
        const double DARDR = (RP*apprc(RP, DS, DC) - RM*apprc(RM, DS, DC))*DRD;
        FXY = Z*(2.*A - DARDR)/(R*R2);
        *BZ = (2.*A* (COST*COST) + DARDR*(SINT*SINT))/R;
    } else {
        const double THETA = atan2(SINT,COST);
        const double TP = THETA + D;
        const double TM = THETA - D;
        const double SINTP = sin(TP);
        const double SINTM = sin(TM);
        const double COSTP = cos(TP);
        const double COSTM = cos(TM);
        const double BR = (SINTP*apprc(R, SINTP, COSTP) - SINTM*apprc(R, SINTM, COSTM))/(R*SINT)*DRD;
        const double BT = (RM*apprc(RM, SINT, COST) - RP*apprc(RP, SINT, COST))/R*DRD;
        FXY = (BR + BT*COST/SINT)/R;
        *BZ = BR*COST - BT*SINT;
    }

    *BX = FXY*X;
    *BY = FXY*Y;
}

double apprc(double R, double SINT, double COST) {
    const double A1 = -80.11202281;
    const double A2 = 12.58246758;
    const double RRC1 = 6.560486035;
    const double DD1 = 1.930711037;
    const double RRC2 = 3.827208119;
    const double DD2 = .7789990504;
    const double P1 = .3058309043;
    const double ALPHA1 = .1817139853;
    const double DAL1 = .1257532909;
    const double BETA1 = 3.422509402;
    const double DG1 = .04742939676;
    const double P2 = -4.800458958;
    const double ALPHA2 = -.02845643596;
    const double DAL2 = .2188114228;
    const double BETA2 = 2.545944574;
    const double DG2 = .00813272793;
    const double BETA3 = .35868244;
    const double P3 = 103.1601001;
    const double ALPHA3 = -.00764731187;
    const double DAL3 = .1046487459;
    const double BETA4 = 2.958863546;
    const double DG3 = .01172314188;
    const double BETA5 = .4382872938;
    const double Q0 = .01134908150;
    const double Q1 = 14.51339943;
    const double ALPHA4 = .2647095287;
    const double DAL4 = .07091230197;
    const double DG4 = .01512963586;
    const double Q2 = 6.861329631;
    const double ALPHA5 = .1677400816;
    const double DAL5 = .04433648846;
    const double DG5 = .05553741389;
    const double BETA6 = .7665599464;
    const double BETA7 = .7277854652;

    int PROX = 0;
    double SINT1 = SINT;
    double COST1 = COST;
    if(SINT1 < 1.E-2){
        SINT1 = 1.E-2;
        COST1 = .99994999875;
        PROX = 1;
    }

    const double ALPHA = (SINT1*SINT1)/ R;
    const double GAMMA = COST1/pow(R,2);
    const double GAMMA_DG1 = GAMMA / DG1;
    const double GAMMA_DG4 = GAMMA / DG4;
    const double ALPHA_ALPHA4_DAL4 = (ALPHA - ALPHA4) / DAL4;
    const double ARG1 = -(GAMMA_DG1*GAMMA_DG1);
    const double ARG2 = -(ALPHA_ALPHA4_DAL4*ALPHA_ALPHA4_DAL4) - (GAMMA_DG4*GAMMA_DG4);

    double DEXP1 = 0.;
    if(ARG1 >= -500.) {
        DEXP1 = exp(ARG1);
    }

    double DEXP2 = 0.;
    if(ARG2 >= -500.) {
        DEXP2 = exp(ARG2);
    }

    const double ALPHA_S = ALPHA*(1. + P1/pow((1. + pow(((ALPHA-ALPHA1)/DAL1),2)),BETA1)*DEXP1 +
                                  P2*(ALPHA - ALPHA2)/pow((1. + pow(((ALPHA-ALPHA2)/DAL2),2)),BETA2)/pow((1. + pow((GAMMA/DG2),2)),BETA3) +
                                  P3*pow((ALPHA - ALPHA3),2)/pow((1. + pow(((ALPHA-ALPHA3)/DAL3),2)),BETA4)/pow((1. + pow((GAMMA/DG3),2)),BETA5));

    const double GAMMA_S = GAMMA*(1. + Q0+Q1*(ALPHA - ALPHA4)*DEXP2 +
                                  Q2*(ALPHA - ALPHA5)/pow((1. + pow(((ALPHA - ALPHA5)/DAL5),2)),BETA6)/pow((1. + pow((GAMMA/DG5),2)),BETA7));

    const double GAMMAS2 = GAMMA_S*GAMMA_S;

    const double ALSQH = (ALPHA_S*ALPHA_S)/2.;
    const double F = 64./27.*GAMMAS2 + ALSQH*ALSQH;
    const double Q = pow((sqrt(F) + ALSQH),(1./3.));
    double C = Q - 4.*pow(GAMMAS2,(1./3.))/(3.*Q);
    if(C < 0.) {
        C = 0.;
    }
    const double G = sqrt(C*C + 4.*pow(GAMMAS2,(1./3.)));
    const double RS = 4./((sqrt(2.*G - C) + sqrt(C))*(G + C));
    const double COSTS = GAMMA_S* (RS*RS);
    const double SINTS = sqrt(1. - COSTS*COSTS);
    const double RHOS = RS*SINTS;
    const double ZS = RS*COSTS;
    const double ZS2 = ZS*ZS;
    const double RRC1_RHOS = RRC1 + RHOS;

    double P = RRC1_RHOS*RRC1_RHOS + ZS2 + DD1*DD1;
    double XK2 = 4.*RRC1*RHOS/P;
    double XK = sqrt(XK2);
    double XKRHO12 = XK*sqrt(RHOS);

    double XK2S = 1. - XK2;
    double DL = log(1./XK2S);
    double ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383 + XK2S*(0.03742563713 + XK2S*0.01451196212))) +
                 DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    double ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
                 DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI1 = ((1. - XK2*0.5)*ELK - ELE)/XKRHO12;
    const double RRC2_RHOS = RRC2 + RHOS;

    P = RRC2_RHOS*RRC2_RHOS + ZS2 + DD2*DD2;
    XK2 = 4.*RRC2*RHOS/P;
    XK = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);

    XK2S = 1. - XK2;
    DL = log(1./XK2S);
    ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383 + XK2S*(0.03742563713 + XK2S*0.01451196212))) +
          DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
          DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI2 = ((1. - XK2*0.5)*ELK - ELE)/XKRHO12;
    double APPRC = A1*APHI1 + A2*APHI2;
    if(PROX) {
        APPRC = APPRC*SINT/SINT1;
    }
    return APPRC;
}

void prc_quad(double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double D = 1.E-4;
    const double DD = 2.E-4;
    const double DS = 1.E-2;
    const double DC = 0.99994999875;

    const double RHO2 = X*X + Y*Y;
    const double R = sqrt(RHO2 + Z*Z);
    const double RHO = sqrt(RHO2);
    const double SINT = RHO/R;
    const double COST = Z/R;
    const double RP = R + D;
    const double RM = R - D;

    if(SINT > DS) {
        const double CPHI = X/RHO;
        const double SPHI = Y/RHO;
        const double BR = br_prc_q(R, SINT, COST);
        const double BT = bt_prc_q(R, SINT, COST);
        const double DBRR = (br_prc_q(RP, SINT, COST) - br_prc_q(RM, SINT, COST))/DD;
        const double THETA = atan2(SINT, COST);
        const double TP = THETA+D;
        const double TM = THETA-D;
        const double SINTP = sin(TP);
        const double COSTP = cos(TP);
        const double SINTM = sin(TM);
        const double COSTM = cos(TM);
        const double DBTT = (bt_prc_q(R, SINTP, COSTP) - bt_prc_q(R, SINTM, COSTM))/DD;
        *BX = SINT*(BR + (BR + R*DBRR + DBTT)*(SPHI*SPHI)) + COST*BT;
        *BY = -SINT*SPHI*CPHI*(BR + R*DBRR + DBTT);
        *BZ = (BR*COST - BT*SINT)*CPHI;
    } else {
        const double ST = DS;
        double CT = DC;
        if(Z < 0.) {
            CT = -DC;
        }
        const double THETA = atan2(ST,CT);
        const double TP = THETA+D;
        const double TM = THETA-D;
        const double SINTP = sin(TP);
        const double COSTP = cos(TP);
        const double SINTM = sin(TM);
        const double COSTM = cos(TM);
        const double BR = br_prc_q(R, ST, CT);
        const double BT = bt_prc_q(R, ST, CT);
        const double DBRR = (br_prc_q(RP, ST, CT) - br_prc_q(RM, ST, CT))/DD;
        const double DBTT = (bt_prc_q(R, SINTP, COSTP) - bt_prc_q(R, SINTM, COSTM))/DD;
        const double FCXY = R*DBRR+DBTT;
        const double R_ST = R * ST;
        const double RST2 = R_ST*R_ST;
        const double Y2 = Y*Y;
        *BX = (BR*(X*X + 2.*Y2) + FCXY*Y2)/ RST2 + BT*COST;
        *BY = -(BR + FCXY)*X*Y/ RST2;
        *BZ = (BR*COST/ST - BT)*X/R;
    }
}

double br_prc_q(double R, double SINT, double COST) {
    const double A1 = -21.2666329;
    const double A2 = 32.24527521;
    const double A3 = -6.062894078;
    const double A4 = 7.515660734;
    const double A5 = 233.7341288;
    const double A6 = -227.1195714;
    const double A7 = 8.483233889;
    const double A8 = 16.80642754;
    const double A9 = -24.63534184;
    const double A10 = 9.067120578;
    const double A11 = -1.052686913;
    const double A12 = -12.08384538;
    const double A13 = 18.61969572;
    const double A14 = -12.71686069;
    const double A15 = 47017.35679;
    const double A16 = -50646.71204;
    const double A17 = 7746.058231;
    const double A18 = 1.531069371;
    const double XK1 = 2.318824273;
    const double AL1 = .1417519429;
    const double DAL1 = .6388013110E-02;
    const double B1 = 5.303934488;
    const double BE1 = 4.213397467;
    const double XK2 = .7955534018;
    const double AL2 = .1401142771;
    const double DAL2 = .2306094179E-01;
    const double B2 = 3.462235072;
    const double BE2 = 2.568743010;
    const double XK3 = 3.477425908;
    const double XK4 = 1.922155110;
    const double AL3 = .1485233485;
    const double DAL3 = .2319676273E-01;
    const double B3 = 7.830223587;
    const double BE3 = 8.492933868;
    const double AL4 = .1295221828;
    const double DAL4 = .01753008801;
    const double DG1 = .01125504083;
    const double AL5 = .1811846095;
    const double DAL5 = .04841237481;
    const double DG2 = .01981805097;
    const double C1 = 6.557801891;
    const double C2 = 6.348576071;
    const double C3 = 5.744436687;
    const double AL6 = .2265212965;
    const double DAL6 = .1301957209;
    const double DRM = .5654023158;

    const double SINT2 = SINT*SINT;
    const double COST2 = COST*COST;
    const double SC = SINT*COST;
    const double ALPHA = SINT2 / R;
    const double GAMMA = COST/ (R*R);

    double F = 0.;
    double FA = 0.;
    double FS = 0.;
    calc_ffs(ALPHA, AL1, DAL1, &F, &FA, &FS);
    const double D1 = SC*pow(F,XK1)/(pow((R/B1),BE1) + 1.);
    const double D2 = D1*COST2;

    calc_ffs(ALPHA, AL2, DAL2, &F, &FA, &FS);
    const double D3 = SC*pow(FS,XK2)/(pow((R/B2),BE2) + 1.);
    const double D4 = D3*COST2;

    calc_ffs(ALPHA, AL3, DAL3, &F, &FA, &FS);
    const double D5 = SC*(pow(ALPHA,XK3))*(pow(FS,XK4))/(pow((R/B3),BE3) + 1.);
    const double D6 = D5*COST2;
    const double ALPHA_AL4_DAL4 = (ALPHA - AL4) / DAL4;
    const double GAMMA_DG1 = GAMMA / DG1;

    double ARGA = ALPHA_AL4_DAL4*ALPHA_AL4_DAL4 + 1.;
    double ARGG = 1. + GAMMA_DG1*GAMMA_DG1;

    const double D7 = SC/ARGA/ARGG;
    const double D8 = D7/ARGA;
    const double D9 = D8/ARGA;
    const double D10 = D9/ARGA;
    const double ALPHA_AL5_DAL5 = (ALPHA - AL5) / DAL5;
    const double GAMMA_DG2 = GAMMA / DG2;

    ARGA = ALPHA_AL5_DAL5*ALPHA_AL5_DAL5 + 1.;
    ARGG = 1. + GAMMA_DG2*GAMMA_DG2;

    const double D11 = SC/ARGA/ARGG;
    const double D12 = D11/ARGA;
    const double D13 = D12/ARGA;
    const double D14 = D13/ARGA;

    const double R4 = R*R*R*R;

    const double D15 = SC/(R4 + C1*C1*C1*C1);
    const double D16 = SC/(R4 + C2*C2*C2*C2)*COST2;
    const double D17 = SC/(R4 + C3*C3*C3*C3)*(COST2*COST2);

    calc_ffs(ALPHA, AL6, DAL6, &F, &FA, &FS);
    const double R_DRM = (R - 1.2) / DRM;
    const double D18 = SC*FS/(1. + R_DRM*R_DRM);

    const double BR_PRC_Q = A1*D1 + A2*D2 + A3*D3 + A4*D4 + A5*D5 + A6*D6 + A7*D7 + A8*D8 + A9*D9 +
                            A10*D10 + A11*D11 + A12*D12 + A13*D13 + A14*D14 + A15*D15 + A16*D16 + A17*D17 + A18*D18;

    return BR_PRC_Q;
}

double bt_prc_q(double R, double SINT, double COST) {
    const double A1 = 12.74640393;
    const double A2 = -7.516393516;
    const double A3 = -5.476233865;
    const double A4 = 3.212704645;
    const double A5 = -59.10926169;
    const double A6 = 46.62198189;
    const double A7 = -.01644280062;
    const double A8 = .1234229112;
    const double A9 = -.08579198697;
    const double A10 = .01321366966;
    const double A11 = .8970494003;
    const double A12 = 9.136186247;
    const double A13 = -38.19301215;
    const double A14 = 21.73775846;
    const double A15 = -410.0783424;
    const double A16 = -69.90832690;
    const double A17 = -848.8543440;
    const double XK1 = 1.243288286;
    const double AL1 = .2071721360;
    const double DAL1 = .05030555417;
    const double B1 = 7.471332374;
    const double BE1 = 3.180533613;
    const double XK2 = 1.376743507;
    const double AL2 = .1568504222;
    const double DAL2 = .02092910682;
    const double BE2 = 1.985148197;
    const double XK3 = .3157139940;
    const double XK4 = 1.056309517;
    const double AL3 = .1701395257;
    const double DAL3 = .1019870070;
    const double B3 = 6.293740981;
    const double BE3 = 5.671824276;
    const double AL4 = .1280772299;
    const double DAL4 = .02189060799;
    const double DG1 = .01040696080;
    const double AL5 = .1648265607;
    const double DAL5 = .04701592613;
    const double DG2 = .01526400086;
    const double C1 = 12.88384229;
    const double C2 = 3.361775101;
    const double C3 = 23.44173897;

    const double SINT2 = SINT*SINT;
    const double COST2 = COST*COST;
    const double ALPHA = SINT2 / R;
    const double GAMMA = COST / (R*R);

    double F = 0.;
    double FA = 0.;
    double FS = 0.;
    calc_ffs(ALPHA, AL1, DAL1, &F, &FA, &FS);
    const double D1 = pow(F,XK1)/(pow((R/B1),BE1) + 1.);
    const double D2 = D1*COST2;

    calc_ffs(ALPHA, AL2, DAL2, &F, &FA, &FS);
    const double D3 = pow(FA,XK2)/pow(R,BE2);
    const double D4 = D3*COST2;

    calc_ffs(ALPHA, AL3, DAL3, &F, &FA, &FS);
    const double D5 = pow(FS,XK3)*pow(ALPHA,XK4)/(pow((R/B3),BE3) + 1.);
    const double D6 = D5*COST2;

    calc_ffs(GAMMA, 0., DG1, &F, &FA, &FS);
    const double ALPHA_AL4_DAL4 = (ALPHA - AL4) / DAL4;
    const double FCC = (1. + ALPHA_AL4_DAL4*ALPHA_AL4_DAL4);
    const double D7 = 1./FCC*FS;
    const double D8 = D7/FCC;
    const double D9 = D8/FCC;
    const double D10 = D9/FCC;
    const double ALPHA_AL5_DAL5 = (ALPHA - AL5) / DAL5;
    const double GAMMA_DG2 = GAMMA / DG2;

    const double ARG = 1. + ALPHA_AL5_DAL5*ALPHA_AL5_DAL5;
    const double D11 = 1./ARG/(1. + GAMMA_DG2*GAMMA_DG2);
    const double D12 = D11/ARG;
    const double D13 = D12/ARG;
    const double D14 = D13/ARG;

    const double R4 = R*R*R*R;

    const double D15 = 1./(R4 + C1*C1);
    const double D16 = COST2/(R4 + C2*C2);
    const double D17 = (COST2*COST2)/(R4 + C3*C3);

    const double BT_PRC_Q = A1*D1 + A2*D2 + A3*D3 + A4*D4 + A5*D5 + A6*D6 + A7*D7 + A8*D8 + A9*D9 +
                            A10*D10 + A11*D11 + A12*D12 + A13*D13 + A14*D14 + A15*D15 + A16*D16 + A17*D17;

    return BT_PRC_Q;
}

void calc_ffs(double A, double A0, double DA, double *F, double *FA, double *FS) {
    const double DA2 = DA*DA;
    const double A_A0 = A + A0;
    const double A_A0_DIF = A - A0;
    const double SQ1 = sqrt(A_A0*A_A0 + DA2);
    const double SQ2 = sqrt(A_A0_DIF*A_A0_DIF + DA2);
    *FA = 2./(SQ1+SQ2);
    *F = *FA * A;
    *FS = 0.5*(SQ1 + SQ2)/(SQ1*SQ2)*(1. - *F * *F);
}

void rc_shield(const double *A, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ) {
    const double X_SC1 = X_SC + 1.;
    const double FAC_SC = X_SC1 * X_SC1 * X_SC1;
    double GX = 0., GY = 0., GZ = 0.;
    int L = 0;

    for(int j = 0; j < 3; j++) {
        const double P = A[72 + j];
        const double YP = Y / P;
        const double SYPI = sin(YP);
        const double CYPI = cos(YP);
        const double P_SQPR = 1./ (P * P);

        for(int k = 0; k < 3; k++) {
            const double R = A[75 + k];
            const double Z1R = Z / R;
            const double SZRK = sin(Z1R);
            const double CZRK = cos(Z1R);
            const double SQPR = sqrt(P_SQPR + 1./ (R * R));
            const double EPR = exp(X*SQPR);

            GX += (-SQPR*EPR*CYPI*SZRK*FAC_SC) * (A[L] + X_SC*A[L + 1] + A[L + 2] + X_SC*A[L + 3]);
            GY += (EPR*SYPI*SZRK/P*FAC_SC) * (A[L] + X_SC*A[L + 1] + A[L + 2] + X_SC*A[L + 3]);
            GZ += (-EPR*CYPI*CZRK/R*FAC_SC) * (A[L] + X_SC*A[L + 1] + A[L + 2] + X_SC*A[L + 3]);
            L += 4;          
        }
    }

    *BX = GX;
    *BY = GY;
    *BZ = GZ;
}
