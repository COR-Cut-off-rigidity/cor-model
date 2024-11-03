#include "trajectory.h"
#include "utils.h"
#include <math.h>
#include <omp.h>

void trajectory_simulation_cutoff(InfileHeader *header, const Model **models, ModelsParams *models_params, FILE *outfile, uint64_t step_limit, double tu_angle, double particle_mass, uint32_t parallel) {
    const int rig_values_size = (int) round((header->rk - header->rig) / header->del);
    int nza = 0;
    double rmx1 = 0., rmx2 = 0., rmi = header->rk;

    #pragma omp parallel for ordered schedule(dynamic, 1) reduction(+:nza) if (parallel)
    for(int i = 0; i <= rig_values_size; i++) {
        int nk0 = header->nk1;
        double thread_rig = header->rig + i * header->del;
        double alength = 0., time = 0.;
        double a = thread_rig*1.0e+09*header->zn*DEF_Q / (particle_mass*DEF_C*DEF_C);

        Vector sph_vel = {
            .r = DEF_C * sqrt(1. - 1./(1. + a*a)),
            .theta = (90. - header->the1) * DEF_TO_RAD,
            .phi = header->fi1 * DEF_TO_RAD
        };

        Vector sph_pos = {
            .r = header->r0 * DEF_RE,
            .theta = (90. - header->the0) * DEF_TO_RAD,
            .phi = header->fi0 * DEF_TO_RAD
        };

        Vector car_pos;
        sph_to_car(&sph_pos, &car_pos);

        Vector gsm_pos;
        geo_to_gsm(&car_pos, &gsm_pos, models_params->A1, models_params->A2, models_params->A3);

        Vector car_field;
        calculate_field(gsm_pos, &car_field, models, models_params);
        double b = 1.E-09*sqrt(car_field.x*car_field.x + car_field.y*car_field.y + car_field.z*car_field.z);

        Vector car_vel;
        sph_to_car(&sph_vel, &car_vel);
        double v = sph_vel.r;
        unsigned int step = 0;

        do {
            step++;

            double vc = v / DEF_C;
            if(vc > 1.){
                nza += 1;
                #pragma omp atomic write
                rmx1 = fmax(rmx1, thread_rig);
                break;
            }

            const double hm = particle_mass / sqrt(1. - vc*vc);
            const double t = (2.*DEF_PI*hm) / (fabs(header->zn)*DEF_Q*b);
            const double h = t/nk0;

            double d = (1.E-09*header->zn*DEF_Q) / hm;
            Vector prev_car_pos = car_pos;
            Vector prev_car_vel = car_vel;
            sph_vel.r = v;

            calculate_position_and_speed(h, d, car_field, &car_vel, &car_pos);
            v = sqrt(car_vel.x*car_vel.x + car_vel.y*car_vel.y + car_vel.z*car_vel.z);

            const double gc = (car_vel.x*prev_car_vel.x + car_vel.y*prev_car_vel.y + car_vel.z*prev_car_vel.z)/(v*sph_vel.r);
            if(gc < 1.) {
                const double gs = sqrt((1.-gc) * (1.+gc));
                const double gu = atan2(gs, gc);
                if(gu >= tu_angle && nk0 <= 50000) {
                    nk0 *= 2;
                    car_pos = prev_car_pos;
                    car_vel = prev_car_vel;
                    step--;
                    continue;
                }
            }

            car_to_sph(&car_pos, &sph_pos);
            geo_to_gsm(&car_pos, &gsm_pos, models_params->A1, models_params->A2, models_params->A3);

            const double pax = (gsm_pos.x/DEF_RE + 25.3)/36.08;
            const double paz = (gsm_pos.y*gsm_pos.y + gsm_pos.z*gsm_pos.z)/(DEF_RE*DEF_RE);
            const double r = sph_pos.r/DEF_RE;
            const double paus = pax*pax + paz/459.6736;

            if(r >= 25. || paus > 1.) {
                if((r >= 25. && fabs(r - 25.) < 0.002) || (r <= 25. && paus > 1. && fabs(paus - 1.) < 0.002)) {
                    const double td = 90. - (sph_pos.theta / DEF_TO_RAD);
                    const double fei = sph_pos.phi / DEF_TO_RAD;
                    
                    const double s = sin(sph_pos.theta);
                    const double c = cos(sph_pos.theta);
                    const double sf = sin(sph_pos.phi);
                    const double cf = cos(sph_pos.phi);
                    const double be = car_vel.x*cf + car_vel.y*sf;
                    Vector v_sp = {
                        .r = be*s + car_vel.z*c,
                        .theta = be*c - car_vel.z*s,
                        .phi = -car_vel.x*sf + car_vel.y*cf
                    };

                    const double va = v_sp.theta*c + v_sp.r*s;
                    const double dm = sqrt(v_sp.theta*v_sp.theta + va*va);
                    double ast = atan2(-v_sp.theta*s + v_sp.r*c, dm) / DEF_TO_RAD;
                    double asf = (sph_pos.phi + atan2(v_sp.phi, va)) / DEF_TO_RAD;

                    if(asf > 360.) {
                        asf -= 360.;
                    }

                    #pragma omp atomic write
                    rmi = fmin(rmi, thread_rig);
                    
                    #pragma omp ordered
                    fprintf(
                        outfile,
                        "%10.6f%15.10f%12.6f%10.3f%10.3f%10.3f%10.3f%12.6f%16.2f\n",
                        thread_rig, v / DEF_C, r, td, fei, ast, asf, time, alength / 1000.
                    );
                    break;
                }

                nk0 *= 2;
                car_pos = prev_car_pos;
                car_vel = prev_car_vel;
                step--;
                continue;
            }

            if(r < 1.) {
                nza += 1;
                #pragma omp atomic write
                rmx1 = fmax(rmx1, thread_rig);
                break;
            }

            if(step > step_limit || sph_pos.phi > 31.4) {
                nza += 1;
                #pragma omp atomic write
                rmx2 = fmax(rmx2, thread_rig);
                break;
            }

            calculate_field(gsm_pos, &car_field, models, models_params);
            b = 1.E-09*sqrt(car_field.x*car_field.x + car_field.y*car_field.y + car_field.z*car_field.z);

            time = time + t/nk0;
            alength += sqrt((car_pos.x-prev_car_pos.x)*(car_pos.x-prev_car_pos.x) + (car_pos.y-prev_car_pos.y)*(car_pos.y-prev_car_pos.y) + (car_pos.z-prev_car_pos.z)*(car_pos.z-prev_car_pos.z));
            nk0 = header->nk1;
        } while(1);
    }

    double rmx = fmax(rmx1, rmx2);
    const double zan = nza - (rmi - header->rig) / header->del;
    double rms = rmi + zan*header->del;
    rmx = rmx + header->del;
    fprintf(outfile, "  CUTOFF with rigidities P(S),P(C),P(M) are:\n%12.5f%12.5f%12.5f\n\n", rmi, rmx, rms);
}

void trajectory_simulation(InfileHeader *header, const Model **models, ModelsParams *models_params, FILE *tragsm, FILE *trasph, uint64_t step_limit, double tu_angle, double particle_mass) {
    const int rig_values_size = (int) round((header->rk - header->rig) / header->del);

    for(int i = 0; i <= rig_values_size; i++) {
        int nk0 = header->nk1;
        double thread_rig = header->rig + i * header->del;
        double alength = 0., time = 0.;
        double a = thread_rig*1.0e+09*header->zn*DEF_Q / (particle_mass*DEF_C*DEF_C);

        Vector sph_vel = {
            .r = DEF_C * sqrt(1. - 1./(1. + a*a)),
            .theta = (90. - header->the1) * DEF_TO_RAD,
            .phi = header->fi1 * DEF_TO_RAD
        };

        Vector sph_pos = {
            .r = header->r0 * DEF_RE,
            .theta = (90. - header->the0) * DEF_TO_RAD,
            .phi = header->fi0 * DEF_TO_RAD
        };

        Vector car_pos;
        sph_to_car(&sph_pos, &car_pos);

        Vector gsm_pos;
        geo_to_gsm(&car_pos, &gsm_pos, models_params->A1, models_params->A2, models_params->A3);

        Vector car_field;
        calculate_field(gsm_pos, &car_field, models, models_params);
        double b = 1.E-09*sqrt(car_field.x*car_field.x + car_field.y*car_field.y + car_field.z*car_field.z);

        Vector car_vel;
        sph_to_car(&sph_vel, &car_vel);
        double v = sph_vel.r;
        unsigned int step = 0;

        do {
            step++;

            double vc = v / DEF_C;
            if(vc > 1.) {
                break;
            }

            const double hm = particle_mass / sqrt(1. - vc*vc);
            const double t = (2.*DEF_PI*hm) / (fabs(header->zn)*DEF_Q*b);
            const double h = t/nk0;

            double d = (1.E-09*header->zn*DEF_Q) / hm;
            Vector prev_pos = car_pos;
            Vector prev_vel = car_vel;
            sph_vel.r = v;

            calculate_position_and_speed(h, d, car_field, &car_vel, &car_pos);
            v = sqrt(car_vel.x*car_vel.x + car_vel.y*car_vel.y + car_vel.z*car_vel.z);

            const double gc = (car_vel.x*prev_vel.x + car_vel.y*prev_vel.y + car_vel.z*prev_vel.z)/(v*sph_vel.r);
            if(gc < 1.) {
                const double gs = sqrt((1.-gc) * (1.+gc));
                const double gu = atan2(gs, gc);
                if(gu >= tu_angle && nk0 <= 50000) {
                    nk0 *= 2;
                    car_pos = prev_pos;
                    car_vel = prev_vel;
                    step--;
                    continue;
                }
            }

            car_to_sph(&car_pos, &sph_pos);
            geo_to_gsm(&car_pos, &gsm_pos, models_params->A1, models_params->A2, models_params->A3);

            const double pax = (gsm_pos.x/DEF_RE + 25.3)/36.08;
            const double paz = (gsm_pos.y*gsm_pos.y + gsm_pos.z*gsm_pos.z)/(DEF_RE*DEF_RE);
            const double r = sph_pos.r/DEF_RE;
            const double paus = pax*pax + paz/459.6736;
            const double td = 90. - (sph_pos.theta / DEF_TO_RAD);
            const double fei = sph_pos.phi / DEF_TO_RAD;

            if(r >= 25. || paus > 1.) {
                if((r >= 25. && fabs(r - 25.) < 0.002) || (r <= 25. && paus > 1. && fabs(paus - 1.) < 0.002)) {
                    break;
                }

                nk0 *= 2;
                car_pos = prev_pos;
                car_vel = prev_vel;
                step--;
                continue;
            }

            if(r < 1. || step > step_limit || sph_pos.phi > 31.4) {
                break;
            }

            calculate_field(gsm_pos, &car_field, models, models_params);
            b = 1.E-09*sqrt(car_field.x*car_field.x + car_field.y*car_field.y + car_field.z*car_field.z);

            time = time + t/nk0;
            alength += sqrt((car_pos.x-prev_pos.x)*(car_pos.x-prev_pos.x) + (car_pos.y-prev_pos.y)*(car_pos.y-prev_pos.y) + (car_pos.z-prev_pos.z)*(car_pos.z-prev_pos.z));
            
            #pragma omp ordered            
            fprintf(tragsm, "%6d%17.5f%17.5f%17.5f%14.2f%15.10f%18.2f\n", step, gsm_pos.x, gsm_pos.y, gsm_pos.z, v, time, alength);

            #pragma omp ordered
            fprintf(trasph, "%6d%17.5f%17.5f%17.5f%14.2f%15.10f%18.2f\n", step, r, td, fei, v, time, alength);
            
            nk0 = header->nk1;
        } while(1);
    }
}

void calculate_position_and_speed(double h, double d, Vector field, Vector *vel, Vector *pos) {
    Vector v = *vel;
    double a1 = h * (d * (v.y*field.z - v.z*field.y));
    double b1 = h * (d * (v.z*field.x - v.x*field.z));
    double c1 = h * (d * (v.x*field.y - v.y*field.x));

    const double y_a2 = v.y + b1/3.;
    const double z_a2 = v.z + c1/3.;
    const double a2 = h*(d*(y_a2*field.z - z_a2*field.y));
    const double y_b2 = v.x + a1/3.;
    const double z_b2 = v.z + c1/3.;
    const double b2 = h*(d*(z_b2*field.x - y_b2*field.z));
    const double y_c2 = v.x + a1/3.;
    const double z_c2 = v.y + b1/3.;
    const double c2 = h*(d*(y_c2*field.y - z_c2*field.x));

    const double y_a3 = v.y + (4.*b1 + 6.*b2)/25.;
    const double z_a3 = v.z + (4.*c1 + 6.*c2)/25.;
    const double a3 = h*(d*(y_a3*field.z - z_a3*field.y));
    const double y_b3 = v.x + (4.*a1 + 6.*a2)/25.;
    const double z_b3 = v.z + (4.*c1 + 6.*c2)/25.;
    const double b3 = h*(d*(z_b3*field.x - y_b3*field.z));
    const double y_c3 = v.x + (4.*a1 + 6.*a2)/25.;
    const double z_c3 = v.y + (4.*b1 + 6.*b2)/25.;
    const double c3 = h*(d*(y_c3*field.y - z_c3*field.x));

    const double y_a4 = v.y + (b1 - 12.*b2+15.*b3)/4.;
    const double z_a4 = v.z + (c1 - 12.*c2 + 15.*c3)/4.;
    const double a4 = h*(d*(y_a4*field.z - z_a4*field.y));
    const double y_b4 = v.x + (a1 - 12.*a2+15.*a3)/4.;
    const double z_b4 = v.z + (c1 - 12.*c2 + 15.*c3)/4.;
    const double b4 = h*(d*(z_b4*field.x - y_b4*field.z));
    const double y_c4 = v.x + (a1 - 12.*a2+15.*a3)/4.;
    const double z_c4 = v.y + (b1 - 12.*b2 + 15.*b3)/4.;
    const double c4 = h*(d*(y_c4*field.y - z_c4*field.x));

    const double y_a5 = v.y + (6.*b1 + 90.*b2 - 50.*b3 + 8.*b4)/81.;
    const double z_a5 = v.z + (6.*c1 + 90.*c2 - 50.*c3 + 8.*c4)/81.;
    const double a5 = h*(d*(y_a5*field.z - z_a5*field.y));
    const double y_b5 = v.x + (6.*a1 + 90.*a2 - 50.*a3 + 8.*a4)/81.;
    const double z_b5 = v.z + (6.*c1 + 90.*c2 - 50.*c3 + 8.*c4)/81.;
    const double b5 = h*(d*(z_b5*field.x - y_b5*field.z));
    const double y_c5 = v.x + (6.*a1 + 90.*a2 - 50.*a3 + 8.*a4)/81.;
    const double z_c5 = v.y + (6.*b1 + 90.*b2 - 50.*b3 + 8.*b4)/81.;
    const double c5 = h*(d*(y_c5*field.y - z_c5*field.x));

    const double y_a6 = v.y + (6.*b1 + 36.*b2 + 10.*b3 + 8.*b4)/75.;
    const double z_a6 = v.z + (6.*c1 + 36.*c2 + 10.*c3 + 8.*c4)/75.;
    const double a6 = h*(d*(y_a6*field.z - z_a6*field.y));
    const double y_b6 = v.x + (6.*a1 + 36.*a2 + 10.*a3 + 8.*a4)/75.;
    const double z_b6 = v.z + (6.*c1 + 36.*c2 + 10.*c3 + 8.*c4)/75.;
    const double b6 = h*(d*(z_b6*field.x - y_b6*field.z));
    const double y_c6 = v.x + (6.*a1 + 36.*a2 + 10.*a3 + 8.*a4)/75.;
    const double z_c6 = v.y + (6.*b1 + 36.*b2 + 10.*b3 + 8.*b4)/75.;
    const double c6 = h*(d*(y_c6*field.y - z_c6*field.x));

    const double aa1 = h*v.x;
    const double bb1 = h*v.y;
    const double cc1 = h*v.z;
    const double aa2 = h*(v.x + aa1/3.);
    const double bb2 = h*(v.y + bb1/3.);
    const double cc2 = h*(v.z + cc1/3.);
    const double aa3 = h*(v.x + (4.*aa1 + 6.*aa2)/25.);
    const double bb3 = h*(v.y + (4.*bb1 + 6.*bb2)/25.);
    const double cc3 = h*(v.z + (4.*cc1 + 6.*cc2)/25.);
    const double aa4 = h*(v.x + (aa1 - 12.*aa2 + 15.*aa3)/4.);
    const double bb4 = h*(v.y + (bb1 - 12.*bb2 + 15.*bb3)/4.);
    const double cc4 = h*(v.z + (cc1 - 12.*cc2 + 15.*cc3)/4.);
    const double aa5 = h*(v.x + (6.*aa1 + 90.*aa2 - 50.*aa3 + 8.*aa4)/81.);
    const double bb5 = h*(v.y + (6.*bb1 + 90.*bb2 - 50.*bb3 + 8.*bb4)/81.);
    const double cc5 = h*(v.z + (6.*cc1 + 90.*cc2 - 50.*cc3 + 8.*cc4)/81.);
    const double aa6 = h*(v.x + (6.*aa1 + 36.*aa2 + 10.*aa3 + 8.*aa4)/75.);
    const double bb6 = h*(v.y + (6.*bb1 + 36.*bb2 + 10.*bb3 + 8.*bb4)/75.);
    const double cc6 = h*(v.z + (6.*cc1 + 36.*cc2 + 10.*cc3 + 8.*cc4)/75.);

    vel->x += (23.*a1 + 125.*a3 - 81.*a5 + 125.*a6)/192.;
    vel->y += (23.*b1 + 125.*b3 - 81.*b5 + 125.*b6)/192.;
    vel->z += (23.*c1 + 125.*c3 - 81.*c5 + 125.*c6)/192.;
    pos->x += (23.*aa1 + 125.*aa3 - 81.*aa5 + 125.*aa6)/192.;
    pos->y += (23.*bb1 + 125.*bb3 - 81.*bb5 + 125.*bb6)/192.;
    pos->z += (23.*cc1 + 125.*cc3 - 81.*cc5 + 125.*cc6)/192.;
}
