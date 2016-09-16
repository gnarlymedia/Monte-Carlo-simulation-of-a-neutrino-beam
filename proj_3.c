#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cpgplot.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

#define MAX_SIZE 10000000
#define ARRAY_SIZE 1000000
double xpi = 3.14159;
double c = 2.99792e8;
double pion_mass = 0.1396;
double pion_life = 26.08e-9;
double kaon_mass = 0.4937;
double kaon_life = 12.37e-9;
double muon_mass = 0.1057;
double K_BR = 0.63;
double pi_BR = 1.0;

double calc_beta(double mom, double mass);

double calc_gamma(double beta);

double calc_rad_pos(double t_mom, double l_mom, double decay_dist);

void calc_meson_mom_lab(int count_max, double *meson_decay_dist_array, double *meson_mom_array, double meson_mass,
                        double meson_neutr_mom_l_lab_array[], double meson_neutr_mom_t_lab_array[],
                        double valid_neutr_energy_lab_array[], int *count_valid_neutr, double decay_prob,
                        double valid_neutr_mom_l_lab[], double valid_neutr_mom_t_lab[], double valid_neutr_decay_dist[],
                        int *count_meson_neutr);

double dgauss(double x)
{
    double xpi = atan(1.0)*4.0;
    return exp(-x*x/(2.0))/sqrt(2.0*xpi);
}

double trnsfrm_stnd_nrml_to_nrml(double number, double mean, double width)
{
    return mean + width * number;
}

void box_muller(double distr_mean, double distr_width, double * value_1, double * value_2)
{
    int max_steps = 1000000;

    double x_1, x_2, u_1, u_2, d, y_1, y_2, sqrt_factor;

    int count;

    for (count = 1; count < max_steps; count++)
    {
        x_1 = drand48();
        x_2 = drand48();

        u_1 = 2.0 * x_1 - 1.0;
        u_2 = 2.0 * x_2 - 1.0;
        d = pow(u_1, 2.0) + pow(u_2, 2.0);

        if (d <= 1.0)
        {
            sqrt_factor = sqrt(-2.0 * log(d)/d);

            y_1 = u_1 * sqrt_factor;
            y_2 = u_2 * sqrt_factor;

            * value_1 = trnsfrm_stnd_nrml_to_nrml(y_1, distr_mean, distr_width);
            * value_2 = trnsfrm_stnd_nrml_to_nrml(y_2, distr_mean, distr_width);

            return;
        }
    }
}

double inverse_trans_meson_decay_time(double lifetime)
{
    return - lifetime * log(drand48());
}

double meson_distance_travelled(double mom, double mass, double time)
{
    return mom * c * time / mass;
}

double calculate_neut_mom_mag(double meson_mass)
{
    return (pow(meson_mass, 2.0) - pow(muon_mass, 2.0)) / (2.0 * meson_mass);
}

double calc_beta(double mom, double mass) {
    return mom / sqrt(pow(mom, 2.0) + pow(mass, 2.0));
}

void disp_histo(int num, double vals[], char * heading, char * x_label)
{
    int j = 0;
    static float fvals[ARRAY_SIZE+1];
    float v_min = 1.0e30;
    float v_max = -1.0e30;
    for(j = 0; j<num; j++)
    {
        fvals[j] = vals[j];
        if(fvals[j] < v_min)
        {
            v_min = fvals[j];
        }
        if(fvals[j] > v_max)
        {
            v_max = fvals[j];
        }
    }
    cpgbbuf();
    /*
     * Now plot a histogram
     */
    cpgsci(15);
    cpghist(num, fvals, v_min, v_max, 200,0);
//    printf("v_min: %f, v_max: %f", v_min, v_max);
    cpglab(x_label, "Counts", heading);
    // cpgsave saves the current graphics
    cpgsave();

}

void disp_histo_custom_v(int num, double vals[], char * heading, char * x_label, float v_min, float v_max)
{
    int j = 0;
    static float fvals[ARRAY_SIZE+1];
    for(j = 0; j < num; j++)
    {
        fvals[j] = vals[j];
    }
    cpgbbuf();

    /*
     * Now plot a histogram
     */
    cpgsci(15);
    cpghist(num, fvals, v_min, v_max, 200,0);
//    printf("v_min: %f, v_max: %f", v_min, v_max);
    cpglab(x_label, "Counts", heading);
    // cpgsave saves the current graphics
    cpgsave();
}

void display_scatter(int npoints, double xp[], double yp[],double xmin, double xmax, double ymin, double ymax, char * heading, char *x_label, char * y_label)
{
    float fxmin,fxmax,fymin,fymax;
    static float fxp[ARRAY_SIZE+1];
    static float fyp[ARRAY_SIZE+1];
    fxmin = (float) xmin;
    fxmax = (float) xmax;
    fymin = (float) ymin;
    fymax = (float) ymax;
    int j = 0;
    for( j=0; j<npoints; j++)
    {
        fxp[j] = (float) xp[j];
        fyp[j] = (float) yp[j];
    }
    //
    // Set up the dispay region
    //
    cpgenv(fxmin,fxmax,fymin,fymax,0,1);
    cpgbbuf();
    /*
     * Now make scatter plot with small dots plotted
     */
    cpgpt(npoints,fxp,fyp,-1);

    // Label the plot

    cpglab(x_label, y_label,heading);
    // cpgsave saves the current graphics
    cpgsave();

}

void display_scatter_two_arrays(int npoints, double xp[], double yp[], double xp2[], double yp2[], double xmin, double xmax, double ymin, double ymax, char * heading, char *x_label, char * y_label)
{
    float fxmin,fxmax,fymin,fymax;
    static float fxp[ARRAY_SIZE+1];
    static float fyp[ARRAY_SIZE+1];
    static float fxp2[ARRAY_SIZE+1];
    static float fyp2[ARRAY_SIZE+1];
    fxmin = (float) xmin;
    fxmax = (float) xmax;
    fymin = (float) ymin;
    fymax = (float) ymax;
    int j = 0;
    for(j=0; j<npoints; j++)
    {
        fxp[j] = (float) xp[j];
        fyp[j] = (float) yp[j];
    }
    int k = 0;
    for(k=0; k<npoints; k++)
    {
        fxp2[j] = (float) xp2[j];
        fyp2[j] = (float) yp2[j];
    }
    //
    // Set up the dispay region
    //
    cpgenv(fxmin,fxmax,fymin,fymax,0,1);
    cpgbbuf();
    /*
     * Now make scatter plot with small dots plotted
     */
    cpgpt(npoints,fxp,fyp,-1);
    cpgpt(npoints,fxp2,fyp2,-1);

    // Label the plot

    cpglab(x_label, y_label,heading);
    // cpgsave saves the current graphics
    cpgsave();

}

int main(void)
{
    int i = 0;
    i = (int) time(NULL);
    //
    // Set a unique random seed
    //
    srand48(i);

    int npoints = ARRAY_SIZE/10;
    static double duni[ARRAY_SIZE+1];
    static double duni2[ARRAY_SIZE+1];
    static double dquad[ARRAY_SIZE+1];
    static double dgau[ARRAY_SIZE+1];
    static double dgauss_mean_arr[ARRAY_SIZE+1];
    static double dgauss_quad_arr[ARRAY_SIZE+1];

    /*
     * Call ppgbeg to initiate PGPLOT and open the output device; cpgbeg
     * will prompt the user to supply the device name and type.
     */
//    if (1 != cpgbeg(0, "?", 1, 1))
    if (1 != cpgbeg(0, "/XWINDOW", 1, 1))
//    if (1 != cpgbeg(0, "proj_3_plot.ps/VCPS", 1, 1))
//    if (1 != cpgbeg(0, "proj_3_plot.ps/CPS", 1, 1))
    {
        exit(EXIT_FAILURE);
    }
    cpgask(1);

    int count_pions = 1;
    int count_kaons = 1;

    // Box-Buller
    int count_meson_mom;
    static double pion_mom_array[ARRAY_SIZE];
    static double kaon_mom_array[ARRAY_SIZE];
    double stnd_val_1, stnd_val_2;

    for (count_meson_mom = 1; count_meson_mom < ARRAY_SIZE; count_meson_mom++)
    {
        box_muller(200, 10, &stnd_val_1, &stnd_val_2);
        // select which type of meson we have
        if (drand48() <+ 0.86) {
            // pion
            pion_mom_array[count_pions] = stnd_val_1;
            count_pions++;
        }
        else {
            // kaon
            kaon_mom_array[count_kaons] = stnd_val_2;
            count_kaons++;
        }
    }

    // plot meson momenta
//    disp_histo(count_pions, pion_mom_array, "Pion momentum distribution", "Momentum");
//    disp_histo(count_kaons, kaon_mom_array, "Kaon momentum distribution", "Momentum");


    // decay distances
    int count_decay_distances;
    static double pion_decay_dist_array[ARRAY_SIZE];
    static double kaon_decay_dist_array[ARRAY_SIZE];

    for (count_decay_distances = 1; count_decay_distances < count_pions; count_decay_distances++)
    {
        pion_decay_dist_array[count_decay_distances] = meson_distance_travelled(pion_mom_array[count_decay_distances], pion_mass, inverse_trans_meson_decay_time(pion_life));
    }

    for (count_decay_distances = 1; count_decay_distances < count_kaons; count_decay_distances++)
    {
        kaon_decay_dist_array[count_decay_distances] = meson_distance_travelled(kaon_mom_array[count_decay_distances], kaon_mass, inverse_trans_meson_decay_time(kaon_life));
    }

//    disp_histo_custom_v(count_decay_distances, pion_decay_dist_array, "Pion decay distances", "Decay distance", 0.0, 10000.0);
//    disp_histo_custom_v(count_decay_distances, kaon_decay_dist_array, "Kaon decay distances", "Decay distance", 0.0, 10000.0);

    // simulate neutrino direction and calculate momentum distribution of pion neutrinos in lab frame
    double neutr_cos_theta_rest, neutr_mom_l_rest, neutr_mom_t_rest, neutr_energy_rest;
    double beta, gamma_var;
    int count_neutr_mom;
    int count_valid_neutr = 1;
    static double valid_neutr_energy_lab_array[ARRAY_SIZE];
    static double valid_neutr_mom_l_lab_array[ARRAY_SIZE];
    static double valid_neutr_mom_t_lab_array[ARRAY_SIZE];
    static double valid_neutr_decay_dist[ARRAY_SIZE];

    // pions
    static double pion_neutr_mom_l_lab_array[ARRAY_SIZE];
    static double pion_neutr_mom_t_lab_array[ARRAY_SIZE];
    int count_pion_neutr = 1;

    calc_meson_mom_lab(count_pions, pion_decay_dist_array, pion_mom_array, pion_mass, pion_neutr_mom_l_lab_array,
                       pion_neutr_mom_t_lab_array, valid_neutr_energy_lab_array, &count_valid_neutr, 1.0,
                       valid_neutr_mom_l_lab_array, valid_neutr_mom_t_lab_array, valid_neutr_decay_dist, &count_pion_neutr);

//    disp_histo(count_pion_neutr, pion_neutr_mom_l_lab_array, "Pion neutrino longitudinal momentum distribution", "Momentum");
//    disp_histo(count_pion_neutr, pion_neutr_mom_t_lab_array, "Pion neutrino transverse momentum distribution", "Momentum");

    // kaons
    static double kaon_neutr_mom_l_lab_array[ARRAY_SIZE];
    static double kaon_neutr_mom_t_lab_array[ARRAY_SIZE];
    int count_kaon_neutr = 1;

    calc_meson_mom_lab(count_kaons, kaon_decay_dist_array, kaon_mom_array, kaon_mass, kaon_neutr_mom_l_lab_array,
                       kaon_neutr_mom_t_lab_array, valid_neutr_energy_lab_array, &count_valid_neutr, 0.64,
                       valid_neutr_mom_l_lab_array, valid_neutr_mom_t_lab_array, valid_neutr_decay_dist, &count_kaon_neutr);

//    disp_histo(count_kaon_neutr, kaon_neutr_mom_l_lab_array, "Kaon neutrino longitudinal momentum distribution", "Momentum");
//    disp_histo(count_kaon_neutr, kaon_neutr_mom_t_lab_array, "Kaon neutrino transverse momentum distribution", "Momentum");

    int count_neutr_radial_pos;
    static double on_target_neutr_energy_array[ARRAY_SIZE];
    static double on_target_neutr_rad_pos_array[ARRAY_SIZE];
    double rad_pos;
    int count_on_target_rad_pos = 1;

    for (count_neutr_radial_pos = 1; count_neutr_radial_pos < count_valid_neutr; count_neutr_radial_pos++)
    {
        rad_pos = calc_rad_pos(valid_neutr_mom_t_lab_array[count_neutr_radial_pos], valid_neutr_mom_l_lab_array[count_neutr_radial_pos], valid_neutr_decay_dist[count_neutr_radial_pos]);
        if (rad_pos <= 1.5)
        {
            on_target_neutr_rad_pos_array[count_on_target_rad_pos] = rad_pos;
            on_target_neutr_energy_array[count_on_target_rad_pos] = valid_neutr_energy_lab_array[count_neutr_radial_pos];
            count_on_target_rad_pos++;
        }
    }

    display_scatter(count_on_target_rad_pos, on_target_neutr_rad_pos_array, on_target_neutr_energy_array, 0.0, 1.5, 0, 250.0, "Neutrino energy (lab frame) vs detector radial position", "Radial position (m)", "Neutrino energy (GeV)");

    cpgend();
}

void calc_meson_mom_lab(int count_max, double *meson_decay_dist_array, double *meson_mom_array, double meson_mass,
                        double meson_neutr_mom_l_lab_array[], double meson_neutr_mom_t_lab_array[],
                        double valid_neutr_energy_lab_array[], int *count_valid_neutr, double decay_prob,
                        double valid_neutr_mom_l_lab[], double valid_neutr_mom_t_lab[], double valid_neutr_decay_dist[],
                        int *count_meson_neutr) {
    // calculate momenta of neutrinos in lab frame
    int count_neutr_mom;
    double beta, gamma_var, neutr_cos_theta_rest, neutr_mom_l_rest, neutr_mom_t_rest, neutr_energy_rest, neutr_mom_l_lab, neutr_mom_t_lab;
    double neutr_mom_mag_rest = calculate_neut_mom_mag(meson_mass);
    for (count_neutr_mom = 1; count_neutr_mom < count_max; count_neutr_mom++)
    {
        if (meson_decay_dist_array[count_neutr_mom] <= 300.0)
        {
            if (drand48() < decay_prob) {
                beta = calc_beta(meson_mom_array[count_neutr_mom], meson_mass);
                gamma_var = calc_gamma(beta);

                // calculate directions
                neutr_cos_theta_rest = drand48() * 2.0 - 1.0;
                neutr_mom_l_rest = neutr_mom_mag_rest * neutr_cos_theta_rest;
                neutr_mom_t_rest = neutr_mom_mag_rest * sin(acos(neutr_cos_theta_rest));

                // transform to lab frame
                neutr_energy_rest = neutr_mom_mag_rest;
                neutr_mom_l_lab = gamma_var * neutr_mom_l_rest + beta * gamma_var * neutr_energy_rest;
                neutr_mom_t_lab = neutr_mom_t_rest;

                meson_neutr_mom_l_lab_array[*count_meson_neutr] = neutr_mom_l_lab;
                meson_neutr_mom_t_lab_array[*count_meson_neutr] = neutr_mom_t_lab;

                valid_neutr_energy_lab_array[*count_valid_neutr] = beta * gamma_var * neutr_mom_l_rest + gamma_var * neutr_energy_rest;
                valid_neutr_mom_l_lab[*count_valid_neutr] = neutr_mom_l_lab;
                valid_neutr_mom_t_lab[*count_valid_neutr] = neutr_mom_t_lab;
                valid_neutr_decay_dist[*count_valid_neutr] = meson_decay_dist_array[count_neutr_mom];

                *count_valid_neutr = *count_valid_neutr + 1;
                *count_meson_neutr = *count_meson_neutr + 1;
            }
        }
    }
}

double calc_rad_pos(double t_mom, double l_mom, double decay_dist) {
    return (t_mom / l_mom) * (700.0 - decay_dist);
}

double calc_gamma(double beta) {
    return 1.0 / sqrt(1.0 - pow(beta, 2.0));
}
