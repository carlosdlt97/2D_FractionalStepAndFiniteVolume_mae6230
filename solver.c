#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

int main()
{
    /* ----------------------------------------------------------------------------------------------------------------------------
    Initializing Variables  -------------------------------------------------------------------------------------------------------
    */
    int i = 0; /* Counters for all "for loops" */
    int j = 0;
    int k = 0;

    double Re = 100; /* Problem parameters */
    double D_t = 0.001;
    int nodes_x = 21;
    int nodes_y = 21;
    double NX = nodes_x;
    double NY = nodes_y;
    double D_x = 1 / (NX - 1);
    double D_y = 1 / (NY - 1);
    double lambda = pow(D_x, -2);
    double f_norm;


    double** u = (double**)calloc(nodes_x, sizeof(double*));    /* Memory allocation for large arrays (velocities, etc.) */
    double** v = (double**)calloc(nodes_y, sizeof(double*));    /* u and v represent barycentric velocities */

    double** u_star = (double**)calloc(nodes_x, sizeof(double*));
    double** v_star = (double**)calloc(nodes_y, sizeof(double*));
    for (i = 0; i < nodes_x; i++) {
        u[i] = (double*)calloc(nodes_y, sizeof(double));
        u_star[i] = (double*)calloc(nodes_y, sizeof(double));
    }
    for (i = 0; i < nodes_y; i++)
    {
        v[i] = (double*)calloc(nodes_x, sizeof(double));
        v_star[i] = (double*)calloc(nodes_x, sizeof(double));
    }

    double** p = (double**)calloc(nodes_x, sizeof(double*));
    double** p_new = (double**)calloc(nodes_x, sizeof(double*));
    double** f = (double**)calloc(nodes_x, sizeof(double*));
    for (i = 0; i < nodes_x; i++)
    {
        p[i] = (double*)calloc(nodes_y, sizeof(double));
        p_new[i] = (double*)calloc(nodes_y, sizeof(double));
        f[i] = (double*)calloc(nodes_y, sizeof(double));
    }

    /* ----------------------------------------------------------------------------------------------------------------------------
    Step 1. -----------------------------------------------------------------------------------------------------------------------
    */

    /* Initializing varibales for this step */
    double Hx;
    double Hy;
    double u_cc; /* cc = Cell cented velocity */
    double v_cc;
    double u_cc_im1; /* im1 = "i minus 1"; cell centered velocity at 1 position back in the i direction */
    double v_cc_jm1; /* jm1 = "j minus 1"; cell centered velocity at 1 position back in the j direction */
    double u_s;  /* s = staggered velocity */
    double v_s;
    double u_s_ip1; /* ip1 = "i plus 1"; staggered velocity at 1 position forward in the i direction */
    double u_s_jp1; /* jp1 = "j plus 1"; staggered velocity at 1 position forward in the j direction */
    double v_s_ip1;
    double v_s_jp1;
    for (j = 0; j < nodes_y; j++) {
        for (i = 0; i < nodes_x; i++) {
            /* Calculating u_star */
            if (i > 0 && j > 0 && i < (nodes_x - 1) && j < (nodes_y - 2)) {  /* At interior points except at the edges */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / D_x) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;

                u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_x, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(D_y, 2)) * u[j][i]) + u[j][i];
            }
            else if (i == 0) {      /* At the left wall */
                u_star[j][i] = 0;
            }
            else if (i == nodes_x - 1) {     /* At the right wall */
                u_star[j][i] = 0;
            }
            else if (j == 0) {      /* Just above the bottom wall */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = 0;
                v_s = 0;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i - 1] + v[j + 1][i]) / 2;

                Hx = (pow(u_cc, 2) + pow(u_cc_im1, 2)) / D_x + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;

                u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_x, 2) + ((1 / D_y) * (u[j + 1][i] - u[j][i]) + (2 / D_y) * (-u[j][i])) / D_y) * u[j][i]) + u[j][i];
            }
            else if (j == nodes_y - 2) {    /* Just below the top lid */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i - 1] + v[j][i]) / 2;
                u_s_jp1 = 1;
                v_s_jp1 = 0;

                Hx = (pow(u_cc, 2) + pow(u_cc_im1, 2)) / D_x + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;

                u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_x, 2) + ((2 / D_y) * (1 - u[j][i]) + (1 / D_y) * (u[j - 1][i] - u[j][i])) / D_y) * u[j][i]) + u[j][i];
            }
            else if (j == nodes_y - 1) {     /* At the top lid */
                u_star[j][i] = 1;
            }
            /* Calculating v_star */
            if (i > 0 && j > 0 && i < (nodes_x - 2) && j < (nodes_y - 1)) {  /* At interior points except at the edges */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / D_y) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / D_x;

                v_star[j][i] = D_t * (Hy + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(D_x, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_y, 2)) * v[j][i]) + v[j][i];
            }
            else if (j == 0) {    /* At the bottom wall */
                v_star[j][i] = 0;
            }
            else if (j == nodes_y - 1) {     /* At the top lid */
                v_star[j][i] = 0;
            }
            else if (i == 0) {      /* Just to the right of the left wall */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = 0;
                v_s = 0;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i] + v[j][i + 1]) / 2;

                Hy = (pow(v_cc, 2) + pow(v_cc_jm1, 2)) / D_y + (u_s_ip1 * v_s_ip1 - u_s * v_s) / D_x;

                v_star[j][i] = D_t * (Hy + (1 / Re) * (((1 / D_x) * (v[j][i + 1] - v[j][i]) + (2 / D_x) * (-v[j][i])) / D_x + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_y, 2)) * v[j][i]) + v[j][i];
            }
            else if (i == nodes_x - 2) {     /* Just the left of the right wall */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = 0;
                v_s_ip1 = 0;

                Hy = (pow(v_cc, 2) + pow(v_cc_jm1, 2)) / D_y + (u_s_ip1 * v_s_ip1 - u_s * v_s) / D_x;

                v_star[j][i] = D_t * (Hy + (1 / Re) * (((2 / D_x) * (-v[j][i]) + (1 / D_x) * (v[j][i - 1] - v[j][i])) / D_x + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_y, 2)) * v[j][i]) + v[j][i];
            }
            else if (i == nodes_x - 1) {      /* At the right wall */
                v_star[j][i] = 0;
            }
        }
    }

    /* ---------------------------------------------------------------------------------------------------
    Step 2 ------------------------------------------------------------------------------------------*/

    for (j = 0; j < nodes_y - 1; j++)
    { /* Compute f on interior points  */
        for (i = 0; i < nodes_x - 1; i++)
        {
            f[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / D_x + (v_star[j + 1][i] - v_star[j][i]) / D_y) / D_t;
        }
    }


    //for (j = 0; nodes_y - 1; j++)
    //{ /*compute f along  */
    //}

    //f_norm = 0; /* compute f_norm */
    //for (j = 0; j < nodes_y; j++)
    //{
    //    for (i = 0; i < nodes_x; i++)
    //    {
    //        f_norm =
    //    }
    //}

    //do
    //{
    //
    //} while ()

        /* ----------------------------------------------------------------------------------------------------------------------------
    Freeing memory ----------------------------------------------------------------------------------------------------------------
    */
    for (i = 0; i < nodes_x; i++)
    {
        free(u[i]);
        free(u_star[i]);
        free(p[i]);
        free(p_new[i]);
        free(f[i]);
    }
    for (i = 0; i < nodes_y; i++)
    {
        free(v[i]);
        free(v_star[i]);
    }
    free(u);
    free(v);
    free(u_star);
    free(v_star);
    free(p);
    free(p_new);
    free(f);
    printf("Done\n");

    return 0;
}