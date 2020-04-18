#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

int main() {
    /* ----------------------------------------------------------------------------------------------------------------------------
    Initializing Variables  -------------------------------------------------------------------------------------------------------
    */
    int i = 0; /* Counters for all "for loops" */
    int j = 0;

    double Re = 100;     /* Problem parameters */
    double D_t = 0.001;
    int nodes_x = 20;
    int nodes_y = 20;
    double NX = nodes_x;
    double NY = nodes_y;
    double D_x = 1 / NX;
    double D_y = 1 / NY;

    double** u = (double**)calloc(nodes_x, sizeof(double*));    /* Memory allocation for large arrays (velocities, etc.) */
    double** v = (double**)calloc(nodes_y, sizeof(double*));    /* u and v represent barycentric velocities */
    double** u_star = (double**)calloc(nodes_x, sizeof(double*));
    double** v_star = (double**)calloc(nodes_y, sizeof(double*));
    for (i = 0; i < nodes_x; i++) {
        u[i] = (double*)calloc(nodes_y, sizeof(double));
        u_star[i] = (double*)calloc(nodes_y, sizeof(double));
    }
    for (i = 0; i < nodes_y; i++) {
        v[i] = (double*)calloc(nodes_x, sizeof(double));
        v_star[i] = (double*)calloc(nodes_x, sizeof(double));
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
            if (i != 0 && j != 0 && i != (nodes_x - 1) && j != (nodes_y - 1)) {  /* NOTE: Need to check and correct all of this for staggering; boarders excluded */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j][i]) / 2;
                v_s_ip1 = (v[j][i + 1] + u[j][i]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2; /**/
                v_s_jp1 = (u[j + 1][i] + u[j + 1][i - 1]) / 2; /**/

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / D_x) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;
                Hy = ((pow(v[j][i], 2) - pow(v[j][i - 1], 2)) / D_y) + (v[j + 1][i] * u[j + 1][i] - v[j][i] * u[j][i]) / D_x;
                u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(D_x, 2) + (u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_y, 2))) + u[j][i];
                v_star[j][i] = D_t * (Hy + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(D_y, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_x, 2))) + v[j][i];
            }
        }
    }



    /* ----------------------------------------------------------------------------------------------------------------------------
    Freeing memory ----------------------------------------------------------------------------------------------------------------
    */
    for (i = 0; i < nodes_x; i++) {
        free(u[i]);
        free(u_star[i]);
    }
    for (i = 0; i < nodes_y; i++) {
        free(v[i]);
        free(v_star[i]);
    }
    free(u);
    free(v);
    free(u_star);
    free(v_star);
    printf("Done\n");

    return 0;
}