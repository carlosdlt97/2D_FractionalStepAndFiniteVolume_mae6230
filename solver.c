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
    int nodes_x = 20;
    int nodes_y = 20;
    double NX = nodes_x;
    double NY = nodes_y;
    double D_x = 1 / NX;
    double D_y = 1 / NY;
    double lambda = pow(D_x, -2);
    double f_norm;

    double **u = (double **)calloc(nodes_x, sizeof(double *)); /* Memory allocation for large arrays (velocities, etc.) */
    double **v = (double **)calloc(nodes_y, sizeof(double *));
    double **u_star = (double **)calloc(nodes_x, sizeof(double *));
    double **v_star = (double **)calloc(nodes_y, sizeof(double *));
    for (i = 0; i < nodes_x; i++)
    {
        u[i] = (double *)calloc(nodes_y, sizeof(double));
        u_star[i] = (double *)calloc(nodes_y, sizeof(double));
    }
    for (i = 0; i < nodes_y; i++)
    {
        v[i] = (double *)calloc(nodes_x, sizeof(double));
        v_star[i] = (double *)calloc(nodes_x, sizeof(double));
    }

    double **p = (double **)calloc(nodes_x, sizeof(double *));
    double **p_new = (double **)calloc(nodes_x, sizeof(double *));
    double **f = (double **)calloc(nodes_x, sizeof(double *));
    for (i = 0; i < nodes_x; i++)
    {
        p[i] = (double *)calloc(nodes_y, sizeof(double));
        p_new[i] = (double *)calloc(nodes_y, sizeof(double));
        f[i] = (double *)calloc(nodes_y, sizeof(double));
    }

    /* ---------------------------------------------------------------------------------------------------
    Initializing variables  --------------------------------------------------------------------------- 
    */

    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY; j++)
        {
            u[j][i] = 0;
            v[j][i] = 0;
        }
    }

    for (i = 0; i < NX; i++)
    { /* top boundary of the box */
        j = NY - 1;
        u[j][i] = 1;
    }

    /* ----------------------------------------------------------------------------------------------------------------------------
    Step 1 -----------------------------------------------------------------------------------------------------------------------
    */
    for (j = 0; j < nodes_y; j++)
    { /* Interior points  */
        for (i = 0; i < nodes_x; i++)
        {
            double Hx;
            double Hy;
            if (i != 0 && j != 0 && i != (nodes_x - 1) && j != (nodes_y - 1))
            { /* NOTE: Need to check and correct all of this for staggering */
                //printf("%d - %d \n", i, j);
                Hx = ((pow(u[i][j], 2) - pow(u[i - 1][j], 2)) / D_x) + (u[i][j + 1] * v[i][j + 1] - u[i][j] * v[i][j]) / D_y;
                Hy = ((pow(v[i][j], 2) - pow(v[i][j - 1], 2)) / D_y) + (v[i + 1][j] * u[i + 1][j] - v[i][j] * u[i][j]) / D_x;
                u_star[i][j] = D_t * (Hx + (1 / Re) * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / pow(D_x, 2) + (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / pow(D_y, 2))) + u[i][j];
                v_star[i][j] = D_t * (Hy + (1 / Re) * ((v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / pow(D_y, 2) + (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / pow(D_x, 2))) + v[i][j];
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

    for (j = 0; nodes_y - 1; j++)
    { /*compute f along  */
    }

    f_norm = 0; /* compute f_norm */
    for (j = 0; j < nodes_y; j++)
    {
        for (i = 0; i < nodes_x; i++)
        {
            f_norm =
        }
    }

    do
    {

    } while ()

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