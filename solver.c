#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>



//function used to print the data to a text file
int print_current_data(int step, double** u, double** v, double** p, int NX, int NY) {
    int i, j;


    char filename1[25] = "step_";
    itoa(step, filename1 + 4, 10);
    strcat(filename1, "_u_data.txt");

    FILE* fpointer1 = fopen(filename1, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer1, "%.20lf ", u[j][i]);
        }
        fprintf(fpointer1, "\n");
    }
    fclose(fpointer1);



    char filename2[25] = "step_";
    itoa(step, filename2 + 4, 10);
    strcat(filename2, "_v_data.txt");

    FILE* fpointer2 = fopen(filename2, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer2, "%.20lf ", v[j][i]);
        }
        fprintf(fpointer2, "\n");
    }
    fclose(fpointer2);




    char filename3[25] = "step_";
    itoa(step, filename3 + 4, 10);
    strcat(filename3, "_p_data.txt");

    FILE* fpointer3 = fopen(filename3, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer3, "%.20lf ", p[j][i]);
        }
        fprintf(fpointer3, "\n");
    }
    fclose(fpointer3);

    return 1;
}













int main()
{
    /* ----------------------------------------------------------------------------------------------------------------------------
    Initializing Variables  -------------------------------------------------------------------------------------------------------
    */
    int i = 0; /* Counters for all "for loops" */
    int j = 0;
    int k = 0;
    int iter = 0;

    double Re = 100; /* Problem parameters */
    double D_t = 0.001;
    int nodes_x = 11;
    int nodes_y = 11;
    double NX = nodes_x;
    double NY = nodes_y;
    double D_x = 1 / (NX - 1);
    double D_y = 1 / (NY - 1);
    double lambda = pow(D_x, -2);
    double f_norm;
    double epsilon = pow(10, -5);
    double laplace_p_minus_f_norm;
    double RHS;
    double num_steps = 1000;



    /* Initializing varibales for this step 1 */
    double Hx;
    double Hy;
    double u_cc; /* cc = Cell cented velocity */
    double v_cc;
    double u_cc_im1; /* im1 = "i minus 1"; cell centered velocity at 1 position back in the i direction */
    double v_cc_jm1; /* jm1 = "j minus 1"; cell centered velocity at 1 position back in the j direction */
    double u_s;  /* s = staggered velocity (i-1/2, j-1/2) from the cc values */
    double v_s;
    double u_s_ip1; /* ip1 = "i plus 1"; staggered velocity at 1 position forward in the i direction */
    double u_s_jp1; /* jp1 = "j plus 1"; staggered velocity at 1 position forward in the j direction */
    double v_s_ip1;
    double v_s_jp1;




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
    double** laplace_p = (double**)calloc(nodes_x, sizeof(double*));
    for (i = 0; i < nodes_x; i++)
    {
        p[i] = (double*)calloc(nodes_y, sizeof(double));
        p_new[i] = (double*)calloc(nodes_y, sizeof(double));
        f[i] = (double*)calloc(nodes_y, sizeof(double));
        laplace_p[i] = (double*)calloc(nodes_y, sizeof(double));
    }


    for (iter = 0; iter < num_steps; iter++) {

        printf("\nstep %d\n", iter + 1);

        /* ----------------------------------------------------------------------------------------------------------------------------
        Step 1. -----------------------------------------------------------------------------------------------------------------------
        */

        for (j = 0; j < nodes_y; j++) {
            for (i = 0; i < nodes_x; i++) {
                /* Calculating u_star */
                if (i > 0 && j > 0 && i < (nodes_x - 1) && j < (nodes_y - 2)) {  /* At interior points except at the edges (CHECKED) */
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                    Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / D_x) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;

                    u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_x, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(D_y, 2))) + u[j][i];
                }
                else if (i == 0) {      /* At the left wall (CHECKED) */
                    u_star[j][i] = 0;
                }
                else if (i == nodes_x - 1) {     /* At the right wall (CHECKED) */
                    u_star[j][i] = 0;
                }
                else if (j == 0) {      /* Just above the bottom wall (CHECKED) */
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = 0;
                    v_s = 0;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i - 1] + v[j + 1][i]) / 2;

                    Hx = (pow(u_cc, 2) - pow(u_cc_im1, 2)) / D_x + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;

                    u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_x, 2) + ((1 / D_y) * (u[j + 1][i] - u[j][i]) + (2 / D_y) * (-u[j][i])) / D_y)) + u[j][i];
                }
                else if (j == nodes_y - 2) {    /* Just below the top lid (CHECKED) */
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i - 1] + v[j][i]) / 2;
                    u_s_jp1 = 1;
                    v_s_jp1 = 0;

                    Hx = (pow(u_cc, 2) - pow(u_cc_im1, 2)) / D_x + (u_s_jp1 * v_s_jp1 - u_s * v_s) / D_y;

                    u_star[j][i] = D_t * (Hx + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(D_x, 2) + ((2 / D_y) * (1 - u[j][i]) + (1 / D_y) * (u[j - 1][i] - u[j][i])) / D_y)) + u[j][i];
                }
                else if (j == nodes_y - 1) {     /* At the top lid (CHECKED) */
                    u_star[j][i] = 1;
                }
                /* Calculating v_star */
                if (i > 0 && j > 0 && i < (nodes_x - 2) && j < (nodes_y - 1)) {  /* At interior points except at the edges (CHECKED) */
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / D_y) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / D_x;

                    v_star[j][i] = D_t * (Hy + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(D_x, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_y, 2))) + v[j][i];
                }
                else if (j == 0) {    /* At the bottom wall (CHECKED) */
                    v_star[j][i] = 0;
                }
                else if (j == nodes_y - 1) {     /* At the top lid (CHECKED) */
                    v_star[j][i] = 0;
                }
                else if (i == 0) {      /* Just to the right of the left wall (CHECKED) */
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;
                    
                    u_s = 0;
                    v_s = 0;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i] + v[j][i + 1]) / 2;

                    Hy = (pow(v_cc, 2) - pow(v_cc_jm1, 2)) / D_y + (u_s_ip1 * v_s_ip1 - u_s * v_s) / D_x;

                    v_star[j][i] = D_t * (Hy + (1 / Re) * (((1 / D_x) * (v[j][i + 1] - v[j][i]) + (2 / D_x) * (-v[j][i])) / D_x + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_y, 2))) + v[j][i];
                }
                else if (i == nodes_x - 2) {     /* Just the left of the right wall (CHECKED) */
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = 0;
                    v_s_ip1 = 0;

                    Hy = (pow(v_cc, 2) - pow(v_cc_jm1, 2)) / D_y + (u_s_ip1 * v_s_ip1 - u_s * v_s) / D_x;

                    v_star[j][i] = D_t * (Hy + (1 / Re) * (((2 / D_x) * (-v[j][i]) + (1 / D_x) * (v[j][i - 1] - v[j][i])) / D_x + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(D_y, 2))) + v[j][i];
                }
                else if (i == nodes_x - 1) {      /* At the right wall (CHECKED) */
                    v_star[j][i] = 0;
                }
            }
        }

        printf("check1\n");

        /* ---------------------------------------------------------------------------------------------------
       Step 2 ------------------------------------------------------------------------------------------*/

        for (j = 0; j < nodes_y - 2; j++)
        { /* Compute f on interior points (and left boundary and bottom boundary in diagram)  */
            for (i = 0; i < nodes_x - 2; i++)
            {
                f[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / D_x + (v_star[j + 1][i] - v_star[j][i]) / D_y) / D_t;
            }
        }


        for (j = 0; j < nodes_y - 2; j++) {    /*compute f along right boundary */
            i = nodes_x - 2;
            f[j][i] = ((0 - u_star[j][i]) / D_x + (v_star[j + 1][i] - v_star[j][i]) / D_y) / D_t;
        }

        for (i = 0; i < nodes_x - 2; i++) {    /*compute f along top boundary */
            j = nodes_y - 2;
            f[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / D_x + (0 - v_star[j][i]) / D_y) / D_t;
        }

        f[nodes_y - 2][nodes_x - 2] = ((0 - u_star[nodes_y - 2][nodes_x - 2]) / D_x + (0 - v_star[nodes_y - 2][nodes_x - 2]) / D_y) / D_t;

         /* PRINTING TO CHECK VALUES (DELETE THIS LATER)
        printf("f\n");
        for (j = nodes_y - 1; j > -1; j--) {
            for (i = 0; i < nodes_x; i++) {
                printf("%f,  ", f[j][i]);
            }
            printf("\n");
        }
        char ch;
        scanf("%c", &ch);
        */

        f_norm = 0; /* compute f_norm and set p = 0 and laplace_p everywhere */
        for (j = 0; j < nodes_y - 1; j++) {
            for (i = 0; i < nodes_x - 1; i++) {
                f_norm = f_norm + pow(f[j][i], 2);
            }
        }

        f_norm = sqrt(f_norm);

        /* PRINTING TO CHECK VALUES (DELETE THIS LATER)
        printf("f_norm = %f\n", f_norm);
        */


        do {


            //update interior values
            for (j = 1; j < nodes_y - 2; j++) {
                for (i = 1; i < nodes_x - 2; i++) {

                    p_new[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                }

            }

            //update left boundary values
            for (j = 1; j < nodes_y - 2; j++) {
                i = 0;
                p_new[j][i] = (p[j][i + 1] + p[j - 1][i] + p[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

            }

            //update right boundary values
            for (j = 1; j < nodes_y - 2; j++) {
                i = nodes_y - 2;
                p_new[j][i] = (p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

            }

            //update bottom boundary values 
            for (i = 1; i < nodes_x - 2; i++) {
                j = 0;
                p_new[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
            }

            //update top boundary values
            for (i = 1; i < nodes_x - 2; i++) {
                j = nodes_x - 2;
                p_new[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
            }

            //update corner points
            p_new[0][0] = (p[1][0] + p[0][1]) / 2 - f[0][0] / (2 * lambda);
            p_new[0][nodes_x - 2] = (p[1][nodes_x - 2] + p[0][nodes_x - 3]) / 2 - f[0][nodes_x - 2] / (2 * lambda);
            p_new[nodes_y - 2][0] = (p[nodes_y - 2][1] + p[nodes_y - 3][0]) / 2 - f[nodes_y - 2][0] / (2 * lambda);
            p_new[nodes_y - 2][nodes_x - 2] = (p[nodes_y - 2][nodes_x - 3] + p[nodes_y - 3][nodes_x - 2]) / 2 - f[nodes_y - 2][nodes_x - 2] / (2 * lambda);

            //update p matrix with p_new values
            for (j = 0; j < nodes_y - 1; j++) {
                for (i = 0; i < nodes_x - 1; i++) {
                    p[j][i] = p_new[j][i];
                }
            }
            
            /* PRINTING TO CHECK VALUES (DELETE THIS LATER)
            printf("p\n");
            for (j = nodes_y - 1; j > -1; j--) {
                for (i = 0; i < nodes_x; i++) {
                    printf("%f,  ", p[j][i]);
                }
                printf("\n");
            }
            char ch;
            scanf("%c", &ch);
            */

            //compute laplace_p matrix
            for (j = 1; j < nodes_y - 2; j++) {
                for (i = 1; i < nodes_x - 2; i++) {
                    laplace_p[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) * lambda - 4 * lambda * p[j][i];
                }
            }
            for (j = 1; j < nodes_y - 2; j++) {
                i = 0;
                laplace_p[j][i] = (p[j][i + 1] + p[j - 1][i] + p[j + 1][i]) * lambda - p[j][i] * (3 * lambda);
            }
            for (j = 1; j < nodes_y - 2; j++) {
                i = nodes_x - 2;
                laplace_p[j][i] = (p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) * lambda - p[j][i] * (3 * lambda);
            }
            for (i = 1; i < nodes_x - 2; i++) {
                j = 0;
                laplace_p[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j + 1][i]) * lambda - p[j][i] * (3 * lambda);
            }
            for (i = 1; i < nodes_x - 2; i++) {
                j = nodes_y - 2;
                laplace_p[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i]) * lambda - p[j][i] * (3 * lambda);
            }
            laplace_p[0][0] = (p[1][0] + p[0][1]) * lambda - p[0][0] * (2 * lambda);
            laplace_p[0][nodes_x - 2] = (p[1][nodes_x - 2] + p[0][nodes_x - 3]) * lambda - p[0][nodes_x - 2] * (2 * lambda);
            laplace_p[nodes_y - 2][0] = (p[nodes_y - 2][1] + p[nodes_y - 3][0]) * lambda - p[nodes_y - 2][0] * (2 * lambda);
            laplace_p[nodes_y - 2][nodes_x - 2] = (p[nodes_y - 2][nodes_x - 3] + p[nodes_y - 3][nodes_x - 2]) * lambda - p[nodes_y - 2][nodes_x - 2] * (2 * lambda);

            /* PRINTING TO CHECK VALUES (DELETE THIS LATER)
            printf("laplace_p\n");
            for (j = nodes_y - 1; j > -1; j--) {
                for (i = 0; i < nodes_x; i++) {
                    printf("%f,  ", laplace_p[j][i]);
                }
                printf("\n");
            }
            char ch;
            scanf("%c", &ch);
            */

            //compute the norm
            laplace_p_minus_f_norm = 0;

            for (j = 0; j < nodes_y - 1; j++) {
                for (i = 0; i < nodes_x - 1; i++) {
                    laplace_p_minus_f_norm = laplace_p_minus_f_norm + pow((laplace_p[j][i] - f[j][i]), 2);
                }
            }

            if (iter + 1 == -1) {
                /* PRINTING TO CHECK VALUES (DELETE THIS LATER)
                printf("f_norm = %f\n", f_norm);
                */

                /* PRINTING TO CHECK VALUES (DELETE THIS LATER) */
                printf("laplace_p_minus_f_norm at each point\n");
                for (j = nodes_y - 1; j > -1; j--) {
                    for (i = 0; i < nodes_x; i++) {
                        printf("%f,  ", laplace_p[j][i] - f[j][i]);
                    }
                    printf("\n");
                }
                /*char ch;
                scanf("%c", &ch);*/
            }

            laplace_p_minus_f_norm = sqrt(laplace_p_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


        } while (laplace_p_minus_f_norm > RHS);

        printf("check2\n");

        /* ----------------------------------------------------------------------------------------------------------------------------
        Step 3 -----------------------------------------------------------------------------------------------------------------------
        */



        //update interior points
        for (j = 1; j < nodes_y - 1; j++) {
            for (i = 1; i < nodes_x - 1; i++) {

                u[j][i] = u_star[j][i] - D_t * (p[j][i] - p[j][i - 1]);

                v[j][i] = v_star[j][i] - D_t * (p[j][i] - p[j - 1][i]);
            }
        }


        //left boundary
        for (j = 1; j < nodes_y - 1; j++) {
            i = 0;
            u[j][i] = 0;
            v[j][i] = v_star[j][i] - D_t * (p[j][i] - p[j - 1][i]);
        }

        //right boundary 
        for (j = 1; j < nodes_y - 1; j++) {
            i = nodes_x - 1;
            u[j][i] = 0;
            v[j][i] = 0;
        }

        //bottom boundary
        for (i = 1; i < nodes_x - 1; i++) {
            j = 0;
            u[j][i] = u_star[j][i] - D_t * (p[j][i] - p[j][i - 1]);
            v[j][i] = 0;
        }

        //top boundary 
        for (i = 1; i < nodes_x - 1; i++) {
            j = nodes_y - 1;
            u[j][i] = 1;
            v[j][i] = 0;
        }


        //corner points
        u[0][0] = 0;
        v[0][0] = 0;

        u[0][nodes_x - 1] = 0;
        v[0][nodes_x - 1] = 0;

        u[nodes_y - 1][0] = 0;
        v[nodes_y - 1][0] = 0;

        u[nodes_y - 1][nodes_x - 1] = 0;
        v[nodes_y - 1][nodes_x - 1] = 0;


    printf("check3\n");
    }


    /* ---------------------------------------------------------------------------------------------------
     Print to a text file ---------------------------------------------------------------------------------------
     */

    int something;
    something = print_current_data(num_steps, u, v, p, nodes_x, nodes_y);



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
        free(laplace_p[i]);
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
    free(laplace_p);
    printf("\n Done\n");

    
    return 0;
}