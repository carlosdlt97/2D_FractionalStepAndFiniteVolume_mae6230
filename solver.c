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


/* ==========================================================================
======================= n-Step Gauss-Seidel Solver ==========================
========================================================================== */

int GS_nstep(double** f, double** phi, int Nx, int Ny, double epsilon, int nGS) {

    /* Initilaizations */
    int i, j;
    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;

    double  f_norm;
    double integral;
    double lambda = pow(D_x, -2);
    double RHS;    
    double laplace_phi_minus_f_norm;
    
    int step = 1;
    int max_num_steps = nGS;

    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    for (j = 0; j < Ny; j++) {
        laplace_phi[j] = (double*)calloc(Nx, sizeof(double));
    }


    /* Solving -------------------------------------------------------------------------------- */

    f_norm = 0; /* compute f_norm */
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            f_norm = f_norm + pow(f[j][i], 2);
        }
    }

    f_norm = sqrt(f_norm);


    if (Nx > 1 && Ny > 1) {

        do {

            if (Nx > 2 && Ny > 2) {

                for (j = 1; j < Ny - 1; j++) {
                    for (i = 1; i < Nx - 1; i++) {
                        phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                    }
                }
                

                //update left boundary values
                for (j = 1; j < Ny - 1; j++) {
                    i = 0;
                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update right boundary values
                for (j = 1; j < Ny - 1; j++) {
                    i = Nx - 1;
                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update bottom boundary values 
                for (i = 1; i < Nx - 1; i++) {
                    j = 0;
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

                //update top boundary values
                for (i = 1; i < Nx - 1; i++) {
                    j = Ny - 1;
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

            }


            //update corner points
            phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
            phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
            phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
            phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

            
            if (Nx > 2 && Ny > 2) {

                //compute laplace_p matrix
                for (j = 1; j < Ny - 1; j++) {
                    for (i = 1; i < Nx - 1; i++) {
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                    }
                }
                for (j = 1; j < Ny - 1; j++) {
                    i = 0;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (j = 1; j < Ny - 1; j++) {
                    i = Nx - 1;
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (i = 1; i < Nx - 1; i++) {
                    j = 0;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (i = 1; i < Nx - 1; i++) {
                    j = Ny - 1;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }

            }

            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
                }
            }

            

            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


            /* break out of the loop if the max number of steps has been reached */
            if (step == max_num_steps) {
                break;
            }


           step += 1;

        } while (laplace_phi_minus_f_norm > RHS);
    }


    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            integral += phi[j][i] * D_x * D_y;
        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            phi[j][i] = phi[j][i] - integral / (Nx * Ny);
        }
    }


    for (i = 0; i < nx; i++) {
        free(laplace_phi[i]);
    }

    free(laplace_phi);

    return 0;

}

/* ====================================================================================================================
=================================== MG_recursion ======================================================================
==================================================================================================================== */


int MG_recursion(double** f, double** phi, int Nx, int Ny, double epsilon, int nGS) {
    
    /* Initializing variables */
    int finished;
    int i, j;

    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;
    double lambda = pow(D_x, -2);

    double half_nx = ceil(nx / 2);
    double half_ny = ceil(ny / 2);
    double half_D_x = 1 / half_nx;
    double half_D_y = 1 / half_nx;
    double half_lambda = pow(half_D_x, -2);
    int half_Nx = half_nx;
    int half_Ny = half_ny;

    double** residual = (double**)calloc(Ny, sizeof(double*));
    double** residual2 = (double**)calloc(half_Ny, sizeof(double*));
    double** error = (double**)calloc(Ny, sizeof(double*));
    double** error2 = (double**)calloc(half_Ny, sizeof(double*));
    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    
    for (i = 0; i < Nx; i++) {
        residual[i] = (double*)calloc(Nx, sizeof(double));
        error[i] = (double*)calloc(Nx, sizeof(double));
        laplace_phi[i] = (double*)calloc(Nx, sizeof(double));
    }
    for (i = 0; i < half_Nx; i++) {
        residual2[i] = (double*)calloc(half_Nx, sizeof(double));
        error2[i] = (double*)calloc(half_Nx, sizeof(double));
    }
    
    /* ================================================================================================*/
    
    /* ------- GETTING RESIDUAL ------- */
    /* compute laplace_p matrix */
    if (Nx > 1 && Ny > 1) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
            }
        }
        for (j = 1; j < Ny - 1; j++) {
            i = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (j = 1; j < Ny - 1; j++) {
            i = Nx - 1;
            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = Ny - 1;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }

        laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
        laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
        laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
        laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                residual[j][i] = f[j][i] - laplace_phi[j][i];

            }
        }
    }

    

    /* ------- RESTRICTION ------- */
    /* Restricting residual by taking average of a point and its surounding 4 points and placing it into residual2 (adjusted for edges & corners) */
    if (Nx > 2 && Ny > 2) { /* For residual sizes greater than 2x2 */

        /* For interior points */
        for (j = 1; j < (half_Ny - 1); j++){
            for(i = 1; i < (half_Nx - 1);  i++) {

                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 - 1] + residual[j * 2 - 1][i * 2]) / 5;
            
            }
        }

        /* For the left boundary */
        for (j = 1; j < (half_Ny - 1); j++) {
            i = 0;
            residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2 - 1][i * 2]) / 4;

        }

        /* For the bottom boundary */
        for (i = 1; j < (half_Nx - 1); j++) {
            j = 0;
            residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2][i * 2 - 1]) / 4;

        }

        /* ODD */
        if (nx / half_nx != 2 && ny / half_ny != 2) { /* For odd-sized square meshes (7x7, 23x23, etc.) */
            /* For the right boundary */
            for (j = 1; j < (half_Ny - 1); j++) {
                i = half_Nx - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2 - 1][i * 2]) / 4;

            }

            /* For the top boundary */
            for (i = 1; j < (half_Nx - 1); j++) {
                j = half_Ny - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2][i * 2 - 1]) / 4;

            }

            /* For the corners */
            residual2[0][0] = (residual[0][0] + residual[0][1] + residual[1][0]) / 3; /* bottom left */
            residual2[0][half_Nx - 1] = (residual[0][(half_Nx - 1) * 2] + residual[0][((half_Nx - 1) * 2) - 1] + residual[1][(half_Nx - 1) * 2]) / 3; /* bottom right */
            residual2[half_Ny - 1][0] = (residual[(half_Ny - 1) * 2][0] + residual[((half_Ny - 1) * 2) - 1][0] + residual[(half_Ny - 1) * 2][1]) / 3; /* top left */
            residual2[half_Ny - 1][half_Nx - 1] = (residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2] +  + residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2 - 1] + residual[(half_Ny - 1) * 2 - 1][(half_Nx - 1) * 2]) / 3; /* top right */
            
        }
        /* EVEN */
        else if (nx / half_nx == 2 && ny / half_ny == 2){ /* For even-sized square meshes (6x6, 24x24, etc.) */
                /* For the right boundary */
            for (j = 1; j < (half_Ny - 1); j++) {
                i = half_Nx - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 - 1] + residual[j * 2 - 1][i * 2]) / 5;

            }

            /* For the top boundary */
            for (i = 1; j < (half_Nx - 1); j++) {
                j = half_Ny - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 - 1] + residual[j * 2 - 1][i * 2]) / 5;
            
            }

            /* For the corners */
            residual2[0][0] = (residual[0][0] + residual[0][1] + residual[1][0]) / 3; /* bottom left */
            residual2[0][half_Nx - 1] = (residual[0][(half_Nx - 1) * 2] + residual[0][((half_Nx - 1) * 2) + 1] + residual[0][((half_Nx - 1) * 2) - 1] + residual[1][(half_Nx - 1) * 2]) / 4; /* bottom right */
            residual2[half_Ny - 1][0] = (residual[(half_Ny - 1) * 2][0] + residual[((half_Ny - 1) * 2) + 1][0] + residual[((half_Ny - 1) * 2) - 1][0] + residual[(half_Ny - 1) * 2][1]) / 4; /* top left */
            residual2[half_Ny - 1][half_Nx - 1] = (residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2] + residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2 + 1] + residual[(half_Ny - 1) * 2 + 1][(half_Nx - 1) * 2] + residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2 - 1] + residual[(half_Ny - 1) * 2 - 1][(half_Nx - 1) * 2]) / 5; /* top right */

        }
        else { /* Making sure the mesh is square, and terminating recursion otherwise (also frees memory) */
            printf("\n\nYOU DID NOT INPUT A SQUARE MESH.\n\n");

            for (j = 0; j < Ny; j++) {
                free(residual[j]);
                free(error[j]);
                free(laplace_phi[j]);
            }
            for (j = 0; j < half_Ny; j++) {
                free(residual2[j]);
                free(error2[j]);
            }
            free(residual);
            free(residual2);
            free(error);
            free(error2);
            free(laplace_phi);                

            return 0;
        }

    }
    else if (Nx == 2 && Ny == 2){ /* For phi & f sizes of 2x2 */
        residual2[0][0] = (residual[0][0] + residual[0][1] + residual[1][0]) / 3;
    }
    /* residual size of 1x1 is not restricted */


    /* ------- RELAXATION ON error ------- */
    GS_nstep(residual2, error2, half_Nx, half_Ny, epsilon, nGS);


    
    /* ------- RECURSION ------- */
    if (half_Nx > 1 && half_Ny > 1) {
        MG_recursion(residual2, error2, half_Nx, half_Ny, epsilon, nGS);
    }

    
    /* ------- PROLONGATION OF ERROR ------- */
    /* EVEN */
    if (nx / half_nx == 2 && ny / half_ny == 2){ /* For even-sized square meshes (6x6, 24x24, etc.) */
        for (j = 0; j < half_Ny; j++){
            for(i = 0; i < half_Nx; i++) {
                error[j * 2][i * 2] = error2[j][i];
                error[j * 2 + 1][i * 2] = error2[j][i];
                error[j * 2][i * 2 + 1] = error2[j][i];
                error[j * 2 + 1][i * 2 + 1] = error2[j][i];
            }
        }
    }
    /* ODD */
    else { /* For odd-sized square meshes (7x7, 23x23, etc.) */
        /* interior points */
        for (j = 0; j < half_Ny - 1; j++){
            for(i = 0; i < half_Nx - 1; i++) {
                error[j * 2][i * 2] = error2[j][i];
                error[j * 2 + 1][i * 2] = error2[j][i];
                error[j * 2][i * 2 + 1] = error2[j][i];
                error[j * 2 + 1][i * 2 + 1] = error2[j][i];
            }
        }
        
        /* edges */
        /* right edge */
        for (j = 0; j < half_Ny - 1; j++) {
            i = half_Nx - 1;
            error[j * 2][i * 2] = error2[j][i];
            error[j * 2 + 1][i * 2] = error2[j][i];
        }

        /* top edge */
        for (i = 0; i < half_Nx - 1; i++) {
            j = half_Ny - 1;
            error[j * 2][i * 2] = error2[j][i];
            error[j * 2][i * 2 + 1] = error2[j][i];
        }

        /* top-right corner */
        error[Ny - 1][Nx - 1] = error2[half_Ny - 1][half_Nx - 1];
    }

    
    /* ------- CORRECTION ------- */
    for (j = 0; j < Ny; j++){
        for(i = 0; i < Nx; i++) {

            phi[j][i] += error[j][i];

        }
    }

    /* Exit GS solving */
    finished = GS_nstep(f, phi, Nx, Ny, epsilon, nGS);


    for (j = 0; j < Ny; j++) {
        free(residual[j]);
        free(error[j]);
        free(laplace_phi[j]);
    }
    for (j = 0; j < half_Ny; j++) {
        free(residual2[j]);
        free(error2[j]);
    }
    free(residual);
    free(residual2);
    free(error);
    free(error2);
    free(laplace_phi);

    
    return 0;
}

/* ====================================================================================================================
=================================== Multigrid solver ==================================================================
==================================================================================================================== */


int Multigrid_solver(double** f, double** phi, int Nx, int Ny, double lambda, double epsilon, double epsilon2, int nGS, double f_norm) {
    
    int i, j;
    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;
    int step = 1;
    int max_num_steps = 2000;
    double laplace_phi_minus_f_norm;
    double RHS;
    double integral;

    
    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));

    for (i = 0; i < Nx; i++) {
        laplace_phi[i] = (double*)calloc(Nx, sizeof(double));
    }


    /* ------- RELAXATION ON phi ------- */
    GS_nstep(f, phi, Nx, Ny, epsilon, nGS);

    do {
        

        /* ------- RECURSION ------- */
        MG_recursion(f, phi, Nx, Ny, epsilon2, nGS);
        
        /* compute laplace_p matrix */
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
            }
        }
        for (j = 1; j < Ny - 1; j++) {
            i = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (j = 1; j < Ny - 1; j++) {
            i = Nx - 1;
            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = Ny - 1;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }

        laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
        laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
        laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
        laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);


        //compute the norm
        laplace_phi_minus_f_norm = 0;

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
            }
        }

        
        laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


        if (f_norm == 0) {
            RHS = epsilon;
        }
        else {
            RHS = epsilon * f_norm;
        }


        if (step == max_num_steps) {
            break;
        }


        step += 1;

    } while (laplace_phi_minus_f_norm > RHS);


    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            integral += phi[j][i] * D_x * D_y;
        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            phi[j][i] = phi[j][i] - integral / (Nx * Ny);
        }
    }

    /* Freeing memory */
    for (i = 0; i < Nx; i++) {
        free(laplace_phi[i]);
    }
    free(laplace_phi);

    
}




/* ============================================================================================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
=================================================  MAIN  ======================================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
============================================================================================================ */

int main()
{

    clock_t t_start = clock();

    /* ----------------------------------------------------------------------------------------------------------------------------
    Initializing Variables  -------------------------------------------------------------------------------------------------------
    */
    int i = 0; /* Counters for all "for loops" */
    int j = 0;
    int k = 0;
    int iter = 0;

    char solver_choice[] = "GS"; /* Must be either "MG" for multigrid method or "GS" for Gauss-Seidel method */

    double Re = 100; /* Problem parameters */
    double D_t = 0.001;
    int nodes_x = 129;
    int nodes_y = 129;
    double NX = nodes_x;
    double NY = nodes_y;
    double D_x = 1 / (NX - 1);
    double D_y = 1 / (NY - 1);
    double lambda = pow(D_x, -2);
    double f_norm;
    double epsilon = pow(10, -5);
    double epsilon2 = epsilon * pow(10, -2);
    double laplace_p_minus_f_norm;
    double RHS;
    double num_steps = 10000;
    int nGS = 7;



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


    double** p = (double**)calloc(nodes_y, sizeof(double*));
    double** p_new = (double**)calloc(nodes_y, sizeof(double*));
    double** f = (double**)calloc(nodes_y, sizeof(double*));
    double** laplace_p = (double**)calloc(nodes_y, sizeof(double*));
    
    for (i = 0; i < nodes_x; i++)
    {
        p[i] = (double*)calloc(nodes_x, sizeof(double));
        p_new[i] = (double*)calloc(nodes_x, sizeof(double));
        f[i] = (double*)calloc(nodes_x, sizeof(double));
        laplace_p[i] = (double*)calloc(nodes_x, sizeof(double));
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


        /* ---------------------------------------------------------------------------------------------------
       Step 2 ------------------------------------------------------------------------------------------*/

       /* COMPUTING f */

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

         
        f_norm = 0; /* compute f_norm and set p = 0 and laplace_p everywhere */
        for (j = 0; j < nodes_y - 1; j++) {
            for (i = 0; i < nodes_x - 1; i++) {
                f_norm = f_norm + pow(f[j][i], 2);
            }
        }

        f_norm = sqrt(f_norm);



        /* ------- Gauss-Seidel Solver ------- */

        if (strcmp(solver_choice, "GS") == 0) {

            do {

                //update interior values
                for (j = 1; j < nodes_y - 2; j++) {
                    for (i = 1; i < nodes_x - 2; i++) {

                        p[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                    }

                }

                //update left boundary values
                for (j = 1; j < nodes_y - 2; j++) {
                    i = 0;
                    p[j][i] = (p[j][i + 1] + p[j - 1][i] + p[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update right boundary values
                for (j = 1; j < nodes_y - 2; j++) {
                    i = nodes_y - 2;
                    p[j][i] = (p[j][i - 1] + p[j - 1][i] + p[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update bottom boundary values 
                for (i = 1; i < nodes_x - 2; i++) {
                    j = 0;
                    p[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

                //update top boundary values
                for (i = 1; i < nodes_x - 2; i++) {
                    j = nodes_x - 2;
                    p[j][i] = (p[j][i + 1] + p[j][i - 1] + p[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

                //update corner points
                p[0][0] = (p[1][0] + p[0][1]) / 2 - f[0][0] / (2 * lambda);
                p[0][nodes_x - 2] = (p[1][nodes_x - 2] + p[0][nodes_x - 3]) / 2 - f[0][nodes_x - 2] / (2 * lambda);
                p[nodes_y - 2][0] = (p[nodes_y - 2][1] + p[nodes_y - 3][0]) / 2 - f[nodes_y - 2][0] / (2 * lambda);
                p[nodes_y - 2][nodes_x - 2] = (p[nodes_y - 2][nodes_x - 3] + p[nodes_y - 3][nodes_x - 2]) / 2 - f[nodes_y - 2][nodes_x - 2] / (2 * lambda);

                

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

                

                //compute the norm
                laplace_p_minus_f_norm = 0;

                for (j = 0; j < nodes_y - 1; j++) {
                    for (i = 0; i < nodes_x - 1; i++) {
                        laplace_p_minus_f_norm = laplace_p_minus_f_norm + pow((laplace_p[j][i] - f[j][i]), 2);
                    }
                }

                
                laplace_p_minus_f_norm = sqrt(laplace_p_minus_f_norm);


                if (f_norm == 0) {
                    RHS = epsilon;
                }
                else {
                    RHS = epsilon * f_norm;
                }


            } while (laplace_p_minus_f_norm > RHS);

        }

        /* ------- Multigrid Solver ------- */
        /*   See Multigrid_solver function  */

        if (strcmp(solver_choice, "MG") == 0) {

            Multigrid_solver(f, p, nodes_x - 1, nodes_y - 1, lambda, epsilon, epsilon2, nGS, f_norm);

        }


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

    clock_t t_end = clock();
    double run_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    printf("run time = %f\n", run_time);

    return 0;
}