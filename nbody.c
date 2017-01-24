#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <plplot/plplot.h>
#include "nbody.h"

int main(){
    printf("\n=============================================================\n");
    printf("This program is able to simulate an n-body problem \n");
    printf("\n");
    printf("============================================================\n");

    //==========================================================================
    //-----------------------INITIALIZATIONS------------------------------------
    //==========================================================================
    srand(time(NULL));

    //==========================================================================
    //-----------------------VARIABLE INITIALIZATIONS---------------------------
    //==========================================================================

    // generic counters
    int i, j;

    // physical system
    environment env;
    double meanMotionResonanceRadius[7];

    // plotting variables
    float plotMass;
    char title[50];
    // color map matrix
    const double cmapi[2] = {0.0, 1.0};                                        // colour map intensity
    const double cmapr[2] = {0.7, 0.0};                                        // red
    const double cmapg[2] = {0.0, 0.0};                                        // green
    const double cmapb[2] = {0.0, 1.0};                                        // blue
    const int alt_hue_path = 1;

    //==========================================================================
    //-----------------------INITIALIZE VARIABLE VALUES-------------------------
    //==========================================================================
    env.nparticles = 1000;
    env.p = (particle*)malloc(sizeof(particle)*env.nparticles);
    env.t = 0.0;
    env.dt = 0.01;
    env.G = 1.0;
    env.bounded = 0;
    env.accretion = 0;

    // initialize particles
    for (i=0; i<env.nparticles; i++){
        initParticleRandom(&env.p[i]);
    }
    // get rid of particles on top of each other
    doAccretion(&env);
    for (i=0; i<env.nparticles; i++){
        particleSetVTangent(&env.p[i], env.G, 1.0, 0.0, 0.0, 0.0);
    }
    // initialize select particles
    // sun
    particleSetR(&env.p[0], 0.0, 0.0, 0.0);
    particleSetV(&env.p[0], 0.0, 0.0, 0.0);
    env.p[0].m = 1.0;
    // jupiter
    particleSetR(&env.p[1], 1.0, 0.0, 0.0);
    particleSetV(&env.p[1], 0.0, pow(env.G*env.p[0].m/env.p[1].r.x, 0.5), 0.0);
    env.p[1].m = 0.01;

    // calculate mean motion resonance for Jupiter
    meanMotionResonanceRadius[0] = env.p[1].r.x;
    meanMotionResonanceRadius[1] = pow(2.0*2.0/(env.p[1].r.x*env.p[1].r.x*env.p[1].r.x), -1.0/3.0);
    meanMotionResonanceRadius[2] = pow(3.0*3.0/(env.p[1].r.x*env.p[1].r.x*env.p[1].r.x), -1.0/3.0);
    meanMotionResonanceRadius[3] = pow((3.0/2.0)*(3.0/2.0)/(env.p[1].r.x*env.p[1].r.x*env.p[1].r.x), -1.0/3.0);
    meanMotionResonanceRadius[4] = pow(4.0*4.0/(env.p[1].r.x*env.p[1].r.x*env.p[1].r.x), -1.0/3.0);
    meanMotionResonanceRadius[5] = pow((4.0/3.0)*(4.0/3.0)/(env.p[1].r.x*env.p[1].r.x*env.p[1].r.x), -1.0/3.0);
    printf("%f, %f, %f, %f\n", meanMotionResonanceRadius[0], meanMotionResonanceRadius[1], meanMotionResonanceRadius[2], meanMotionResonanceRadius[3]);


    // asteroid
    particleSetR(&env.p[2], meanMotionResonanceRadius[4], 0.0, 0.0);
    particleSetV(&env.p[2], 0.0, pow(env.G*env.p[0].m/env.p[2].r.x, 0.5), 0.0);
    env.p[1].m = 0.01;

    // radii for histogram
    double *radii;
    radii = (double*)malloc(sizeof(double)*env.nparticles);



    //==========================================================================
    //-----------------------INITIALIZE PLOTTER---------------------------------
    //==========================================================================
    plsdev("xcairo");
    plsetopt("geometry", "1800x900");
    plinit();
    plscmap1l(1, 2, cmapi, cmapr, cmapg, cmapb, &alt_hue_path);
    //plssub(2,1);
    //pladv(1);
    //plvpor(0.0, 1.0, 0.0, 1.0);
    //plwind(-0.5, 1.5, -0.5, 1.5);
    plenv(-1.5, 1.5, -1.5, 1.5, 1, 0);
    pllab("x", "y", "Title");
    //plenv(0.0, 1.0, 0.0, 1.0, 1, 0);
    //pladv(2);
    //plwind(0.0, 1.0, 0.0, 1.0);
    //plenv(0.0, 1.0, 0.0, 1.0, 1, 0);



    //==========================================================================
    //-----------------------MAIN LOOP------------------------------------------
    //==========================================================================
    i=0;
    plclear();
    while(1) {
        env.t += env.dt;
        i++;

        orbitLeapfrogStep(&env);
        //calculate radii for histogram
        for (j=0; j<env.nparticles; j++){
            radii[j] = pow((env.p[j].r.x)*(env.p[j].r.x) + (env.p[j].r.y)*(env.p[j].r.y), 0.5);
        }

        // xy plot
        sprintf(title, "t=%f; i=%d", env.t, i);
        plvpor(0.1, 0.45, 0.2, 0.9);
        plwind(-1.5, 1.5, -1.5, 1.5);
        plclear();
        plcol0(1);
        pllab("x", "y", title);
        plbox("abcnt", 0.0, 0.0, "abcnt", 0.0, 0.0);
        for (j=0; j<env.nparticles; j++){
            if (env.p[j].m != 0.0) {
                plotMass = env.p[j].m / 0.01;
                if (plotMass > 1.0) plotMass = 1.0;
                plcol1(plotMass);
                pljoin(env.p[j].rold.x, env.p[j].rold.y, env.p[j].r.x, env.p[j].r.y);
                plpoin(1, &env.p[j].r.x, &env.p[j].r.y, 1);
            }
        }
        // histogram
        //pladv(2);
        //plenv(0.0, 1.0, 0.0, 100.0, 0, 1);
        plvpor(0.55, 0.9, 0.2, 0.9);
        plwind(0.0, 1.50, 0.0, 15.0);
        plcol0(3);
        plbox("abcnt", 0.0, 0.0, "abcnt", 0.0, 0.0);
        pllab("r", "nparticles", "Kirkwood gaps?");
        plhist(env.nparticles, radii, 0.0, 3.5, 300, PL_HIST_NOSCALING);
        plcol0(5);
        pljoin(meanMotionResonanceRadius[0], 0.0,meanMotionResonanceRadius[0], 40.0);
        plcol0(6);
        pljoin(meanMotionResonanceRadius[1], 0.0,meanMotionResonanceRadius[1], 40.0);
        plcol0(7);
        pljoin(meanMotionResonanceRadius[2], 0.0,meanMotionResonanceRadius[2], 40.0);
        pljoin(meanMotionResonanceRadius[3], 0.0,meanMotionResonanceRadius[3], 40.0);
        plcol0(8);
        pljoin(meanMotionResonanceRadius[4], 0.0,meanMotionResonanceRadius[4], 40.0);
        pljoin(meanMotionResonanceRadius[5], 0.0,meanMotionResonanceRadius[5], 40.0);



        // reset sun
        particleSetR(&env.p[0], 0.0, 0.0, 0.0);
        particleSetV(&env.p[0], 0.0, 0.0, 0.0);
        //env.p[0].m = 100.0;

        plflush();

    }

    plend();
}














float dvxdtGrav(float t, float x, float y){
    float r = pow(x*x + y*y, 0.5);
    return (-1.0 * x) / (r*r*r);
}
float dvydtGrav(float t, float x, float y){
    float r = pow(x*x + y*y, 0.5);
    return (-1.0 * y) / (r*r*r);
}

void orbitLeapfrogStep(environment *env){
    int i,j;
    // t and dt for easier access
    float t = env->t; float dt = env->dt;
    // vector distances
    float vecX, vecY;

    // accrete - do at the beginning to get rid of particles on top of each other and avoid division by 0 error;
    if (env->accretion != 0) doAccretion(env);

    // archive positions
    for (i=0; i<env->nparticles; i++){
        env->p[i].rold.x = env->p[i].r.x;
        env->p[i].rold.y = env->p[i].r.y;
        env->p[i].rold.z = env->p[i].r.z;
    }

    // kick velocity half a timestep
    for (i=0; i<env->nparticles; i++){
        //for (j=0; j<env->nparticles; j++){
        for (j=0; j<2; j++){
            if (i != j){
                // measure vector distances
                vecX = env->p[i].r.x - env->p[j].r.x;
                vecY = env->p[i].r.y - env->p[j].r.y;
                if (env->bounded != 0) {
                    while (vecX > 1.0) vecX -= 2.0;
                    while (vecX < -1.0) vecX += 2.0;
                    while (vecY > 1.0) vecY -= 2.0;
                    while (vecY < -1.0) vecY += 2.0;
                }

                // do the kicking, from all other particles
                env->p[i].v.x += dvxdtGrav(t, vecX, vecY) * env->p[j].m * env->G * dt/2.0;
                env->p[i].v.y += dvydtGrav(t, vecX, vecY) * env->p[j].m * env->G * dt/2.0;
            }
        }
    }

    // drift position forward by full timestep
    for (i=0; i<env->nparticles; i++){
        env->p[i].r.x += env->p[i].v.x * dt;
        env->p[i].r.y += env->p[i].v.y * dt;
        if (env->bounded != 0){
            while (env->p[i].r.x > 1.0) {
                env->p[i].r.x -= 2.0;
                env->p[i].rold.x = env->p[i].r.x;
            }
            while (env->p[i].r.x < -1.0) {
                env->p[i].r.x += 2.0;
                env->p[i].rold.x = env->p[i].r.x;
            }
            while (env->p[i].r.y > 1.0) {
                env->p[i].r.y -= 2.0;
                env->p[i].rold.y = env->p[i].r.y;
            }
            while (env->p[i].r.y < -1.0) {
                env->p[i].r.y += 2.0;
                env->p[i].rold.y = env->p[i].r.y;
            }
        }
    }

    // kick velocity the rest of the timestep
    for (i=0; i<env->nparticles; i++){
        //for (j=0; j<env->nparticles; j++){
        for (j=0; j<2; j++){
            if (i != j){
                // measure vector distances
                vecX = env->p[i].r.x - env->p[j].r.x;
                vecY = env->p[i].r.y - env->p[j].r.y;
                if (env->bounded != 0){
                    while (vecX > 1.0) vecX -= 2.0;
                    while (vecX < -1.0) vecX += 2.0;
                    while (vecY > 1.0) vecY -= 2.0;
                    while (vecY < -1.0) vecY += 2.0;
                }
                // do the kicking, from all other particles
                env->p[i].v.x += dvxdtGrav(t, vecX, vecY) * env->p[j].m * env->G * dt/2.0;
                env->p[i].v.y += dvydtGrav(t, vecX, vecY) * env->p[j].m * env->G * dt/2.0;
            }
        }
    }

    /*
    // dampen velocity
    for (i=0; i<env->nparticles; i++){
        if((env->p[i].v.x*env->p[i].v.x + env->p[i].v.y*env->p[i].v.y ) >= 10.0){
            env->p[i].v.x *= 0.2;
            env->p[i].v.y *= 0.2;
        }
    }
    */

}

void particleDeleteShift(environment *env, int n){
    int i;
    for (i=n; i<env->nparticles; i++){
        memcpy(&env->p[i], &env->p[i+1], sizeof(particle));
    }
    env->nparticles --;
}

void particleCombine(particle *p1, particle *p2){
    p1->v.x = (p1->m*p1->v.x + p2->m*p2->v.x)/(p1->m + p2->m);
    p1->v.y = (p1->m*p1->v.y + p2->m*p2->v.y)/(p1->m + p2->m);
    p1->v.z = (p1->m*p1->v.z + p2->m*p2->v.z)/(p1->m + p2->m);
    p1->m += p2->m;
    p2->m = 0.0;
}

void particleSetR(particle *p, double x, double y, double z){
    p->r.x = x; p->r.y = y; p->r.z = z;
}
void particleSetV(particle *p, double x, double y, double z){
    p->v.x = x; p->v.y = y; p->v.z = z;
}

void doAccretion(environment *env){
    int i,j;
    double vecX, vecY;

    for (i=0; i<env->nparticles; i++){
        for (j=0; j<env->nparticles; j++){
            if (i != j){
                vecX = env->p[i].r.x - env->p[j].r.x;
                vecY = env->p[i].r.y - env->p[j].r.y;
                if ((vecX*vecX + vecY*vecY) < 0.001){
                    particleCombine(&env->p[i], &env->p[j]);
                    particleDeleteShift(env, j);
                    printf("total: %d; acrete %d with %d; mass: %f\n", env->nparticles, i, j, env->p[i].m);
                }
            }
        }
    }
}

void initParticleZeroes(particle *p){
    p->r.x = 0.0;    p->r.y = 0.0;    p->r.z = 0.0;      // position
    p->v.x = 0.0;    p->v.y = 0.0;    p->v.z = 0.0;      // velocity
    p->a.x = 0.0;    p->a.y = 0.0;    p->a.z = 0.0;      // acceleration

    p->rold.x = 0.0; p->rold.y = 0.0; p->rold.z = 0.0;
    p->m = 1.0;
}

void initParticleRandom(particle *p){
    double r = rand()%1000/1000.0;
    double theta = rand()%6283/1000.0;
    p->r.x = r*cos(theta);                  p->r.y = r*sin(theta); p->r.z = 0.0;
    p->v.x = (rand()%1000)/10.0-50.0;    p->v.y = (rand()%1000)/10.0-50.0;    p->v.z = 0.0;      // velocity
    p->a.x = 0.0;    p->a.y = 0.0;    p->a.z = 0.0;      // acceleration

    p->rold.x = p->r.x; p->rold.y = p->r.y; p->rold.z = p->r.z;
    //p->m = (rand()%1000)/1000000.0;
    p->m = 0.001;
}


void particleSetVTangent(particle *p, double G, double M, double centerX, double centerY, double centerZ){
    double x = p->r.x - centerX;
    double y = p->r.y - centerY;
    double r = pow(x*x+y*y,0.5);
    p->v.x = -1.0*y/r;
    p->v.y = x/r;
    p->v.x *= pow(G*M/r, 0.5);
    p->v.y *= pow(G*M/r, 0.5);
    //p->v.x += (rand()%1000)/10.0 - 50.0;
    //p->v.y += (rand()%1000)/10.0 - 50.0;
}

void errorCase(const int errorCode){
    system("cat nagato");
    switch(errorCode){
        case ERR_INVALID_INPUT:
			printf("Error: invalid input\n");
			exit(-1);
		case ERR_MALLOC_FAIL:
			printf("Error: out of memory\n");
			exit(-1);
		case ERR_FILE_OPEN:
			printf("Error: file cannot be opened\n");
			exit(-1);
		case ERR_PGPLOT:
			printf("Error: cannot open pgplot window\n");
			exit(-1);
    }
}
