#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "epic.h"
#include "image.h"
#include "array_types.h"
#include "epic_aux.h"

#include "omp.h"

#include "io.h"


/* given a set of matches, return the set of points in the first image where a match exists */
static int_image matches_to_seeds(image_t** flow, const int n_thread){
    int a=0;
    int y,x;
    //printf("%d \n", flow[0]->height);
    //printf("%d \n", flow[0]->width);
    for (y = 0; y < flow[0]->height; y++){
        for (x = 0; x < flow[0]->width; x++) {
            //printf("%f \n", flow[0]->data[y*flow[0]->stride+x]);
            //printf("%f \n", flow[1]->data[y*flow[1]->stride+x]);
            if ((flow[0]->data[y*flow[0]->stride+x] != 0) && (flow[1]->data[y*flow[1]->stride+x] != 0)) {
                a++;
            }
        }
        //printf("%ds \n", a);
    }
    //printf("%ds \n", a);
    int_image res = empty_image(int, 2, a);
    int i=0;
    for (y = 0; y < flow[0]->height; y++){
        for (x = 0; x < flow[0]->width; x++) {
            if ((flow[0]->data[y*flow[0]->stride+x] != 0) && (flow[1]->data[y*flow[1]->stride+x] != 0)) {
                res.pixels[2*i] = x;
                res.pixels[2*i+1] = y;
                i++;
            }
        }
    }
    return res;
}

/* given a set of matches, return the set of vecotrs of the matches*/
static float_image matches_to_vects(image_t** flow, const int n_thread){
    int a=0;
    int y,x;
    for (y = 0; y < flow[0]->height; y++){
        for (x = 0; x < flow[0]->width; x++) {
            if ((flow[0]->data[y*flow[0]->stride+x] != 0) && (flow[1]->data[y*flow[1]->stride+x] != 0)) {
                a++;
            }
        }
    }
    float_image res = empty_image(float, 2, a);
    int i=0;
    for (y = 0; y < flow[0]->height; y++){
        for (x = 0; x < flow[0]->width; x++) {
            if ((flow[0]->data[y*flow[0]->stride+x] != 0) && (flow[1]->data[y*flow[1]->stride+x] != 0)) {
                res.pixels[2*i] = flow[0]->data[y*flow[0]->stride+x];
                res.pixels[2*i+1] = flow[1]->data[y*flow[1]->stride+x];
                i++;
            }
        }
    }
    return res;
}


/* set params to default value */
void epic_params_default(epic_params_t* params){
    strcpy(params->method, "LA");
    params->saliency_th = 0.045f;
    params->pref_nn = 25;
    params->pref_th = 5.0f;
    params->nn = 100;
    params->coef_kernel = 0.8f;
    params->euc = 0.001f;
    params->verbose = 0;
}

/* main function for edge-preserving interpolation of correspondences
    flowx                  x-component of the flow (output)
    flowy                  y-component of the flow (output)
    input_matches          input matches with each line a match and the first four columns containing x1 y1 x2 y2 
    im                     first image (in lab colorspace)
    edges                  edges cost (can be modified)
    params                 parameters
    n_thread               number of threads
*/
void epic(image_t *flowx, image_t *flowy, const color_image_t *im, image_t** match, float_image* edges, const epic_params_t* params, const int n_thread){

    // prepare variables
    printf("1 \n");
    int_image seeds = matches_to_seeds(match, n_thread);
    printf("2 \n");
    float_image vects = matches_to_vects(match, n_thread);
    printf("3 \n");
    printf("%d \n", vects.ty);
    const int nns = MIN(params->nn, vects.ty);
    if( nns < params->nn ) fprintf(stderr, "Warning: not enough matches for interpolating\n");
    if( params->verbose ) printf("Computing %d nearest neighbors for each match\n", nns);
    printf("4 \n"); 
    // compute nearest matches for each seed
    int_image nnf = empty_image( int, nns, vects.ty);
    float_image dis = empty_image( float, nns, vects.ty);
    int_image labels = empty_image( int, edges->tx, edges->ty);
    dist_trf_nnfield_subset( &nnf, &dis, &labels, &seeds, edges, NULL, &seeds, n_thread);  
    printf("5 \n");       
    // apply kernel to the distance
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(int i=0 ; i<dis.tx*dis.ty ; i++){
        dis.pixels[i] = expf(-params->coef_kernel*dis.pixels[i])+1e-08;
    }
    printf("6 \n");
    // interpolation
    if( params->verbose ) printf("Interpolation of matches using %s\n", params->method);
    float_image newvects = empty_image( float, 2, im->width*im->height);
    if( !strcmp( params->method, "LA") ){
        float_image seedsaffine = empty_image(float, 6, vects.ty);
        fit_localaffine(&seedsaffine, &nnf, &dis, &seeds, &vects);
        apply_localaffine(&newvects, &seedsaffine, &labels, n_thread);
        free(seedsaffine.pixels);
    } else if ( !strcmp( params->method, "NW") ){
        float_image seedsvects = empty_image(float, 2, vects.ty);
        fit_nadarayawatson(&seedsvects, &nnf, &dis, &vects, n_thread);
        apply_nadarayawatson(&newvects, &seedsvects, &labels, n_thread);        
        free(seedsvects.pixels);
    } else {
        fprintf(stderr, "method %s not recognized\n", params->method);
        exit(EXIT_FAILURE);
    }
    printf("7 \n");
    // copy result to the output
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(int i=0 ; i<im->height ; i++){
        for( int j=0 ; j<im->width ; j++){
            flowx->data[i*im->stride+j] = newvects.pixels[2*(i*im->width+j)];
            flowy->data[i*im->stride+j] = newvects.pixels[2*(i*im->width+j)+1];
        }
    }  
    
    // free memory
    free(seeds.pixels);
    free(vects.pixels);
    free(nnf.pixels);
    free(dis.pixels);
    free(labels.pixels);
    free(newvects.pixels);
    //free(match.pixels);
}

