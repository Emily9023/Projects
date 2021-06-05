#include "seamcarving.h"
#include <stdio.h>
#include <math.h>

double find_min_dbl(double a, double b, double c){
    double min = a;
    if (min > b) {
        min = b;
    }
    if (min > c) {
        min = c;
    }
    return min;
}

int find_min_int(int a, int b, int c){
    int min = a;
    if (min > b) {
        min = b;
    }
    if (min > c) {
        min = c;
    }
    return min;
}

void calc_energy(struct rgb_img *im, struct rgb_img **grad) {

    create_img(grad, im -> height, im -> width);

    int y_max = im->height - 1;
    int x_max = im->width - 1;

    int delta_x, delta_y;
    int R_diff_x, R_diff_y;
    int G_diff_x, G_diff_y;
    int B_diff_x, B_diff_y;
    uint8_t energy;

    //corners

    // (0,0)
    R_diff_x = get_pixel(im, 0, 1, 0) - get_pixel(im, 0, x_max, 0);
    G_diff_x = get_pixel(im, 0, 1, 1) - get_pixel(im, 0, x_max, 1);
    B_diff_x = get_pixel(im, 0, 1, 2) - get_pixel(im, 0, x_max, 2);

    R_diff_y = get_pixel(im, 1, 0, 0) - get_pixel(im, y_max, 0, 0);
    G_diff_y = get_pixel(im, 1, 0, 1) - get_pixel(im, y_max, 0, 1);
    B_diff_y = get_pixel(im, 1, 0, 2) - get_pixel(im, y_max, 0, 2);

    delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
    delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

    energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
    set_pixel(*grad, 0, 0, energy, energy, energy);

    // (0, x_max)
    R_diff_x = get_pixel(im, 0, 0, 0) - get_pixel(im, 0, x_max - 1, 0);
    G_diff_x = get_pixel(im, 0, 0, 1) - get_pixel(im, 0, x_max - 1, 1);
    B_diff_x = get_pixel(im, 0, 0, 2) - get_pixel(im, 0, x_max - 1, 2);

    R_diff_y = get_pixel(im, 1, x_max, 0) - get_pixel(im, y_max, x_max, 0);
    G_diff_y = get_pixel(im, 1, x_max, 1) - get_pixel(im, y_max, x_max, 1);
    B_diff_y = get_pixel(im, 1, x_max, 2) - get_pixel(im, y_max, x_max, 2);

    delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
    delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

    energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
    set_pixel(*grad, 0, x_max, energy, energy, energy);

    // (y_max, 0)
    R_diff_x = get_pixel(im, y_max, 1, 0) - get_pixel(im, y_max, x_max, 0);
    G_diff_x = get_pixel(im, y_max, 1, 1) - get_pixel(im, y_max, x_max, 1);
    B_diff_x = get_pixel(im, y_max, 1, 2) - get_pixel(im, y_max, x_max, 2);

    R_diff_y = get_pixel(im, 0, 0, 0) - get_pixel(im, y_max - 1, 0, 0);
    G_diff_y = get_pixel(im, 0, 0, 1) - get_pixel(im, y_max - 1, 0, 1);
    B_diff_y = get_pixel(im, 0, 0, 2) - get_pixel(im, y_max - 1, 0, 2);

    delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
    delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

    energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
    set_pixel(*grad, y_max, 0, energy, energy, energy);

    // (y_max, x_max)
    R_diff_x = get_pixel(im, y_max, 0, 0) - get_pixel(im, y_max, x_max-1, 0);
    G_diff_x = get_pixel(im, y_max, 0, 1) - get_pixel(im, y_max, x_max-1, 1);
    B_diff_x = get_pixel(im, y_max, 0, 2) - get_pixel(im, y_max, x_max-1, 2);

    R_diff_y = get_pixel(im, 0, x_max, 0) - get_pixel(im, y_max - 1, x_max, 0);
    G_diff_y = get_pixel(im, 0, x_max, 1) - get_pixel(im, y_max - 1, x_max, 1);
    B_diff_y = get_pixel(im, 0, x_max, 2) - get_pixel(im, y_max - 1, x_max, 2);

    delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
    delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

    energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
    set_pixel(*grad, y_max, x_max, energy, energy, energy);

    // edges 

    for (int m = 1; m < x_max; m++) {

        R_diff_x = get_pixel(im, 0, m+1, 0) - get_pixel(im, 0, m-1, 0);
        G_diff_x = get_pixel(im, 0, m+1, 1) - get_pixel(im, 0, m-1, 1);
        B_diff_x = get_pixel(im, 0, m+1, 2) - get_pixel(im, 0, m-1, 2);

        R_diff_y = get_pixel(im, 1, m, 0) - get_pixel(im, y_max, m, 0);
        G_diff_y = get_pixel(im, 1, m, 1) - get_pixel(im, y_max, m, 1);
        B_diff_y = get_pixel(im, 1, m, 2) - get_pixel(im, y_max, m, 2);

        delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
        delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

        energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
        set_pixel(*grad, 0, m, energy, energy, energy);

        R_diff_x = get_pixel(im, y_max, m+1, 0) - get_pixel(im, y_max, m-1, 0);
        G_diff_x = get_pixel(im, y_max, m+1, 1) - get_pixel(im, y_max, m-1, 1);
        B_diff_x = get_pixel(im, y_max, m+1, 2) - get_pixel(im, y_max, m-1, 2);

        R_diff_y = get_pixel(im, 0, m, 0) - get_pixel(im, y_max - 1, m, 0);
        G_diff_y = get_pixel(im, 0, m, 1) - get_pixel(im, y_max - 1, m, 1);
        B_diff_y = get_pixel(im, 0, m, 2) - get_pixel(im, y_max - 1, m, 2);

        delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
        delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

        energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
        set_pixel(*grad, y_max, m, energy, energy, energy);
    }

    for (int n = 1; n < y_max; n++) {
        R_diff_x = get_pixel(im, n, 1, 0) - get_pixel(im, n, x_max, 0);
        G_diff_x = get_pixel(im, n, 1, 1) - get_pixel(im, n, x_max, 1);
        B_diff_x = get_pixel(im, n, 1, 2) - get_pixel(im, n, x_max, 2);

        R_diff_y = get_pixel(im, n+1, 0, 0) - get_pixel(im, n-1, 0, 0);
        G_diff_y = get_pixel(im, n+1, 0, 1) - get_pixel(im, n-1, 0, 1);
        B_diff_y = get_pixel(im, n+1, 0, 2) - get_pixel(im, n-1, 0, 2);

        delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
        delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

        energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
        set_pixel(*grad, n, 0, energy, energy, energy);

        R_diff_x = get_pixel(im, n, 0, 0) - get_pixel(im, n, x_max - 1, 0);
        G_diff_x = get_pixel(im, n, 0, 1) - get_pixel(im, n, x_max - 1, 1);
        B_diff_x = get_pixel(im, n, 0, 2) - get_pixel(im, n, x_max - 1, 2);

        R_diff_y = get_pixel(im, n+1, x_max, 0) - get_pixel(im, n-1, x_max, 0);
        G_diff_y = get_pixel(im, n+1, x_max, 1) - get_pixel(im, n-1, x_max, 1);
        B_diff_y = get_pixel(im, n+1, x_max, 2) - get_pixel(im, n-1, x_max, 2);

        delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
        delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

        energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
        set_pixel(*grad, n, x_max, energy, energy, energy);
    }

    for (int i = 1; i < y_max; i++) {
        for (int j = 1; j < x_max; j++) {
            R_diff_x = get_pixel(im, i, j+1, 0) - get_pixel(im, i, j-1, 0);
            G_diff_x = get_pixel(im, i, j+1, 1) - get_pixel(im, i, j-1, 1);
            B_diff_x = get_pixel(im, i, j+1, 2) - get_pixel(im, i, j-1, 2);

            R_diff_y = get_pixel(im, i+1, j, 0) - get_pixel(im, i-1, j, 0);
            G_diff_y = get_pixel(im, i+1, j, 1) - get_pixel(im, i-1, j, 1);
            B_diff_y = get_pixel(im, i+1, j, 2) - get_pixel(im, i-1, j, 2);

            delta_x = pow(R_diff_x, 2) + pow(G_diff_x, 2) + pow(B_diff_x, 2);
            delta_y = pow(R_diff_y, 2) + pow(G_diff_y, 2) + pow(B_diff_y, 2);

            energy = (uint8_t) (pow(delta_x + delta_y, 0.5)/10);
            set_pixel(*grad, i, j, energy, energy, energy);
        }
    }
}

void dynamic_seam(struct rgb_img *grad, double **best_arr){
    
    int width = grad -> width;
    int height = grad -> height;

    *best_arr = (double *)malloc(height*width*sizeof(double));

    double root, left, right, middle, min;

    
    //inputs numbers in first row into the matrix
    for (int a = 0; a< width; a++){

        (*best_arr)[a] = get_pixel(grad, 0, a, 0);

        //printf("%f  ", (*best_arr)[a]);

    }

                

    //goes through the rows
    for (int i = 1; i<height; i++){
    
        //goes through the columns
        for (int j = 0; j< width; j++){


            left = 100000;
            right = 100000;

            root = get_pixel(grad, i, j, 0);

            //compute the sum of each path
            if (j == 0){
                middle = (*best_arr)[(i-1)*width+(j)];
                right = (*best_arr)[(i-1)*width+(j+1)];

            }
            else if (j == width-1){
                left = (*best_arr)[(i-1)*width+(j-1)];
                middle = (*best_arr)[(i-1)*width+(j)];

            }
            else {
                left = (*best_arr)[(i-1)*width+(j-1)];
                middle = (*best_arr)[(i-1)*width+(j)];
                right = (*best_arr)[(i-1)*width+(j+1)];
            }

            min = find_min_dbl(left, middle, right); //retrieves smallest path

            //assigns the lowest calculated cost
            (*best_arr)[i*width+j] = min + root;  
            //printf("%f  ", (*best_arr)[i*width+j]);

        }
        
    }

}


void recover_path(double *best, int height, int width, int **path){

    *path = (int*)malloc(height*sizeof(int));
    int size = height*width;
    double min = 100000; 
    int x = 0; 

    for (int i = 1; i<=width; i++){

        //checks for smallest cost from bottom row
        if (min>best[size-i]){
            min = best[size-i];
            x = width-i;  
        }
    }

    int ind = height-1; //index where the path value goes

    (*path)[ind] = x; //sets the last index = the index of the lowest cost


    double left, middle, right;

    //starting at row second from the top
    for (int j = height-2; j>=0; j--){
        left = 100000;
        right = 100000;
        middle = 100000;


        //compute the sum of each path
        if (x == 0){
            middle = (best)[j*width + x]; 
            right = (best)[j*width + x+1];
            //printf("why are you not working");

        }
        else if (x == width-1){
            left = (best)[j*width + x-1]; 
            middle = (best)[j*width + x]; 
            //printf("this is happening");

        }
        else {
            left = (best)[j*width + x-1]; 
            middle = (best)[j*width + x]; 
            right = (best)[j*width + x+1];
            //printf("bruh");
        }

        min = find_min_dbl(left, middle, right);

        if (min == left){
            x = x-1;
        } else if (min == right){
            x = x+1;
        }
        //printf("fuck my life");

        ind--; // goes to next index in path

        (*path)[ind] = x;
        //printf("stoopid index: %d", x);
    }
}

void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path){
    int height = src -> height;
    int width = src -> width - 1;
    create_img(dest, height, width);
    int dest_index = 0;
    int r, g, b;

    for (int i = 0; i < height; i++) {
        dest_index = 0;

        for (int j = 0; j < (src -> width); j++){

            if (j == path[i]) {
                j++;
            }
            
            if (j != (src -> width)) {
                r = get_pixel(src, i, j, 0);
                g = get_pixel(src, i, j, 1);
                b = get_pixel(src, i, j, 2);
                set_pixel(*dest, i, dest_index, r, g, b);
                dest_index++;
            }
       
        }

    }

}

void print_my_array(double a[]){
    int i;
    for(i=0;i<*a;i++){
        printf("%f  ",a[i]);

    }
}

/*

int main(void) {

    char filename[] = "6x5.bin";

    struct rgb_img *ima;
    
    read_in_img(&ima, filename);

    struct rgb_img *grad;

    create_img(&grad, ima -> height, ima -> width);

    calc_energy(ima, &grad);
  
    print_grad(grad);

    double *best_arr;

    dynamic_seam(grad, &best_arr);

    //print_my_array(&best_arr[0]);

    int height = ima -> height ;
    int width = ima -> width ;
    int **path;
    recover_path(best_arr, height, width, &path);

    int loop;
    for(loop = 0; loop < height; loop++){
        printf("%d \n", path[loop]);
        
    }

    struct rgb_img *im;
    struct rgb_img *cur_im;
    struct rgb_img *grad;
    double *best;
    int *path;

    read_in_img(&im, "HJoceanSmall.bin");
    
    for(int i = 0; i < 5; i++){
        printf("i = %d\n", i);
        calc_energy(im,  &grad);
        dynamic_seam(grad, &best);
        recover_path(best, grad->height, grad->width, &path);
        remove_seam(im, &cur_im, path);

        char filename[200];
        sprintf(filename, "img%d.bin", i);
        write_img(cur_im, filename);


        destroy_image(im);
        destroy_image(grad);
        free(best);
        free(path);
        im = cur_im;
    }
    destroy_image(im);
}
*/