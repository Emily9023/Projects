#include "c_img.h"
#include "seamcarving.h"


int main(void) {
    struct rgb_img *im;
    struct rgb_img *cur_im;
    struct rgb_img *grad;
    double *best;
    int *path;

    read_in_img(&im, "test8.bin");

    //reduce by 400 pixels
    int loop = 400;
    for(int i = 0; i < loop; i++){
        printf("i = %d\n", i);
        calc_energy(im,  &grad);
        dynamic_seam(grad, &best);
        recover_path(best, grad->height, grad->width, &path);
        remove_seam(im, &cur_im, path);

        //once the it reduces all the pixels, writes it to a bin file
        if (i == loop-1){
            write_img(cur_im, "img8.bin");

        }
        destroy_image(im);
        destroy_image(grad);
        free(best);
        free(path);
        im = cur_im;
    }
    destroy_image(im);
 
    

    return 0;

}
