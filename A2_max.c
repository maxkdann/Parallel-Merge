#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef struct {
    int *local_a_sizes;
    int *local_b_sizes;
    int *local_a_indices;
    int *local_b_indices;
} array_info;

/**
 * Prints an array for debugging purposes
 * Parameters:
 * arr: array to be printed
 * size: size of arr
*/
void print_arr(const int arr[], const int size){
    for(int i=0;i<size;i++){
        printf("%d, ", arr[i]);
    }
    printf("\n");
}

/**
 * find the sizes of the local_a arrays, they will differ by at most 1
 * Parameters:
 * local_a_sizes: pointer to array storing the size of each local array
 * num_procs: number of processors
 * ARR_SIZE: total number of elements in array A (and array B)
*/
void find_local_a_sizes(int *local_a_sizes, int num_procs, int ARR_SIZE){
    //load balancing
    int remainder = ARR_SIZE % num_procs;
    for(int i =num_procs-1;i>-1;i--){
        local_a_sizes[i] = ARR_SIZE / num_procs;
        if(remainder>0){
            local_a_sizes[i]+=1;
            --remainder;
        }
    }
}

/**
 * find the indices of A for each local_a array
 * local_a_indices: pointer to array containing the indices for the start of each local_a in a
 * local_a_sizes: pointer to array storing the size of each local array
 * num_procs: number of processors
 * 
*/
void find_local_a_indices(int *local_a_indices,int *local_a_sizes, int num_procs){
    int start = 0;
    int index;
    local_a_indices[0] = start;
    for (int i = 1; i < num_procs; i++) {
      index = local_a_sizes[i-1] + local_a_indices[i-1];
      local_a_indices[i] = index;
    }
}

/**
 * Perform binary search on the array to find the greatest index i such that arr[i] <= x
 * Parameters:
 * arr: array to be searched
 * l: left-most element of arr
 * r: right-most element of arr
 * x: maximum element
*/
int binary_search(int arr[], int l, int r, int x) {
    int m;
    while (l <= r) {
        m = l + (r-l)/2;
        if (arr[m] == x) {
            return m;
        } else if (arr[m] < x) {
            l = m + 1;
        } else {
            r = m - 1;
        }
    }
    return r; // Return the index of the greatest element that is less than or equal to x
}


/**
 * find the starting indices of each local_b array within the larger b array as well
 * as the size of each of those local_b arrays
 * Parameters:
 * local_a_sizes: pointer to array storing the size of each local array
 * num_procs: number of processors
 * ARR_SIZE: total number of elements in array A (and array B)
 * a: the total a array
 * b: the total b array
 * local_a_indices: pointer to array containing the indices for the start of each local_a in a
 * local_b_indices: pointer to array containing the indices for the start of each local_b in b
 * local_b_sizes: pointer to array storing the size of each local array
*/
void find_b_info(int num_procs,int *local_a_indices,int *local_a_sizes,int *local_b_indices,int *local_b_sizes,int ARR_SIZE, int *a, int *b){
    int key_index;
    int key;
    int size;
    int curr=0;
    int index;
    for(int i = 0;i<num_procs;i++){
        key_index = local_a_sizes[i] + local_a_indices[i]-1;
        key = a[key_index];
        //find the index of the greatest element that is less than or equal to key
        index = binary_search(b,0,ARR_SIZE-1,key);
        size = index +1 - curr;
        local_b_indices[i] = curr;
        local_b_sizes[i] = size;
        curr = index +1;
        index+=local_a_sizes[i+1];
    }
}


/**
 * partitions arrays a and b for parallel merging
 * local_a_sizes: pointer to array storing the size of each local array
 * num_procs: number of processors
 * ARR_SIZE: total number of elements in array A (and array B)
 * a: the total a array
 * b: the total b array
 * local_a_indices: pointer to array containing the indices for the start of each local_a in a
 * local_b_indices: pointer to array containing the indices for the start of each local_b in b
 * local_b_sizes: pointer to array storing the size of each local array
*/
void partition(int ARR_SIZE, int* a, int* b, int num_procs, array_info *data) {
    //initialize arrays to return data
    int *arr_a_size, *arr_b_size,*arr_a_indices,*arr_b_indices;

    //malloc space for each array
    arr_a_size = (int*)malloc(num_procs * sizeof(int));
    arr_b_size = (int*)malloc(num_procs * sizeof(int));
    arr_a_indices = (int*)malloc(num_procs * sizeof(int));
    arr_b_indices = (int*)malloc(num_procs * sizeof(int));
    
    //find sizes of local a arrays
    find_local_a_sizes(arr_a_size,num_procs,ARR_SIZE);

    //determine a_indices
    find_local_a_indices(arr_a_indices,arr_a_size,num_procs);
    
    //find b info
    find_b_info(num_procs,arr_a_indices,arr_a_size,arr_b_indices,arr_b_size,ARR_SIZE,a,b);

    size_t size_of_array = (sizeof *data->local_a_sizes) * num_procs;
   
    // malloc memory for each attribute of array_data struct
    data->local_a_sizes = malloc(size_of_array);
    data->local_b_sizes = malloc(size_of_array);
    data->local_a_indices = malloc(size_of_array);
    data->local_b_indices = malloc(size_of_array);
    
    // initialize array_data struct with appropriate data
    data->local_a_sizes = arr_a_size;
    data->local_b_sizes = arr_b_size;
    data->local_a_indices = arr_a_indices;
    data->local_b_indices = arr_b_indices;
    
}


/**
 * Merges two arrays local_a and local_b and stores the result in local_c
 * Parameters:
 * local_a: i-th partition of array A
 * local_b: i-th partition of array B
 * local_c: array to store the resulting merged elements
 * size_a: the size of local_a
 * size_b: the size of local_b
*/
void seq_merge(int local_a[], int local_b[], int local_c[], int size_a, int size_b) {
    int i = 0, j = 0, k = 0;

    while (i < size_a && j < size_b) {
        if (local_a[i] < local_b[j]) {
            local_c[k++] = local_a[i++];
        }
        else {
            local_c[k++] = local_b[j++];
        }
    }

    while (i < size_a) {
        local_c[k++] = local_a[i++];
    }

    while (j < size_b) {
        local_c[k++] = local_b[j++];
    }

}

/**
 * generates two sorted arrays for testing
 * n: size of the arrays
 * arr1: first empty array
 * arr2: second empty array
*/
void generate_sorted_arrays(int n, int arr1[], int arr2[]) {
    int i;
    srand(time(NULL)); // Initialize the random seed

    // Generate the first array
    arr1[0] = rand() % 10; // First element is random
    for (i = 1; i < n; i++) {
        arr1[i] = arr1[i-1] + (rand() % 10 + 1); // Next element is +1 to +10 of the previous element
    }

    // Generate the second array
    arr2[0] = rand() % 10; // First element is random
    for (i = 1; i < n; i++) {
        arr2[i] = arr2[i-1] + (rand() % 10 + 1); // Next element is +1 to +10 of the previous element
    }
}




int main(int argc, char *argv[]){
    //initialize variables
    int ARR_SIZE = 10;
    int *a;
    int *b;
    MPI_Status status;
    MPI_Request request;
    int process_rank;
    int num_procs;

    array_info *array_data = (array_info *)malloc(sizeof(array_info));

    //start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double start_time = 0.0;
    double end_time = 0.0;

    if(process_rank==0){
        //allocate memory for two arrays
        a = malloc(ARR_SIZE*sizeof(int));
        b = malloc(ARR_SIZE*sizeof(int));
        //generate sorted arrays to be stored in a and b
        generate_sorted_arrays(ARR_SIZE,a,b);
        //partition the data
        partition(ARR_SIZE,a,b,num_procs,array_data);
    }
    else{
        //allocate memory
        array_data->local_a_sizes = (int*)malloc(num_procs * sizeof(int));
        array_data->local_b_sizes = (int*)malloc(num_procs * sizeof(int));
        array_data->local_a_indices = (int*)malloc(num_procs * sizeof(int));
        array_data->local_b_indices= (int*)malloc(num_procs * sizeof(int));
    }
    //wait for all procs to be ready
    MPI_Barrier(MPI_COMM_WORLD);
    //Broadcast info to all other processes
    MPI_Bcast(array_data->local_a_sizes,num_procs,MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(array_data->local_a_indices,num_procs,MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(array_data->local_b_sizes,num_procs,MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(array_data->local_b_indices,num_procs,MPI_INT, 0, MPI_COMM_WORLD);

    //get the length of the local subarrays
    int size_of_local_a = array_data->local_a_sizes[process_rank];
    int size_of_local_b = array_data->local_b_sizes[process_rank];

    //start the clock
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    //allocate space for the arrays on each processor depending on size of local a and b
    int *local_a = (int *) malloc(sizeof(int) * size_of_local_a);
    int *local_b = (int *) malloc(sizeof(int) * size_of_local_b);

    //scatter A and B across all processes, noting that sizes can differ
    MPI_Scatterv(a,array_data->local_a_sizes,array_data->local_a_indices,MPI_INT,local_a,size_of_local_a,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Scatterv(b,array_data->local_b_sizes,array_data->local_b_indices,MPI_INT,local_b,size_of_local_b,MPI_INT,0,MPI_COMM_WORLD);

    //merge local_a and local_b into local_c
    int size_of_local_c = size_of_local_a+size_of_local_b;
    int *local_c = (int *)malloc(sizeof(int) * size_of_local_c);
    seq_merge(local_a,local_b,local_c,size_of_local_a,size_of_local_b);
    
    //allocate memory for c - merged total array of a and b
    int *c = (int *)malloc(sizeof(int)*(ARR_SIZE*2));
    int *local_c_sizes = (int *)malloc(sizeof(int)*num_procs);
    int *local_c_indices = (int *)malloc(sizeof(int) * num_procs);

    //populate each instance of local_c
    for(int i=0;i<num_procs;i++){
        local_c_sizes[i] = array_data->local_a_sizes[i]+array_data->local_b_sizes[i];
        local_c_indices[i] = array_data->local_a_indices[i]+array_data->local_b_indices[i];
    }

    //gather all local_c instances back to the root process, keeping in mind they can be of different sizes
    MPI_Gatherv(local_c,size_of_local_c,MPI_INT, c, local_c_sizes,local_c_indices,MPI_INT,0,MPI_COMM_WORLD);
    //stop the clock
    end_time = MPI_Wtime();
    
    //print output
    if(process_rank==0){
        printf("Array A: ");
        print_arr(a, ARR_SIZE);
        
        printf("Array B: ");
        print_arr(b, ARR_SIZE);
        
        printf("Array C: ");
        print_arr(c, ARR_SIZE * 2);
    
        printf("Wallclock time elapsed: %.9lf seconds\n", end_time - start_time);
    
        //free memory
        free(a);
        free(b);
    }
    free(array_data->local_a_sizes);
    free(array_data->local_a_indices);  
    free(local_a);
    free(array_data->local_b_sizes);
    free(array_data->local_b_indices);
    free(local_b);
    free(local_c_sizes);
    free(local_c_indices);
    free(local_c);
    free(c);

    MPI_Finalize();
    return 0;
}
