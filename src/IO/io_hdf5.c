#define IO_PRIVATE_DEFS
#define H5FILE_NAME     "testout.h5"
#define DATASETNAME 	"DoubleArray" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/Grid.h"
#include "../Headers/Cell.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

void io_hdf5_out(struct IO *io_pointer,struct Grid * theGrid, char * output_filename){
  int NUM_Q = grid_NUM_Q(theGrid);
  int Ncells = grid_Ncells(theGrid);
  double prim_data[Ncells][NUM_Q+3];
  int i,q;
  for (i=0;i<Ncells;++i){
      prim_data[i][0] = io_pointer->primitives[i][0];
      prim_data[i][1] = io_pointer->primitives[i][1];
      prim_data[i][2] = io_pointer->primitives[i][2];
    for (q=0;q<NUM_Q;++q){
      prim_data[i][q+3] = io_pointer->primitives[i][q+3];
    }
  }

  // HDF5 APIs definitions
  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t     dimsf[2];                 /* dataset dimensions */
  hsize_t     chunk_dims[2];            /* chunk dimensions */
  hsize_t	count[2];	          /* hyperslab selection parameters */
  hsize_t	stride[2];
  hsize_t	block[2];
  hsize_t	offset[2];
  hid_t	plist_id;                 /* property list identifier */
  herr_t	status;

  // MPI variables
  MPI_Info info  = MPI_INFO_NULL;

  // Set up file access property list with parallel I/O access
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, grid_comm, info);

  // Create a new file collectively and release property list identifier.
  file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // Create the dataspace for the dataset.
  dimsf[0] = grid_Ncells_global(theGrid);
  dimsf[1] = NUM_Q+3;
  chunk_dims[0] = grid_Ncells(theGrid);
  chunk_dims[1] = NUM_Q+3;
  filespace = H5Screate_simple(RANK, dimsf, NULL); 
  memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 

  // Create chunked dataset.
  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  //H5Pset_chunk(plist_id, RANK, chunk_dims);
  //dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  dset_id = H5Dcreate1(file_id, DATASETNAME, H5T_NATIVE_DOUBLE, filespace,plist_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);

  // Each process defines dataset in memory and writes it to the hyperslab in the file.
  count[0] = 1;
  count[1] = 1;
  stride[0] = 1;
  stride[1] = 1;
  block[0] = chunk_dims[0];
  block[1] = chunk_dims[1];
  offset[0] = grid_offset(theGrid);
  offset[1] = 0;

  // Select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, prim_data);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Fclose(file_id);
}    

void io_hdf5_in(struct IO *io_pointer,struct Grid * theGrid){
  int NUM_Q = grid_NUM_Q(theGrid);
  int Ncells = grid_Ncells(theGrid);
  double prim_data[Ncells][NUM_Q+3];
   // double **chunk_out = io_pointer->primitives;
  hid_t       file;                        /* handles */
  hid_t       dataset;  
  hid_t       filespace;                   
  hid_t       memspace;                  
  hsize_t     dims[2];                     /* dataset and chunk dimensions*/ 
  hsize_t     chunk_dims[2];
  hsize_t     count[2];
  hsize_t     offset[2];
  hsize_t     block[2];
  hsize_t     stride[2];
  herr_t      status, status_n;                             
  int         rank, rank_chunk;
  //int         i, j,k;

  // Open the file and the dataset.
  file = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen1(file, "DoubleArray");

  // Get dataset rank and dimension.
  filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
  rank      = H5Sget_simple_extent_ndims(filespace);
  status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);

  chunk_dims[0] = grid_Ncells(theGrid);
  chunk_dims[1] = NUM_Q+3;

  // Define the memory space to read a chunk.
  memspace = H5Screate_simple(RANK,chunk_dims,NULL);

  count[0] = 1;
  count[1] = 1;
  stride[0] = 1;
  stride[1] = 1;
  block[0] = chunk_dims[0];
  block[1] = chunk_dims[1];

  offset[0] = grid_offset(theGrid);
  offset[1] = 0; 

  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

  // Read chunk 
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, prim_data);

  int i,q;
  for (i=0;i<grid_Ncells(theGrid);++i){
    io_pointer->primitives[i][0] = prim_data[i][0] ;
    io_pointer->primitives[i][1] = prim_data[i][1] ;
    io_pointer->primitives[i][2] = prim_data[i][2] ;
    for (q=0;q<NUM_Q;++q){
      io_pointer->primitives[i][q+3] = prim_data[i][q+3] ;
    }
  }


  // Close/release resources.
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Fclose(file);
  int index;
  if (1==0){
    for (index=0;index<grid_N_r(theGrid);++index){
      printf("io_pointer->primitives[index][0]: %f\n",io_pointer->primitives[index][0]);
    }
  }
}
