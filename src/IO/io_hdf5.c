#define IO_PRIVATE_DEFS
#define DATASETNAME 	"Data" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

#ifdef CHECKPOINTING
void io_hdf5_out(struct IO *theIO,struct Sim * theSim,struct TimeStep * theTimeStep){
  int NUM_Q = sim_NUM_Q(theSim);
  int Ncells = sim_Ncells(theSim);
  int i,q;
  // HDF5 APIs definitions
  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t     dimsf1[1];                 /* 1d dataset dimensions */
  hsize_t     dimsf2[2];                 /* 2d dataset dimensions */
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
  H5Pset_fapl_mpio(plist_id, sim_comm, info);

  // Create a new file collectively and release property list identifier.
  file_id = H5Fcreate(theIO->filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // **********************
  // First we save the time
  // **********************

  // we will only be saving one element
  dimsf1[0] = 1;

  // Create the dataspace for the dataset.
  filespace = H5Screate_simple(1, dimsf1, NULL); 
  memspace  = H5Screate_simple(1, dimsf1, NULL); 

  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  dset_id = H5Dcreate1(file_id, "T", H5T_NATIVE_DOUBLE, filespace,plist_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);

  // Select hyperslab in the file.
  filespace = H5Dget_space(dset_id);

  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, prim_data);
  double time = timestep_get_t(theTimeStep);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &(time));

  // Close/release resources.
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Dclose(dset_id);

  // *****************************************
  // Then we save the info of the GravMasses
  // *****************************************

  dimsf2[0] = 2;  // number of masses
  dimsf2[1] = 11; //number of quantities per mass DD!!!  4;  
  // Create the dataspace for the dataset.
  filespace = H5Screate_simple(2, dimsf2, NULL); 
  memspace  = H5Screate_simple(2, dimsf2, NULL); 

  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  dset_id = H5Dcreate1(file_id, "GravMass", H5T_NATIVE_DOUBLE, filespace,plist_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);

  filespace = H5Dget_space(dset_id);

  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, prim_data);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, theIO->GravMassBuffer);

  // Close/release resources.
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Dclose(dset_id);

  // ***************************
  // Now we save the actual data
  // ***************************

  // Create the dataspace for the dataset.
  dimsf2[0] = sim_Ncells_global(theSim);
  dimsf2[1] = NUM_Q+3;
  chunk_dims[0] = sim_Ncells(theSim);
  chunk_dims[1] = NUM_Q+3;
  filespace = H5Screate_simple(RANK, dimsf2, NULL); 
  memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 

  // Create chunked dataset.
  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  //H5Pset_chunk(plist_id, RANK, chunk_dims);
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
  offset[0] = sim_offset(theSim);
  offset[1] = 0;

  // Select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, prim_data);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &(theIO->buffer[0][0]));

  // Close/release resources.
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Fclose(file_id);

  theIO->tcheck += theIO->dtcheck;
  ++(theIO->nfile);
  sprintf(theIO->filename,"checkpoint_%04d.h5",io_nfile(theIO));

}    

void io_hdf5_in(struct IO *theIO,struct Sim * theSim,struct TimeStep * theTimeStep){
  int NUM_Q = sim_NUM_Q(theSim);
  int Ncells = sim_Ncells(theSim);
  hid_t       file;                        /* handles */
  hid_t       dataset;  
  hid_t       filespace;                   
  hid_t       memspace;                  
  hsize_t     dims1[1];                     // 1d dataset dimensions 
  hsize_t     dims2[2];                     // 2d dataset dimensions 
  hsize_t     chunk_dims[2];
  hsize_t     count[2];
  hsize_t     offset[2];
  hsize_t     block[2];
  hsize_t     stride[2];
  herr_t      status, status_n;                             
  int         rank, rank_chunk;
  //int         i, j,k;

  // Open the file and the dataset.
  file = H5Fopen("input.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  // **********************
  // First read in the time
  // **********************

  dataset = H5Dopen1(file,"T");
  dims1[0] = 1;
  // Get dataset rank and dimension.
  filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
  rank      = H5Sget_simple_extent_ndims(filespace);

  // Define the memory space to read a chunk.
  memspace = H5Screate_simple(rank,dims1,NULL);

  double time = 0.0;
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT,&time);

  timestep_set_t(theTimeStep,time);

  // Close/release resources.
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);

  // ********************************************
  // Then we read in the info of the 1st GravMass
  // ********************************************

  dataset = H5Dopen1(file,"GravMass");
  dims2[0] = 2;
  dims2[1] = 9; //DD!!!  4;

  // Get dataset rank and dimension.
  filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
  rank      = H5Sget_simple_extent_ndims(filespace);

  // Define the memory space to read a chunk.
  memspace = H5Screate_simple(rank,dims2,NULL);

  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT,theIO->GravMassBuffer);

  // Close/release resources.
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);



  // **********************
  // Then we read in the actual data
  // **********************

  dataset = H5Dopen1(file, DATASETNAME);

  // Get dataset rank and dimension.
  filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
  rank      = H5Sget_simple_extent_ndims(filespace);
  status_n  = H5Sget_simple_extent_dims(filespace, dims2, NULL);

  chunk_dims[0] = sim_Ncells(theSim);
  chunk_dims[1] = NUM_Q+3;

  // Define the memory space to read a chunk.
  memspace = H5Screate_simple(RANK,chunk_dims,NULL);

  count[0] = 1;
  count[1] = 1;
  stride[0] = 1;
  stride[1] = 1;
  block[0] = chunk_dims[0];
  block[1] = chunk_dims[1];

  offset[0] = sim_offset(theSim);
  offset[1] = 0; 

  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

  // Read chunk 
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT,&(theIO->buffer[0][0]));

  // Close/release resources.
  H5Dclose(dataset);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Fclose(file);
}
#endif
