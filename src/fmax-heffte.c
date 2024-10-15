/*****************************************************************
 *                        PINOCCHIO  V4.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco
 Copyright (C) 2016
 
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

//#define DEBUG

#ifdef DEBUG
  FILE *results;
  char filename[300];
#endif


double greens_function(double *, double, int, int);


int set_one_grid(int ThisGrid)
{
  int const Nx = MyGrids[ThisGrid].GSglobal_x;
  int const Ny = MyGrids[ThisGrid].GSglobal_y;
  int const Nz = MyGrids[ThisGrid].GSglobal_z;

  // Distribute the planes as evenly as possible across the MPI tasks.
  // The smaller task numbers will get one extra plane at most.
  // Note sure this is correct for Ny != Nz or if the z-component in real
  // space actually corresponds to the y-component in Fourier space when
  // using HeFFTe, as it does with FFTW.
  int nz, z0;
  int nk, k0;
  int extra_z = Nz % NTasks;
  int extra_k = Ny % NTasks;
  if (ThisTask < extra_z) {
    nz = Nz / NTasks + 1;
    z0 = ThisTask * nz;
  } else {
    nz = Nz / NTasks;
    z0 = extra_z * (nz + 1) + (ThisTask - extra_z) * nz;
  }
  if (ThisTask < extra_k) {
    nk = Ny / NTasks + 1;
    k0 = ThisTask * nk;
  } else {
    nk = Ny / NTasks;
    k0 = extra_k * (nk + 1) + (ThisTask - extra_k) * nk;
  }

  MyGrids[ThisGrid].GSstart_x = 0;
  MyGrids[ThisGrid].GSstart_y = 0;
  MyGrids[ThisGrid].GSstart_z = z0;
  MyGrids[ThisGrid].GSlocal_x = Nx;
  MyGrids[ThisGrid].GSlocal_y = Ny;
  MyGrids[ThisGrid].GSlocal_z = nz;

  MyGrids[ThisGrid].GSstart_k_x = 0;
  MyGrids[ThisGrid].GSstart_k_y = k0;
  MyGrids[ThisGrid].GSstart_k_z = 0;
  MyGrids[ThisGrid].GSlocal_k_x = Nx;
  MyGrids[ThisGrid].GSlocal_k_y = nk;
  MyGrids[ThisGrid].GSlocal_k_z = Nz;

  // This gives the number of unused x-grid points in rvector_fft.
  if (Nx % 2)
    MyGrids[ThisGrid].off = 1;
  else
    MyGrids[ThisGrid].off = 2;

  MyGrids[ThisGrid].norm = 1.0 / (Nx * Ny * Nz);
  MyGrids[ThisGrid].CellSize = MyGrids[ThisGrid].BoxSize / Nx;

  int const input_indices_start[3]  = {   0,    0, z0};
  int const input_indices_end[3]    = {Nx-1, Ny-1, nz-1};
  int const input_index_order[3]    = {0, 1, 2};
  int const output_indices_start[3] = {   0,   k0, 0};
  int const output_indices_end[3]   = {Nx/2, nk-1, Nz-1};
  int const output_index_order[3]   = {0, 2, 1};
  int const r2c_direction = 0;
  int const fft_backend = Heffte_BACKEND_STOCK; // Heffte_BACKEND_CUFFT, Heffte_BACKEND_FFTW
  heffte_plan_options plan_options;
  heffte_set_default_options(fft_backend, &plan_options);
  // plan_options.use_reorder = 0;

  heffte_plan_create_r2c(
    fft_backend,
    input_indices_start, input_indices_end,
    input_index_order,
    output_indices_start, output_indices_end,
    output_index_order,
    r2c_direction,
    MPI_COMM_WORLD,
    &plan_options,
    &MyGrids[ThisGrid].fft_plan);

  MyGrids[ThisGrid].total_local_size     = heffte_size_inbox(MyGrids[ThisGrid].fft_plan);
  MyGrids[ThisGrid].total_local_size_fft = heffte_size_outbox(MyGrids[ThisGrid].fft_plan);

  return 0;
}


double forward_transform(int ThisGrid)
{
  double time = MPI_Wtime();
  heffte_forward_d2z(
    MyGrids[ThisGrid].fft_plan,
    rvector_fft[ThisGrid], cvector_fft[ThisGrid],
    Heffte_SCALE_NONE);
  return MPI_Wtime() - time;
}


double reverse_transform(int ThisGrid)
{
  double time = MPI_Wtime();
  heffte_backward_z2d(
    MyGrids[ThisGrid].fft_plan,
    cvector_fft[ThisGrid], rvector_fft[ThisGrid],
    Heffte_SCALE_NONE);
  for (int i = 0; i < MyGrids[ThisGrid].total_local_size_fft; i++)
      rvector_fft[ThisGrid][i] *= MyGrids[ThisGrid].norm;
#ifdef DEBUG
  for (int i = 0; i < MyGrids[ThisGrid].total_local_size_fft; i++)
    fprintf(results, " %d  %g \n", i, rvector_fft[ThisGrid][i]);
  fclose(results);
#endif
  return MPI_Wtime() - time;
}


int finalize_fft()
{
  for (int igrid = Ngrids - 1; igrid >= 0; igrid--)
    heffte_plan_destroy(MyGrids[igrid].fft_plan);
  for (int igrid = Ngrids - 1; igrid >= 0; igrid--)
    if (deallocate_fft_vectors(igrid))
      return 1;
  return 0;
}


int compute_derivative(int ThisGrid, int first_derivative, int second_derivative)
{
  int swap, local_x,local_y,local_z,ixx,iyy,izz,index,nxhalf,nyhalf,nzhalf;
  double kx,ky,kz,kxnorm,kynorm,kznorm,green,k_squared,smoothing,tmp,time;
  double diff_comp[4];

#ifdef DEBUG
  sprintf(filename,"results.%d-%d.%d",first_derivative,second_derivative,ThisTask);
  results=fopen(filename,"w");
#endif

  /* k vectors */
  kxnorm = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_x;
  kynorm = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_y;
  kznorm = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_z;

  /* Nyquist frequencies */
  nxhalf = MyGrids[ThisGrid].GSglobal_x/2;
  nyhalf = MyGrids[ThisGrid].GSglobal_y/2;
  nzhalf = MyGrids[ThisGrid].GSglobal_z/2;

  /* for first derivatives the real and imaginary parts must be swapped */
  swap=((first_derivative==0 && second_derivative> 0) ||
	(first_derivative> 0 && second_derivative==0));

/*
  NB: ix, iy and iz are global coordinates of the box (from 0 to N-1)
      ix: [0,N/2]
      iy, iz: [0,N-1]
      iyy and izz are unfolded, they run from -N/2+1 to N/2 (ixx = ix)
      local_x, local_y and local_z are local coordinates of the slice
      (for the x and z coordinates they are the same as ix and iz)
*/

/* loop over k-space indices */
/* This loop is correct for transposed order in FFTW */

  for (local_z = 0; local_z < MyGrids[ThisGrid].GSlocal_k_z; local_z++)
    {
      izz = local_z + MyGrids[ThisGrid].GSstart_k_z;
      if (local_z > nzhalf)
	izz -= MyGrids[ThisGrid].GSglobal_z;
      kz  = kznorm*izz;

      for (local_y = 0; local_y < MyGrids[ThisGrid].GSlocal_k_y; local_y++)
	{
	  iyy = local_y + MyGrids[ThisGrid].GSstart_k_y;
	  if (iyy > nyhalf)
	    iyy -= MyGrids[ThisGrid].GSglobal_y;
	  ky  = kynorm*iyy;

	  for (local_x = 0; local_x <= nxhalf; local_x++)
	    {
	      ixx = local_x;
	      kx  = kxnorm*ixx;

              k_squared  = kx*kx + ky*ky + kz*kz;

	      /* corresponding index of real part in vector (imaginary in index + 1) */
	      index = 1 + 2*local_x + (MyGrids[ThisGrid].GSglobal_x+MyGrids[ThisGrid].off)
		*(local_z + local_y* MyGrids[ThisGrid].GSglobal_z);

              if (k_squared != 0.)
		{

		  /* Gaussian smoothing window */
		  smoothing = exp(-0.5 * k_squared * Rsmooth * Rsmooth);

		  /* the components are stored in the vectors that control the differentiation */
		  diff_comp[0] = 1.0;
		  diff_comp[1] = kx;
		  diff_comp[2] = ky;
		  diff_comp[3] = kz;

		  green = greens_function(diff_comp, k_squared, first_derivative, second_derivative);

		  (cvector_fft[ThisGrid][index/2]).re *= green * smoothing;
		  (cvector_fft[ThisGrid][index/2]).im *= green * smoothing;
		}

              if (swap)
		{
		  tmp                                 = (cvector_fft[ThisGrid][index/2]).im;
		  (cvector_fft[ThisGrid][index/2]).im = (cvector_fft[ThisGrid][index/2]).re;
		  (cvector_fft[ThisGrid][index/2]).re = -tmp;
		}
	    }
	}
    }

  if (!ThisTask)
    printf("[%s] compute_derivative: starting fft\n",fdate());

  time=reverse_transform(ThisGrid);

  if (!ThisTask)
    printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);

  cputime.fft+=time;

  return 0;
}


double greens_function(double *diff_comp, double k_squared, int first_derivative, int second_derivative)
{

  /* in this case the greens_function is simply 1 */
  if (first_derivative == -1 && second_derivative == -1)
    return 1.0;

  if (first_derivative==0 && second_derivative==0)
    return -diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;
  else
    return diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;

}


void write_in_cvector(int ThisGrid, double *vector)
{
  for (int i = 0; i < MyGrids[ThisGrid].total_local_size_fft; i++)
    *((double*)cvector_fft[ThisGrid]+i)=*(vector+i);

}


void write_from_cvector(int ThisGrid, double *vector)
{
  for (int i = 0; i < MyGrids[ThisGrid].total_local_size_fft; i++)
    *(vector+i)=*((double*)cvector_fft[ThisGrid]+i);

}


void write_in_rvector(int ThisGrid, double *vector)
{
  for (int lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (int ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (int lx=0; lx<MyGrids[ThisGrid].GSlocal_x; lx++)
	*(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y)) =
	    *(vector + lx + MyGrids[ThisGrid].GSlocal_x *(ly + lz* MyGrids[ThisGrid].GSlocal_y));
  for (int lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (int ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (int lx=MyGrids[ThisGrid].GSlocal_x; lx<MyGrids[ThisGrid].GSlocal_x+MyGrids[ThisGrid].off; lx++)
	*(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y)) = 0.0;

}


void write_from_rvector(int ThisGrid, double *vector)
{
  for (int lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (int ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (int lx=0; lx<MyGrids[ThisGrid].GSlocal_x; lx++)
	*(vector + lx + MyGrids[ThisGrid].GSlocal_x *(ly + lz* MyGrids[ThisGrid].GSlocal_y)) =
          *(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y));
}
