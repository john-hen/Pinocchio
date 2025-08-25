/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan, 
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025
 
 github: https://github.com/pigimonaco/Pinocchio
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
#include "def_splines.h"
#include "initial_density.h"

void abort_code(void);
void write_cputimes(void);

int main(int argc, char **argv, char **envp)
{

  /* Initialize MPI */
  int got_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_level); // Hybrid MPI and OPENMP parallel
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); // gives you your rank i.e. ID
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks); // size of your pool 

#ifdef _OPENMP
  /* initialization of OpemMP */
#pragma omp parallel
  {
#pragma omp master
    internal.nthreads_omp = omp_get_num_threads();
  }
#endif

  /* timing of the code */
  cputime.total=MPI_Wtime();
  greetings();

  /* checks that the parameter file is given in the command line */
  if (argc<2)
    {
      if (!ThisTask)
	printf("Usage: pinocchio.x parameterfile\n");
      MPI_Finalize();
      return 0;
    }

  /* exit now if a snapshot is required and SNAPSHOT is not set */
#ifndef SNAPSHOT
  if (argc>=3 && atoi(argv[2])>1)
    {
      if (!ThisTask)
	printf("Sorry but you have to use the SNAPSHOT directive in compilation to write a snapshot\n");
      MPI_Finalize();
      return 0;
    }
#endif

  /* exit now if a snapshot is required and SNAPSHOT is not set */
#ifndef TABULATED_CT
  if (argc>=3 && atoi(argv[2])==1)
    {
      if (!ThisTask)
	printf("Sorry but you have to use the TABULATED_CT directive in compilation to write the collapse time table\n");
      MPI_Finalize();
      return 0;
    }
#endif


  /* initialization */
  memset(&params, 0, sizeof(param_data));
  strcpy(params.ParameterFile,argv[1]);
  if (initialization())
    abort_code();

  /*****************************************/
  /*********** Special behaviour ***********/
  /*****************************************/
  /* called as "pinocchio.x parameterfile 1" it computes and writes collapse time table, then exit */
  if (argc>=3 && atoi(argv[2])==1)
    {
#ifdef TABULATED_CT
      if (!ThisTask)
	{
	  printf("In this configuration pinocchio will only compute a table of collapse times\n");
	}

      /* CYCLE ON SMOOTHING RADII */
      for (int ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
	{
	  double cputmp=MPI_Wtime();

	  if (!ThisTask)
	    printf("\n[%s] Starting smoothing radius %d of %d (R=%9.5f, sigma=%9.5f)\n", 
		   fdate(), ismooth+1, Smoothing.Nsmooth, Smoothing.Radius[ismooth],
		   sqrt(Smoothing.Variance[ismooth]) );

	  if (initialize_collapse_times(ismooth,1))
	    return 1;

	  if (!ThisTask)
	    printf("[%s] Collapse times computed, cpu time =%f s\n",fdate(),cputmp);

	}

      if (!ThisTask)
	printf("Pinocchio done!\n");

      MPI_Finalize();

#endif
      return 0;
    }

  /* On request, it writes the density field */
  if (argc>=3 && atoi(argv[2])==2)
    {
#ifdef SNAPSHOT
      /* called as "pinocchio.x parameterfile 2" it writes the density field
	 in configuration space and exits */
      if (!ThisTask)
	printf("In this configuration pinocchio only writes the linear density field in a snapshot\n");

      for (int ThisGrid=0; ThisGrid<Ngrids; ThisGrid++)
	{
	  write_in_cvector(ThisGrid, kdensity[ThisGrid]);
	  double time=reverse_transform(ThisGrid);
	  if (!ThisTask)
	    printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);
	  write_from_rvector(ThisGrid, density[ThisGrid]);
	  if (write_density(ThisGrid))
	    abort_code();
	}
#else
      if (!ThisTask)
	printf("Please compile the code with SNAPSHOT directive to use this option\n");
      
#endif
      if (!ThisTask)
	{
	  write_cputimes();
	  printf("Pinocchio done!\n");
	}

      MPI_Finalize();
      return 0;
    }

  /* called as "pinocchio.x parameterfile 3" it computes displacements and writes a standard snapshot, 
     then exit */
  if (argc>=3 && atoi(argv[2])==3)
    {
#ifdef SNAPSHOT

      if (!ThisTask)
	      {
	        printf("In this configuration pinocchio will only produce a GADGET snapshot\n");
	        printf("at the first redshift specified by the %s file (z=%f)\n", params.OutputList,outputs.z[0]);
        }

      /* compute displacements for all particles */
      double cputmp=MPI_Wtime();
      if (!ThisTask)
	printf("\n[%s] Computing displacements\n",fdate());

      if (compute_displacements(1,1,outputs.z[0]))
          abort_code();

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done computing displacements, cpu time = %f s\n",fdate(),cputmp);

      /* write the snapshot */
      cputmp=MPI_Wtime();
      if (!ThisTask)
	printf("\n[%s] Writing the snapshot\n",fdate());

      if (write_LPT_snapshot())
          abort_code();

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done snapshot, cpu time = %f s\n",fdate(),cputmp);
#endif

      if (!ThisTask) 
        printf("Pinocchio done!\n");

      MPI_Finalize();

      return 0; 
    } 

  // Output initial density, then exit, when called as:
  //   pinocchio parameter_file --write-initial-density
  // The output file name can be specified as an optional third command-line
  // argument that defaults to "initial_density.dat".
  if (argc >= 3 && strcmp(argv[2], "--write-initial-density") == 0) {
    char* const file_name = (argc >= 4) ? argv[3] : "initial_density.dat";
    if (ThisTask == 0)
      printf("Writing initial density to file \"%s\".\n", file_name);
    int ThisGrid = 0;         // Not implemented for multiple grids.
    int exit_code = write_initial_density(
      kdensity, file_name, ThisGrid, MyGrids
    );
    MPI_Finalize();
    return exit_code;
  }

  // Abort on any other second command-line argument.
  if (argc == 3) {
    printf("Error: Unknown command-line argument \"%s\".\n", argv[2]);
    MPI_Finalize();
    return 1;
  }

  // Abort if more arguments are passed.
  if (argc > 3) {
    printf("Error: Too many command-line argument.\n");
    MPI_Finalize();
    return 1;
  }


  /******************************************/
  /*********** Standard behaviour ***********/
  /******************************************/


  if (params.ReadProductsFromDumps)
    {
      /* on request, read products from dump files and skip the first part */
      if (read_dumps())
	abort_code();
    }
  else
    {
      /* computation of collapse times and displacements */
      if (compute_fmax())
	abort_code();

      /* on request, dump products to files for skipping fmax */
      if (params.DumpProducts)
	if (dump_products())
	  abort_code();
    }

  /* fragmentation of the collapsed medium */
  if (fragment_driver())
    abort_code();

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* output detailed cpu times */
  cputime.total=MPI_Wtime()-cputime.total;
  if (!ThisTask)
    write_cputimes();

  /* done */
  if (!ThisTask)
    printf("Pinocchio done!\n");
  MPI_Finalize();

  return 0;
}


void abort_code(void)
{
  printf("Task %d aborting...\n",ThisTask);
  MPI_Abort(MPI_COMM_WORLD,1);
}


void write_cputimes()
{
  printf("Total:            %14.6f\n", cputime.total);
  printf("Initialization:   %14.6f (%5.2f%%)\n", cputime.init, 100.*cputime.init/cputime.total);
  printf("  Density in PS:  %14.6f (%5.2f%%)\n", cputime.dens, 100.*cputime.dens/cputime.total);
  printf("fmax:             %14.6f (%5.2f%%)\n", cputime.fmax, 100.*cputime.fmax /cputime.total);
#ifdef TWO_LPT
  printf("  LPT:            %14.6f (%5.2f%%)\n", cputime.lpt,  100.*cputime.lpt  /cputime.total);
#endif
  printf("  Derivatives:    %14.6f (%5.2f%%)\n", cputime.deriv,  100.*cputime.deriv  /cputime.total);
  printf("    Mem transfer: %14.6f (%5.2f%%)\n", cputime.mem_transf, 100.*cputime.mem_transf  /cputime.total);
  printf("    FFTs:         %14.6f (%5.2f%%)\n", cputime.fft,  100.*cputime.fft  /cputime.total);
  printf("  Collapse times: %14.6f (%5.2f%%)\n", cputime.coll, 100.*cputime.coll /cputime.total);
  printf("    inv.collapse: %14.6f (%5.2f%%)\n", cputime.invcoll, 100.*cputime.invcoll /cputime.total);
  printf("    ellipsoid:    %14.6f (%5.2f%%)\n", cputime.ell, 100.*cputime.ell /cputime.total);
  printf("  Velocities:     %14.6f (%5.2f%%)\n", cputime.vel,  100.*cputime.vel  /cputime.total);
  printf("Fragmentation:    %14.6f (%5.2f%%)\n", cputime.frag, 100.*cputime.frag /cputime.total);
  printf("  Redistribution: %14.6f (%5.2f%%)\n", cputime.distr,100.*cputime.distr/cputime.total);
  printf("  Sorting:        %14.6f (%5.2f%%)\n", cputime.sort, 100.*cputime.sort /cputime.total);
#ifdef PLC
  printf("  Groups total:   %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
  printf("  Groups PLC:     %14.6f (%5.2f%%)\n", cputime.plc,100.*cputime.plc/cputime.total);
#else
  printf("  Groups:         %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
#endif
  printf("Total I/O:        %14.6f (%5.2f%%)\n", cputime.io,   100.*cputime.io   /cputime.total);
}

