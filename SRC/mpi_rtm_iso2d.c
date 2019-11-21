#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "fixed_parameters.h"
#include "rtm_iso2d.h"

int main(void) {
	
	//----------------------------------------------------------------
	// MPI Variables...
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int world_size, rank, name_len;
	float simulations_per_node;
	double start, end;
	
	int simulation_index;
	int first_simulation, last_simulation;
	//----------------------------------------------------------------
	
	
	MPI_Init(NULL, NULL);
	
	// Get the number of processes...
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	// Get the rank of the process...
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// Get the name of the processor...
	MPI_Get_processor_name(processor_name, &name_len);
	
// 	printf("\n\n.................................... STARTING 2D ISOTROPIC REVERSE TIME MIGRATION !!! ..........\n");
	printf("[Host %s] [MPI rank %d] %s \n", processor_name, rank, "STARTING 2D ISOTROPIC REVERSE TIME MIGRATION !!!");
	
	
	// Number of simulations per node...
	simulations_per_node = ( (float)numberOfSimulations ) / world_size;
	
	// Simulation ranges...
	first_simulation = (int) (ceil((simulations_per_node * rank) + 1)) ;
	last_simulation  = (int) (ceil( simulations_per_node * (rank + 1)));
	
	
	start = MPI_Wtime();
	for (simulation_index = first_simulation; simulation_index <= last_simulation; simulation_index++) {
            
		rtm_routine(simulation_index);
	
	}
	end = MPI_Wtime();
// 	
	printf("--> The process %d took %lf seconds to run.\n", rank, end - start);
	
	MPI_Finalize();
        
        return 0;
}
