#include <stdio.h>
// #include <mpi.h>
#include "main.h"

void mpi_sample(int myid, int nodes, MPI_Status istatus) {
  int i, j, type, dest, source;

	type = 777;
	if( myid == 0 )
		for (i=1; i < nodes; i++){
			j = 2*i; dest = i;
			MPI_Send( 	&j, 
					1, 
					MPI_INT, 
					dest, 
					type, 
					MPI_COMM_WORLD);
			printf(" Node %d sent a message %d to %d\n", myid, j, dest );	
		}
	else {
                source = 0;
		MPI_Recv( 	&j, 
				1,
				MPI_INT, 
				source,
				MPI_ANY_TAG,
				MPI_COMM_WORLD, 
				&istatus);
		printf(" Node %d received a message %d from %d\n", myid, j, source );	
	}
}