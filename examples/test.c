#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>

#include <gsl/gsl_rng.h>

#include <mpi.h>

sem_t wait_count;
sem_t finish;

int flag = 0;
int total;

/* Prints x’s to stderr.  The parameter is unused.  Does not return.
*/
void* print_xs (void* unused)
{
    int i=0;
    int number;
    int rank, numtasks;
    
    while(1)
    {
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	
	if (rank == 0)
	    sem_wait (&wait_count);

	MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (flag == 1)
	    break;

	number = 0;
	i = 0;
	while (i<rank)
	{
	    number++;
	    printf("%d", number);
	    fflush(stdout);
	    i++;
	}
	
	fflush(stderr);
	MPI_Reduce(&number, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
	    sem_post(&finish);
	
    }
    return NULL;
}
/* The main program.
*/
int main ()
{
    gsl_rng *r;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    sem_init (&wait_count, 0, 0);
    sem_init (&finish, 0, 0);

    pthread_t thread_id;
    /* Create a new thread. The new thread will run the print_xs
       function. */
    pthread_create (&thread_id, NULL, &print_xs, NULL);

    if (rank==0)
    {
	/* Print o’s continuously to stderr. */
	int i = 0;
	while (i<10)
	{
	    fputc ('o', stderr);
	    fflush(stderr);
	    sem_post(&wait_count);
	    sem_wait(&finish);
	    printf("t: %d\n", total);
	    fflush(stdout);
	    i++;
	}
	flag = 1;
	sem_post(&wait_count);
    }
    
    pthread_join(thread_id, NULL);
    fputc ('\n', stderr);
    MPI_Finalize();
    return 0;


}
