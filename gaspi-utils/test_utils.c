#include <GASPI.h>
#include <stdlib.h>

#include "gaspi_utils.h"

void test_bcast(gaspi_rank_t iproc, gaspi_rank_t nproc)
{
    gaspi_size_t       seg_size = sizeof(int) * 2;
    gaspi_segment_id_t seg_id = 200;
    gaspi_pointer_t    seg_ptr = NULL;

/*    if(iproc == 1)
    {
	seg_id = 201;
    }
*/  
    UITLS_CHECK_ERROR(gaspi_segment_create(seg_id, seg_size, GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_MEM_UNINITIALIZED));	

    UITLS_CHECK_ERROR(gaspi_segment_ptr(seg_id, &seg_ptr));
    if(iproc == 0)
    {
	int * seg_iptr = (int *) seg_ptr;
	seg_iptr[0] = 42;
	seg_iptr[1] = 1337;
    }
    
    UITLS_CHECK_ERROR(gaspi_bcast_asym(seg_id, 0UL, seg_size, 100, 0));

    int * seg_iptr = (int *) seg_ptr;
    gaspi_printf("[0] %d\n", seg_iptr[0]);
    gaspi_printf("[1] %d\n", seg_iptr[1]);

    UITLS_CHECK_ERROR(gaspi_segment_delete(seg_id));
}


int main(int argc, char ** argv)
{
    gaspi_rank_t iProc;
    gaspi_rank_t nProc;
    gaspi_segment_id_t seg_id = 0;
    gaspi_pointer_t seg_ptr = NULL;
    int * seg_iptr = NULL;
    const gaspi_size_t seg_size = 2 * sizeof(int);
    const gaspi_notification_id_t notify_id = 0;
    const gaspi_notification_t notify_val = 42;
    gaspi_queue_id_t queue_id = 0;
    gaspi_number_t isInitialized = GASPI_ERROR;

	
    UITLS_CHECK_ERROR(gaspi_proc_init( GASPI_BLOCK ));
    UITLS_CHECK_ERROR(gaspi_initialized(&isInitialized));
    gaspi_printf("GASPI init was called %i\n", isInitialized);
    
    UITLS_CHECK_ERROR(gaspi_proc_rank(&iProc));
    UITLS_CHECK_ERROR(gaspi_proc_num(&nProc));

    UITLS_CHECK_ERROR(create_segment(seg_size, &seg_id));

    UITLS_CHECK_ERROR(gaspi_segment_ptr(seg_id, &seg_ptr));
    seg_iptr = (int *) seg_ptr;
    seg_iptr[1] = iProc;

    gaspi_rank_t const iNext = (iProc + 1 + nProc) % nProc;

    check_queue_size(queue_id);
    wait_for_queue_entries(&queue_id, 1);

    UITLS_CHECK_ERROR(gaspi_write_notify(seg_id, sizeof(int), iNext, seg_id, 0UL, sizeof(int), notify_id, notify_val, queue_id, GASPI_BLOCK));

    gaspi_notification_id_t first_id;
    gaspi_notification_t    old_val;
    blocking_waitsome(notify_id, 1, &first_id, &old_val, seg_id);

    if(old_val != notify_val)
    {
	gaspi_printf("Get wrong notify value\n");
    }
    gaspi_printf("prev %d\n",seg_iptr[0]);

    UITLS_CHECK_ERROR(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
    test_bcast(iProc, nProc);

    
    delete_all_segments();

    UITLS_CHECK_ERROR(gaspi_proc_term(GASPI_BLOCK));

    return EXIT_SUCCESS;
}
