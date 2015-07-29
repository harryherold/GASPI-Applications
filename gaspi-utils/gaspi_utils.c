#include "gaspi_utils.h"
#include <string.h>
static gaspi_segment_id_t gaspi_utils_seg_counter;

void
delete_all_segments()
{
    while(gaspi_utils_seg_counter-- > 0)
    {
	CHECK_GASPI_ERRORS(gaspi_segment_delete(gaspi_utils_seg_counter));
    }
}

gaspi_return_t
create_segment(const gaspi_size_t size, gaspi_segment_id_t * seg_id)
{
    gaspi_number_t seg_max = 0;
    CHECK_GASPI_ERRORS(gaspi_segment_max(&seg_max));

    if(seg_max < gaspi_utils_seg_counter)
    {
	gaspi_printf("Error: Can't create segment, reached max. of segments\n");
	return GASPI_ERROR;
    }
    
    CHECK_GASPI_ERRORS(gaspi_segment_create(gaspi_utils_seg_counter,
					    size,
					    GASPI_GROUP_ALL,
					    GASPI_BLOCK,
					    GASPI_MEM_UNINITIALIZED));
    *seg_id = gaspi_utils_seg_counter; 
    
    gaspi_utils_seg_counter++;

    return GASPI_SUCCESS;
}

void
check_queue_size( gaspi_queue_id_t queue )
{
    gaspi_number_t queue_size = 0;
    CHECK_GASPI_ERRORS(gaspi_queue_size( queue , &queue_size));

    gaspi_number_t queue_size_max = 0;
    CHECK_GASPI_ERRORS(gaspi_queue_size_max ( &queue_size_max ));

    if(  queue_size >= queue_size_max )
    {
	CHECK_GASPI_ERRORS(gaspi_wait( queue , GASPI_BLOCK ));
    }
}

void
wait_for_queue_entries ( gaspi_queue_id_t* queue, gaspi_number_t wanted_entries )
{
    gaspi_number_t queue_size_max;
    gaspi_number_t queue_size;
    gaspi_number_t queue_num;

    CHECK_GASPI_ERRORS(gaspi_queue_size_max( &queue_size_max));
    CHECK_GASPI_ERRORS(gaspi_queue_size(*queue, &queue_size));
    CHECK_GASPI_ERRORS(gaspi_queue_num(&queue_num));

    if (! (queue_size + wanted_entries <= queue_size_max))
    {
	*queue = (*queue + 1) % queue_num;
	CHECK_GASPI_ERRORS(gaspi_wait(*queue, GASPI_BLOCK));
    }
}

void
blocking_waitsome(gaspi_notification_id_t id_begin, gaspi_notification_id_t id_count,
                  gaspi_notification_id_t* id_available, gaspi_notification_t* notify_val, gaspi_segment_id_t seg)
{
    CHECK_GASPI_ERRORS(gaspi_notify_waitsome(seg, id_begin, id_count, id_available, GASPI_BLOCK));
    CHECK_GASPI_ERRORS(gaspi_notify_reset(seg, *id_available, notify_val));
}

void
flush_queues(gaspi_queue_id_t queue_begin, gaspi_queue_id_t queue_count)
{
    gaspi_number_t queue_size = 0;

    for(gaspi_queue_id_t queue = queue_begin; queue < (queue_count + queue_begin) ; ++queue)
    {
        CHECK_GASPI_ERRORS(gaspi_queue_size(queue, &queue_size));

        if(queue_size > 0)
        {
            CHECK_GASPI_ERRORS(gaspi_wait(queue, GASPI_BLOCK));
        }
    }
}

/************************************************************************************
 *
 * collective operations helper functions
 *
 ************************************************************************************/
static unsigned int
npot( unsigned int v )
{
    v--;
    v |= v >>  1;
    v |= v >>  2;
    v |= v >>  4;
    v |= v >>  8;
    v |= v >> 16;
    
    return v + 1;
}

unsigned int
gaspi_utils_get_bino_num( int rank, int root, int rank_count )
{
    rank -= root;
    if ( rank < 0 )
    {
        rank += rank_count;
    }
    return rank;
}


unsigned int
gaspi_utils_get_rank( int relative_rank, int root, int rank_count )
{
    relative_rank += root;
    if ( relative_rank >= rank_count )
    {
        relative_rank -= rank_count;
    }
    return relative_rank;
}

int
gaspi_utils_compute_comms( int*  parent,
                           int** children,
                           int   me,
                           int   root )
{
    gaspi_rank_t size;

    CHECK_GASPI_ERRORS(gaspi_proc_num( &size ));
    unsigned int size_pot = npot( size );

    unsigned int d;
    unsigned int number_of_children = 0;

    /* Be your own parent, ie. the root, by default */
    *parent = me;

    me = gaspi_utils_get_bino_num( me, root, size );

    /* Calculate the number of children for me */
    for ( d = 1; d; d <<= 1 )
    {
        /* Actually break condition */
        if ( d > size_pot )
        {
            break;
        }

        /* Check if we are actually a child of someone */
        if ( me & d )
        {
            /* Yes, set the parent to our real one, and stop */
            *parent = me ^ d;
            *parent = gaspi_utils_get_rank( *parent, root, size );
            break;
        }

        /* Only count real children, of the virtual hypercube */
        if ( ( me ^ d ) < size )
        {
            number_of_children++;
        }
    }

    /* Put the ranks of all children into a list and return */
    *children = malloc( sizeof( **children ) * number_of_children );
    unsigned int child = number_of_children;

    d >>= 1;
    while ( d )
    {
        if ( ( me ^ d ) < size )
        {
            ( *children )[ --child ] = me ^ d;
            ( *children )[ child ]   = gaspi_utils_get_rank( ( *children )[ child ], root, size );
        }
        d >>= 1;
    }

    return number_of_children;
}
/************************************************************************************
 *
 * collective operations helper functions
 *
 ************************************************************************************/
/**
 * blocking operation for GASPI_GROUP_ALL
 * @param seg_id:   same segement for all participating processes
 * @param offset:   same offset for all participating processes
 * @param root:     Rank number of the root processes
 * @param bytesize: Size in bytes of the data 
 */
gaspi_return_t
gaspi_bcast( gaspi_segment_id_t    seg_id,
	     gaspi_offset_t        offset,
	     gaspi_size_t          bytesize,
             gaspi_rank_t          root )
{
    int                     children_count;
    int                     child;
    gaspi_rank_t            rank;
    gaspi_rank_t            rankcount;
    gaspi_return_t          retval    = GASPI_SUCCESS;
    gaspi_pointer_t         p_segment = NULL;
    gaspi_notification_id_t notify_id = 0;
    gaspi_queue_id_t        queue     = 0;
    gaspi_notification_id_t first_id;
    gaspi_notification_t    old_value;

    int    parent;
    int*   children = NULL;

    CHECK_GASPI_ERRORS(gaspi_proc_rank( &rank ));
    CHECK_GASPI_ERRORS(gaspi_proc_num( &rankcount ));

    CHECK_GASPI_ERRORS(gaspi_segment_ptr( seg_id, &p_segment ));
    p_segment = p_segment + offset;
    
    children_count = gaspi_utils_compute_comms( &parent, &children, rank, root );

    CHECK_GASPI_ERRORS(gaspi_barrier( GASPI_GROUP_ALL, GASPI_BLOCK ));
    /*
     * parents + children wait for upper parents data
     */
    if ( rank != parent )
    {
	blocking_waitsome(notify_id, 1, &first_id, &old_value, seg_id);	
    }
    /*
     * write to all childs
     */
    for ( child = 0; child < children_count; child++ )
    {
	check_queue_size(queue);
        CHECK_GASPI_ERRORS(gaspi_write_notify(
            seg_id,
            offset,
            children[ child ],
            seg_id,
            offset,
            bytesize,
            notify_id,
            42,
            queue,
            GASPI_BLOCK ));
    }

    free( children );
    CHECK_GASPI_ERRORS(gaspi_barrier( GASPI_GROUP_ALL, GASPI_BLOCK ));
    return retval;
}

gaspi_return_t
gaspi_bcast_asym( gaspi_segment_id_t    seg_id,
		  gaspi_offset_t        offset,
		  gaspi_size_t          bytesize,
		  gaspi_segment_id_t    transfer_seg_id,
		  gaspi_rank_t          root)

{
    gaspi_rank_t iProc;
    gaspi_pointer_t transfer_ptr = NULL;
    gaspi_pointer_t src_ptr = NULL;
    CHECK_GASPI_ERRORS(gaspi_proc_rank(&iProc));
    CHECK_GASPI_ERRORS(gaspi_segment_create(transfer_seg_id,
					    bytesize,
					    GASPI_GROUP_ALL,
					    GASPI_BLOCK,
					    GASPI_MEM_UNINITIALIZED));

    CHECK_GASPI_ERRORS(gaspi_segment_ptr(transfer_seg_id, &transfer_ptr));
    CHECK_GASPI_ERRORS(gaspi_segment_ptr(seg_id, &src_ptr));

    if(iProc == root)
    {
	memcpy(transfer_ptr, src_ptr + offset, bytesize);
    }

    CHECK_GASPI_ERRORS(gaspi_bcast(transfer_seg_id, 0UL, bytesize, root));

    if(iProc != root)
    {
	memcpy(src_ptr + offset, transfer_ptr, bytesize);
    }
    
    CHECK_GASPI_ERRORS(gaspi_segment_delete(transfer_seg_id));

    return GASPI_SUCCESS;
}
