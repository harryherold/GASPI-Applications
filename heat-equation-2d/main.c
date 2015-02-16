
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <GASPI.h>

#define IDX(I,J,W) ((I)*(W)+(J))
#define TIMEOUT 2000
#define NO_NEIGHBOUR        -1
#define NO_JUMP              0
#define NOTIFY_ID_NORTH_HALO 1
#define NOTIFY_ID_EAST_HALO  2
#define NOTIFY_ID_SOUTH_HALO 3
#define NOTIFY_ID_WEST_HALO  4

#define NOTIFY_VAL_NORTH_HALO 41
#define NOTIFY_VAL_EAST_HALO  42
#define NOTIFY_VAL_SOUTH_HALO 43
#define NOTIFY_VAL_WEST_HALO  44


typedef int64_t neighbour_rank_t;

gaspi_segment_id_t uid   = 0;
gaspi_segment_id_t vid   = 1;
gaspi_segment_id_t aid   = 2;
gaspi_segment_id_t hid_u = 3;
gaspi_segment_id_t hid_v = 4;

double heat_u_min= 10.0;
double heat_u_max= 100.0;
double heat_alpha= 22.0;

double * u_double_ptr = NULL;
double * v_double_ptr = NULL;
double * a_double_ptr = NULL;
double * hu_double_ptr = NULL;
double * hv_double_ptr = NULL;

uint64_t m              = 0;
uint64_t n              = 0;
uint64_t rank_col_count = 2;
uint64_t rank_row_count = 2;

size_t typesize         = 0;

gaspi_rank_t rank;
gaspi_rank_t rankcount;

typedef struct halo_offset
{
	gaspi_offset_t halo_offset;
	gaspi_offset_t jump_offset;
}halo_offset_t;

typedef struct neighbour{
	neighbour_rank_t        rank_no;
	gaspi_notification_id_t notify_id;
	gaspi_notification_t    notify_value;
	gaspi_offset_t			data_offset;
	gaspi_offset_t			jump_offset;
}neighbour_t;

halo_offset_t north_halo;
halo_offset_t east_halo;
halo_offset_t south_halo;
halo_offset_t west_halo;

neighbour_t north_rank;
neighbour_t east_rank;
neighbour_t south_rank;
neighbour_t west_rank;

typedef enum seg_type
{
	SEG_U,
	SEG_V
}seg_type_t;

void
check_errors( gaspi_return_t retval, char * function, char * msg )
{
  if( retval == GASPI_TIMEOUT ) {
    gaspi_printf("[[ function : %s : reached timeout: %d ms ]]\n",function,TIMEOUT);
  }
  else if( retval == GASPI_ERROR ) {
    gaspi_printf("[[ function : %s : get an GASPI-ERROR ]]\n",function);
  }
  if( msg != NULL && retval != GASPI_SUCCESS )
    gaspi_printf("%s\n", msg );
}

void
check_dma_requests( gaspi_queue_id_t queue )
{
  gaspi_number_t queue_size = 0;
  gaspi_queue_size( queue , &queue_size);

  gaspi_number_t queue_size_max = 0;
  gaspi_queue_size_max ( &queue_size_max );

  if(  queue_size >= queue_size_max )
  {
    if( GASPI_ERROR == gaspi_wait( queue , GASPI_BLOCK ) )
    {
      gaspi_printf("wait failed\n");
    }
  }
}

void
init_map(uint64_t m, uint64_t n, size_t typesize )
{
	gaspi_return_t ret = GASPI_SUCCESS;
	gaspi_pointer_t u_ptr  = NULL;
	gaspi_pointer_t v_ptr  = NULL;
	gaspi_pointer_t a_ptr  = NULL;

	ret = gaspi_segment_create(uid, m * n * typesize, GASPI_GROUP_ALL,
			GASPI_BLOCK, GASPI_MEM_INITIALIZED);

	check_errors(ret, "gaspi_segement_create", "uid");

	ret = gaspi_segment_create(vid, m * n * typesize, GASPI_GROUP_ALL,
			GASPI_BLOCK, GASPI_MEM_INITIALIZED);

	check_errors(ret, "gaspi_segement_create", "vid");

	ret = gaspi_segment_create(aid, (m - 1) * (n - 1) * typesize,
			GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_MEM_INITIALIZED);

	check_errors(ret, "gaspi_segement_create", "aid");

	gaspi_segment_ptr(uid, &u_ptr);
	gaspi_segment_ptr(vid, &v_ptr);
	gaspi_segment_ptr(aid, &a_ptr);

	u_double_ptr = (double *) u_ptr;
	v_double_ptr = (double *) v_ptr;

	for(size_t i = 0 ; i < n ; ++i)
	{
		for(size_t j = 0 ; j < m ; ++j)
		{
			if(i < ((m - 1) / 2))
			{
				u_double_ptr[IDX(i,j,m)] = heat_u_max;
				v_double_ptr[IDX(i,j,m)] = heat_u_max;
			}
			else
			{
				u_double_ptr[IDX(i,j,m)] = heat_u_min;
				u_double_ptr[IDX(i,j,m)] = heat_u_min;
			}
		}
	}

	a_double_ptr = (double *) a_ptr;

	for(size_t i = 0 ; i < n - 1 ; ++i)
	{
		for(size_t j = 0 ; j < m - 1 ; ++j)
		{
			a_double_ptr[IDX(i,j,m - 1)] = heat_alpha;
		}
	}
}

void
init_halos(uint64_t m, uint64_t n, size_t typesize)
{
	north_halo.halo_offset = typesize;
	north_halo.jump_offset = NO_JUMP;

	south_halo.halo_offset = (((n - 1) * m) + 1) * typesize;
	south_halo.jump_offset = NO_JUMP;

	west_halo.halo_offset = m * typesize;
	west_halo.jump_offset = m * typesize;

	east_halo.halo_offset = ((2 * m) - 1) * typesize;
	east_halo.jump_offset = m * typesize;
}

void
init_neightbours(uint64_t m, uint64_t n, size_t typesize)
{
	int64_t col_idx = rank % rank_col_count;
	int64_t row_idx = rank / rank_col_count;

	if(col_idx + 1 < rank_col_count)
	{
		east_rank.rank_no       = rank + 1;
		east_rank.notify_id     = NOTIFY_ID_WEST_HALO;
		east_rank.notify_value  = NOTIFY_VAL_WEST_HALO;
		east_rank.data_offset   = (m + 1) * typesize;
		east_rank.jump_offset   = m * typesize;
		east_rank.element_count = n - 2;
		gaspi_printf("east %d\n",east_rank.rank_no);
	}
	else
		east_rank.rank_no       = NO_NEIGHBOUR;

	if(col_idx - 1 >= 0)
	{
		west_rank.rank_no       = rank - 1;
		west_rank.notify_id     = NOTIFY_ID_EAST_HALO;
		west_rank.notify_value  = NOTIFY_VAL_EAST_HALO;
		west_rank.data_offset   = ((2 * m) - 2) * typesize;
		west_rank.jump_offset   = m * typesize;
		west_rank.element_count = n - 2;
		gaspi_printf("west %d\n",west_rank.rank_no);
	}
	else
		west_rank.rank_no       = NO_NEIGHBOUR;

	if(row_idx - 1 >= 0)
	{
		north_rank.rank_no       = rank - rank_col_count;
		north_rank.notify_id     = NOTIFY_ID_SOUTH_HALO;
		north_rank.notify_value  = NOTIFY_VAL_SOUTH_HALO;
		north_rank.data_offset   = (((n - 2) * m) + 1) * typesize;
		north_rank.jump_offset   = NO_JUMP;
		north_rank.element_count = m - 2;

		gaspi_printf("north %d\n",north_rank.rank_no);
	}
	else
		north_rank.rank_no      = NO_NEIGHBOUR;

	if(row_idx + 1 < rank_row_count)
	{
		south_rank.rank_no       = rank + rank_col_count;
		south_rank.notify_id     = NOTIFY_ID_NORTH_HALO;
		south_rank.notify_value  = NOTIFY_VAL_NORTH_HALO;
		south_rank.data_offset   = (m + 1) * typesize;
		south_rank.jump_offset   = NO_JUMP;
		south_rank.element_count = m - 2;

		gaspi_printf("south %d\n",south_rank.rank_no);
	}
	else
		south_rank.rank_no      = NO_NEIGHBOUR;
}

gaspi_return_t
read_data(neighbour_t * neighbour, halo_offset_t halo_offset, seg_type_t seg_type)
{
	return GASPI_SUCCESS;
}

void
finish_map()
{
	gaspi_segment_delete(uid);
	gaspi_segment_delete(vid);
	gaspi_segment_delete(aid);
}

int
main(int argc, char ** argv)
{
  if(gaspi_proc_init( GASPI_BLOCK ) != GASPI_SUCCESS )
  {
        fprintf(stderr,"startup failed\n");
        gaspi_proc_term( GASPI_BLOCK );
        exit(EXIT_FAILURE);
  }

  typesize = sizeof(double);

  gaspi_proc_rank(&rank);
  gaspi_proc_num(&rankcount);

  m = (12 / rankcount) + 1;
  n = (12 / rankcount) + 1;

  init_neightbours(m, n, typesize);

  init_map(m, n, typesize);

  finish_map(uid, vid, aid);

  gaspi_proc_term( GASPI_BLOCK );

  return EXIT_SUCCESS;
}
