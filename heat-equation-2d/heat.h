#ifndef HEAT_H_
#define HEAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <GASPI.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define IDX(I,J,W) ((I)*(W)+(J))
#define TIMEOUT GASPI_BLOCK
#define NO_NEIGHBOUR        -1
#define NO_JUMP              0UL

#define NOTIFY_ID_BEGIN 1

#define NOTIFY_ID_NORTH_HALO NOTIFY_ID_BEGIN
#define NOTIFY_ID_EAST_HALO  ((NOTIFY_ID_BEGIN) + 1)
#define NOTIFY_ID_SOUTH_HALO ((NOTIFY_ID_BEGIN) + 2)
#define NOTIFY_ID_WEST_HALO  ((NOTIFY_ID_BEGIN) + 3)

#define NOTIFY_ID_END NOTIFY_ID_WEST_HALO

#define NOTIFY_VAL_NORTH_HALO 41
#define NOTIFY_VAL_EAST_HALO  42
#define NOTIFY_VAL_SOUTH_HALO 43
#define NOTIFY_VAL_WEST_HALO  44

#define SEG_ID(seg_type) ((seg_type == SEG_U) ? uid : vid)

typedef int64_t neighbour_rank_t;

typedef struct halo_offset
{
	gaspi_offset_t halo_offset;
	gaspi_offset_t jump_offset;
}halo_t;

typedef struct neighbour{
	neighbour_rank_t        rank_no;
	gaspi_notification_id_t send_id;
	gaspi_notification_t    send_value;
	gaspi_notification_id_t recv_id;
	gaspi_notification_t    recv_value;
	gaspi_offset_t			data_offset;
	gaspi_offset_t			jump_offset;
	gaspi_size_t            element_count;
}neighbour_t;

typedef enum seg_type
{
	SEG_U,
	SEG_V
}seg_type_t;

extern const gaspi_segment_id_t uid;
extern const gaspi_segment_id_t vid;
extern const gaspi_segment_id_t aid;
extern gaspi_segment_id_t hid_u;
extern gaspi_segment_id_t hid_v;

extern double * u_double_ptr;
extern double * v_double_ptr;
extern double * a_double_ptr;
extern double * hu_double_ptr;
extern double * hv_double_ptr;

extern uint64_t m;
extern uint64_t n;
extern const int64_t rank_col_count;
extern const int64_t rank_row_count;

extern size_t typesize;

extern gaspi_rank_t rank;
extern gaspi_rank_t rankcount;
extern uint64_t     neighbour_count;

extern halo_t north_halo;
extern halo_t east_halo;
extern halo_t south_halo;
extern halo_t west_halo;

extern neighbour_t north_rank;
extern neighbour_t east_rank;
extern neighbour_t south_rank;
extern neighbour_t west_rank;

void init_map(uint64_t m, uint64_t n, size_t typesize );
void check_dma_requests( gaspi_queue_id_t queue );
void print_map(FILE * out, seg_type_t seg);
FILE * get_file_handle(void);
void init_halos(uint64_t m, uint64_t n, size_t typesize);
void init_neightbours(uint64_t m, uint64_t n, size_t typesize);
gaspi_return_t read_halo(neighbour_t neighbour, halo_t halo_offset, seg_type_t seg_type, size_t typesize);
gaspi_return_t read_all_halos(seg_type_t type, size_t type_size);

gaspi_return_t notify_neighbour(neighbour_t n, seg_type_t seg_type, gaspi_queue_id_t qid);
gaspi_return_t notify_all_neighbours(seg_type_t seg_type, gaspi_queue_id_t qid);


bool
check_notify(neighbour_t n,
			 gaspi_notification_id_t recv_id,
			 gaspi_notification_t recv_val);

gaspi_return_t blocking_waitsome(  gaspi_notification_id_t   id_begin,
				  				   gaspi_notification_id_t   id_count,
                                   gaspi_notification_id_t*  id_available,
								   gaspi_notification_t*     notify_val,
						           gaspi_segment_id_t        seg);
void check_errors( gaspi_return_t retval, char * function, char * msg );
void finish_map(void);

#endif /* HEAT_H_ */
