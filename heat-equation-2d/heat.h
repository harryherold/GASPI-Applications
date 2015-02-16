#ifndef HEAT_H_
#define HEAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <GASPI.h>

#define IDX(I,J,W) ((I)*(W)+(J))
#define TIMEOUT 2000
#define NO_NEIGHBOUR        -1
#define NO_JUMP              0UL
#define NOTIFY_ID_NORTH_HALO 1
#define NOTIFY_ID_EAST_HALO  2
#define NOTIFY_ID_SOUTH_HALO 3
#define NOTIFY_ID_WEST_HALO  4

#define NOTIFY_VAL_NORTH_HALO 41
#define NOTIFY_VAL_EAST_HALO  42
#define NOTIFY_VAL_SOUTH_HALO 43
#define NOTIFY_VAL_WEST_HALO  44


typedef int64_t neighbour_rank_t;

typedef struct halo_offset
{
	gaspi_offset_t halo_offset;
	gaspi_offset_t jump_offset;
}halo_t;

typedef struct neighbour{
	neighbour_rank_t        rank_no;
	gaspi_notification_id_t notify_id;
	gaspi_notification_t    notify_value;
	gaspi_offset_t			data_offset;
	gaspi_offset_t			jump_offset;
	gaspi_size_t            element_count;
}neighbour_t;

typedef enum seg_type
{
	SEG_U,
	SEG_V
}seg_type_t;

extern gaspi_segment_id_t uid;
extern gaspi_segment_id_t vid;
extern gaspi_segment_id_t aid;
extern gaspi_segment_id_t hid_u;
extern gaspi_segment_id_t hid_v;

extern double heat_u_min;
extern double heat_u_max;
extern double heat_alpha;

extern double * u_double_ptr;
extern double * v_double_ptr;
extern double * a_double_ptr;
extern double * hu_double_ptr;
extern double * hv_double_ptr;

extern uint64_t m;
extern uint64_t n;
extern int64_t rank_col_count;
extern int64_t rank_row_count;

extern size_t typesize;

extern gaspi_rank_t rank;
extern gaspi_rank_t rankcount;

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
void check_errors( gaspi_return_t retval, char * function, char * msg );
void finish_map(void);

#endif /* HEAT_H_ */
