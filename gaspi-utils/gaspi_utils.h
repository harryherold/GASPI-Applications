#ifndef GASPI_UTILS_H
#define GASPI_UTILS_H

#include <stdlib.h>
#include <GASPI.h>

gaspi_return_t gaspi_bcast(gaspi_segment_id_t seg_id, gaspi_offset_t offset,
                           gaspi_size_t bytesize, gaspi_rank_t root);

gaspi_return_t gaspi_bcast_asym(gaspi_segment_id_t seg_id,
                                gaspi_offset_t offset, gaspi_size_t bytesize,
                                gaspi_segment_id_t transfer_seg_id,
                                gaspi_rank_t root);

gaspi_return_t create_segment(const gaspi_size_t size,
                              gaspi_segment_id_t *seg_id);

void delete_all_segments();

void check_queue_size(gaspi_queue_id_t queue);

void wait_for_queue_entries(gaspi_queue_id_t *queue,
                            gaspi_number_t wanted_entries);

void blocking_waitsome(gaspi_notification_id_t id_begin,
                       gaspi_notification_id_t id_count,
                       gaspi_notification_id_t *id_available,
                       gaspi_notification_t *notify_val,
                       gaspi_segment_id_t seg);

void flush_queues(gaspi_queue_id_t queue_begin, gaspi_queue_id_t queue_count);

#define UITLS_CHECK_ERROR(func...)                                            \
  do {                                                                         \
    const gaspi_return_t retval = func;                                        \
    if (retval != GASPI_SUCCESS) {                                             \
      gaspi_printf("ERROR in %s : %s on line %i return value %i\n", #func,     \
                   __FILE__, __LINE__, retval);                                \
    }                                                                          \
  }while (0)

#define UITLS_CHECK_ERROR_RET(ret, func...)                                         \
    do {                                                                         \
      ret = func;                                                                \
      if (retval != GASPI_SUCCESS) {                                             \
        gaspi_printf("ERROR in %s : %s on line %i return value %i\n", #func,     \
                     __FILE__, __LINE__, retval);                                \
      }                                                                          \
    }while (0)

typedef struct reportentry{
  unsigned int id;
  char ret;
  int line_num;
}report_entry_t;

extern unsigned long long test_case_count;
extern gaspi_pointer_t test_seg_ptr;
extern report_entry_t * report_ptr;

void report_to_master(unsigned long long idx);

#define TEST_EQUAL(id, expected, actual)           \
  do{                                              \
    const unsigned int idx = __COUNTER__;          \
    report_ptr[idx].id       = id;                 \
    report_ptr[idx].line_num = __LINE__;           \
    report_ptr[idx].ret      = expected == actual; \
    report_to_master(idx);                         \
  }while(0)

void init_gaspi_test(gaspi_segment_id_t seg_num, gaspi_rank_t master_rank);
void finish_gaspi_test(void);



#endif // GASPI_UTILS_H
