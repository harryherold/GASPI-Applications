#include "heat.h"
#include "test.h"

const int64_t            rank_col_count = 2;
const int64_t            rank_row_count = 2;
const size_t             matrix_width   = 40;
const size_t             matrix_height  = 40;
const gaspi_segment_id_t uid            = 0;
const gaspi_segment_id_t vid            = 1;
const gaspi_segment_id_t aid            = 2;

double * u_double_ptr = NULL;
double * v_double_ptr = NULL;
double * a_double_ptr = NULL;

uint64_t m = 0;
uint64_t n = 0;

gaspi_rank_t rank;
gaspi_rank_t rankcount;
uint64_t neighbour_count = 0;

halo_t north_halo;
halo_t east_halo;
halo_t south_halo;
halo_t west_halo;

neighbour_t north_rank;
neighbour_t east_rank;
neighbour_t south_rank;
neighbour_t west_rank;

void
check_errors(gaspi_return_t retval, char * function, char * msg)
{
    if (retval == GASPI_TIMEOUT)
    {
        gaspi_printf("[[ function : %s : reached timeout: %d ms ]]\n", function,
        TIMEOUT);
    }
    else if (retval == GASPI_ERROR)
    {
        gaspi_printf("[[ function : %s : get an GASPI-ERROR ]]\n", function);
    }
    if (msg != NULL && retval != GASPI_SUCCESS)
        gaspi_printf("%s\n", msg);
}

void
check_dma_requests(gaspi_queue_id_t queue)
{
    gaspi_number_t queue_size = 0;
    gaspi_queue_size(queue, &queue_size);

    gaspi_number_t queue_size_max = 0;
    gaspi_queue_size_max(&queue_size_max);

    if (queue_size >= queue_size_max)
    {
        if (GASPI_ERROR == gaspi_wait(queue, GASPI_BLOCK))
        {
            gaspi_printf("wait failed\n");
        }
    }
}

void
init_map(uint64_t m, uint64_t n, size_t typesize)
{
    gaspi_return_t ret    = GASPI_SUCCESS;
    gaspi_pointer_t u_ptr = NULL;
    gaspi_pointer_t v_ptr = NULL;
    gaspi_pointer_t a_ptr = NULL;

    init_halos(m, n, typesize);
    init_neightbours(m, n, typesize);

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

    double heat_u_min = 10.0;
    double heat_u_max = 100.0;
    double heat_alpha = 22.0;

    for (size_t i = 0; i < n; ++i)
    {
        bool is_hot = (i < (m / 2));
        for (size_t j = 0; j < m; ++j)
        {
            if (is_hot)
            {
                u_double_ptr[IDX(i, j, m)] = heat_u_max;
                v_double_ptr[IDX(i, j, m)] = heat_u_max;
            }
            else
            {
                u_double_ptr[IDX(i, j, m)] = heat_u_min;
                v_double_ptr[IDX(i, j, m)] = heat_u_min;
            }
        }
    }

    a_double_ptr = (double *) a_ptr;

    for (size_t i = 0; i < n - 1; ++i)
    {
        for (size_t j = 0; j < m - 1; ++j)
        {
            a_double_ptr[IDX(i, j, m - 1)] = heat_alpha;
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

    if (col_idx + 1 < rank_col_count)
    {
        east_rank.rank_no = rank + 1;
        east_rank.send_id = NOTIFY_ID_WEST_HALO;
        east_rank.send_value = NOTIFY_VAL_WEST_HALO;
        east_rank.recv_id = NOTIFY_ID_EAST_HALO;
        east_rank.recv_value = NOTIFY_VAL_EAST_HALO;
        east_rank.data_offset = (m + 1) * typesize;
        east_rank.jump_offset = m * typesize;
        east_rank.element_count = n - 2;
        ++neighbour_count;
    }
    else
        east_rank.rank_no = NO_NEIGHBOUR;

    if (col_idx - 1 >= 0)
    {
        west_rank.rank_no = rank - 1;
        west_rank.send_id = NOTIFY_ID_EAST_HALO;
        west_rank.send_value = NOTIFY_VAL_EAST_HALO;
        west_rank.recv_id = NOTIFY_ID_WEST_HALO;
        west_rank.recv_value = NOTIFY_VAL_WEST_HALO;
        west_rank.data_offset = ((2 * m) - 2) * typesize;
        west_rank.jump_offset = m * typesize;
        west_rank.element_count = n - 2;
        ++neighbour_count;
    }
    else
        west_rank.rank_no = NO_NEIGHBOUR;

    if (row_idx - 1 >= 0)
    {
        north_rank.rank_no = rank - rank_col_count;
        north_rank.send_id = NOTIFY_ID_SOUTH_HALO;
        north_rank.send_value = NOTIFY_VAL_SOUTH_HALO;
        north_rank.recv_id = NOTIFY_ID_NORTH_HALO;
        north_rank.recv_value = NOTIFY_VAL_NORTH_HALO;
        north_rank.data_offset = (((n - 2) * m) + 1) * typesize;
        north_rank.jump_offset = NO_JUMP;
        north_rank.element_count = m - 2;
        ++neighbour_count;
    }
    else
        north_rank.rank_no = NO_NEIGHBOUR;

    if (row_idx + 1 < rank_row_count)
    {
        south_rank.rank_no = rank + rank_col_count;
        south_rank.send_id = NOTIFY_ID_NORTH_HALO;
        south_rank.send_value = NOTIFY_VAL_NORTH_HALO;
        south_rank.recv_id = NOTIFY_ID_SOUTH_HALO;
        south_rank.recv_value = NOTIFY_VAL_SOUTH_HALO;
        south_rank.data_offset = (m + 1) * typesize;
        south_rank.jump_offset = NO_JUMP;
        south_rank.element_count = m - 2;
        ++neighbour_count;
    }
    else
        south_rank.rank_no = NO_NEIGHBOUR;
}

gaspi_return_t
blocking_waitsome(gaspi_notification_id_t id_begin, gaspi_notification_id_t id_count,
                  gaspi_notification_id_t* id_available, gaspi_notification_t* notify_val, gaspi_segment_id_t seg)
{
    gaspi_return_t retval;
    gaspi_notification_id_t first_id;

    retval = gaspi_notify_waitsome(seg, id_begin, id_count, (gaspi_notification_id_t * const ) &first_id,
    TIMEOUT);
    *id_available = first_id;
    gaspi_notify_reset(seg, first_id, notify_val);

    return retval;
}

gaspi_return_t
read_halo(neighbour_t neighbour, halo_t halo_offset, seg_type_t seg_type, size_t typesize)
{
    gaspi_return_t ret = GASPI_SUCCESS;
    gaspi_segment_id_t segid = (seg_type == SEG_U) ? uid : vid;
    gaspi_queue_id_t queue_id = 0;
    gaspi_size_t i_max = (neighbour.jump_offset != NO_JUMP) ? neighbour.element_count : 1;
    gaspi_size_t read_count = (neighbour.jump_offset != NO_JUMP) ? 1UL : neighbour.element_count;

    for (gaspi_size_t i = 0; i < i_max; ++i)
    {
        check_dma_requests(queue_id);

        ret = gaspi_read(segid, halo_offset.halo_offset + (i * halo_offset.jump_offset), neighbour.rank_no, segid,
                neighbour.data_offset + (i * neighbour.jump_offset), read_count * typesize, queue_id,
                GASPI_BLOCK);

        check_errors(ret, "read_halo", "gaspi_read");
    }

    return ret;
}

bool
check_notify(neighbour_t n, gaspi_notification_id_t recv_id, gaspi_notification_t recv_val)
{
    return (n.rank_no != NO_NEIGHBOUR && n.recv_id == recv_id && n.recv_value == recv_val);
}

gaspi_return_t
read_all_halos(seg_type_t seg_type, size_t type_size)
{
    gaspi_return_t ret = GASPI_SUCCESS;
    gaspi_segment_id_t seg = SEG_ID(seg_type);
    for (uint64_t i = 0; i < neighbour_count; ++i)
    {
        gaspi_notification_id_t notify_id;
        gaspi_notification_t notify_val;

        ret = blocking_waitsome(NOTIFY_ID_BEGIN, NOTIFY_ID_END, &notify_id, &notify_val, seg);
        check_errors(ret, "read_all_halos", "blocking_waitsome");

        if (check_notify(north_rank, notify_id, notify_val))
        {
            ret = read_halo(north_rank, north_halo, seg_type, type_size);
            check_errors(ret, "read_all_halos", "read north halo");
        }
        if (check_notify(south_rank, notify_id, notify_val))
        {
            ret = read_halo(south_rank, south_halo, seg_type, type_size);
            check_errors(ret, "read_all_halos", "read south halo");
        }
        if (check_notify(east_rank, notify_id, notify_val))
        {
            ret = read_halo(east_rank, east_halo, seg_type, type_size);
            check_errors(ret, "read_all_halos", "read east halo");
        }
        if (check_notify(west_rank, notify_id, notify_val))
        {
            ret = read_halo(west_rank, west_halo, seg_type, type_size);
            check_errors(ret, "read_all_halos", "read west halo");
        }
    }
    return ret;
}

gaspi_return_t
notify_neighbour(neighbour_t n, seg_type_t seg_type, gaspi_queue_id_t qid)
{
    gaspi_segment_id_t seg_id = SEG_ID(seg_type);
    gaspi_return_t ret = GASPI_SUCCESS;
    check_dma_requests(qid);
    ret = gaspi_notify(seg_id, n.rank_no, n.send_id, n.send_value, qid,
    TIMEOUT);
    check_errors(ret, "notify_neighbour", NULL);

    return ret;
}

gaspi_return_t
notify_all_neighbours(seg_type_t seg_type, gaspi_queue_id_t qid)
{
    gaspi_return_t ret = GASPI_SUCCESS;
    if (north_rank.rank_no != NO_NEIGHBOUR)
    {
        ret = notify_neighbour(north_rank, seg_type, qid);
        check_errors(ret, "notify_all_neighbours", "to north rank");
    }
    if (south_rank.rank_no != NO_NEIGHBOUR)
    {
        ret = notify_neighbour(south_rank, seg_type, qid);
        check_errors(ret, "notify_all_neighbours", "to south rank");
    }
    if (east_rank.rank_no != NO_NEIGHBOUR)
    {
        ret = notify_neighbour(east_rank, seg_type, qid);
        check_errors(ret, "notify_all_neighbours", "to east rank");
    }
    if (west_rank.rank_no != NO_NEIGHBOUR)
    {
        ret = notify_neighbour(west_rank, seg_type, qid);
        check_errors(ret, "notify_all_neighbours", "to west rank");
    }
    return ret;
}

void
finish_map(void)
{
    gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
    gaspi_segment_delete(uid);
    gaspi_segment_delete(vid);
    gaspi_segment_delete(aid);
}

FILE *
get_file_handle(void)
{
    char buffer[128];
    int n = sprintf(buffer, "/home/cherold/info-rank-%d.txt", rank);
    buffer[n] = '\0';
    return fopen(buffer, "w");
}

void
print_map(FILE * out, seg_type_t seg)
{
    double * ptr = (seg == SEG_U) ? u_double_ptr : v_double_ptr;
    for (uint64_t i = 0; i < n; ++i)
    {
        for (uint64_t j = 0; j < m; ++j)
        {
            fprintf(out, "%lf ", ptr[IDX(i, j, m)]);
        }
        fprintf(out, "\n");
    }
}

double
compute_heat(double *a, double *u, double *v)
{
    double ln = 2.0;
    double h = ln / (n - 1);
    double hi = 1.0 / h;
    double hi2 = hi * hi;
    double residual = 0.0;
    double alpha = 22;
    double du;
    double dt = h * h / 4 / alpha;

    int w = m;
    for (int l = 1; l < (n - 1); l++)
    {
        for (int k = 1; k < (m - 1); k++)
        {
            du = (a[IDX(l, k - 1, w)] * u[IDX(l, k - 1, w)] + a[IDX(l, k, w)] * u[IDX(l, k + 1, w)]
                    + a[IDX(l - 1, k, w)] * u[IDX(l - 1, k, w)] + a[IDX(l, k, w)] * u[IDX(l + 1, k, w)]
                    - (a[IDX(l, k - 1, w)] + a[IDX(l, k, w)] + a[IDX(l - 1, k, w)] + a[IDX(l, k, w)]) * u[IDX(l, k, w)])
                    * dt * hi2;
            v[IDX(l, k, w)] = u[IDX(l, k, w)] + du;
            du = MAX(du, -du);
            residual = MAX(residual, du);
        }
    }
    return residual;
}

int
main(int argc, char ** argv)
{
    if (gaspi_proc_init( GASPI_BLOCK) != GASPI_SUCCESS)
    {
        fprintf(stderr, "startup failed\n");
        gaspi_proc_term( GASPI_BLOCK);
        exit(EXIT_FAILURE);
    }

    size_t typesize = sizeof(double);

    gaspi_proc_rank(&rank);
    gaspi_proc_num(&rankcount);

    m = (matrix_width  / rankcount) + 1;
    n = (matrix_height / rankcount) + 1;

    init_map(m, n, typesize);

    double residual = 0.0;
    double eps = 1.0e-4;

    double * t = NULL;
    double * a = a_double_ptr;
    double * u = u_double_ptr;
    double * v = v_double_ptr;

    gaspi_return_t ret = GASPI_SUCCESS;
    seg_type_t seg_type = SEG_V;
    uint64_t step = 0;

    do
    {
        gaspi_wait(0, GASPI_BLOCK);
        residual = compute_heat(a, u, v);
        notify_all_neighbours(seg_type, 0);
        read_all_halos(seg_type, typesize);

        seg_type = (seg_type == SEG_U) ? SEG_V : SEG_U;

        t = u;
        u = v;
        v = t;

        ret = gaspi_allreduce(&residual, &residual, 1, GASPI_OP_MAX, GASPI_TYPE_DOUBLE, GASPI_GROUP_ALL, TIMEOUT);

        check_errors(ret, "gaspi_allreduce", NULL);

        if (step >= 100000)
        {
            gaspi_printf("res %lf\n", residual);
            gaspi_printf("Too many steps\n");
            break;
        }
        ++step;

    } while (residual > eps);

    gaspi_printf("Finished computation with %d steps\n", step);
    gaspi_printf("res %lf\n", residual);

    finish_map();

    gaspi_proc_term( GASPI_BLOCK);

    return EXIT_SUCCESS;
}
