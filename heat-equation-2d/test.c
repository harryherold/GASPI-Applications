#include "test.h"

void
test_read_halo(uint64_t m, uint64_t n, size_t typesize)
{
    gaspi_queue_id_t queue_id = 0;
	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

	if(west_rank.rank_no != NO_NEIGHBOUR)
	{
		read_halo(west_rank, west_halo, SEG_V, typesize, queue_id);
		gaspi_wait(queue_id, GASPI_BLOCK);

		gaspi_pointer_t seg_ptr;
		gaspi_segment_ptr(vid, &seg_ptr);
		for(int i = 0; i < west_rank.element_count ; ++i)
		{
			double * v = (double *)(seg_ptr + west_halo.halo_offset + (west_halo.jump_offset * i));
			if(*v != (west_rank.rank_no + 1))
			{
				gaspi_printf("Wrong value in element %lu\n",i);
			}
		}
	}

	if(east_rank.rank_no != NO_NEIGHBOUR)
	{
		read_halo(east_rank, east_halo, SEG_V, typesize, queue_id);
		gaspi_wait(queue_id, GASPI_BLOCK);

		gaspi_pointer_t seg_ptr;
		gaspi_segment_ptr(vid, &seg_ptr);

		for (int i = 0; i < east_rank.element_count; ++i)
		{
			double * v = (double *) (seg_ptr + east_halo.halo_offset + (east_halo.jump_offset * i));
			if (*v != (east_rank.rank_no + 1))
			{
				gaspi_printf("Wrong value in element %lu\n", i);
			}
		}
	}

	if(north_rank.rank_no != NO_NEIGHBOUR)
	{
		read_halo(north_rank, north_halo, SEG_V, typesize, queue_id);
		gaspi_wait(queue_id, GASPI_BLOCK);

		gaspi_pointer_t seg_ptr;
		gaspi_segment_ptr(vid, &seg_ptr);

		for (int i = 0; i < north_rank.element_count; ++i)
		{
			double * v = (double *) (seg_ptr + north_halo.halo_offset + (north_halo.jump_offset * i));
			if (*v != (north_rank.rank_no + 1))
			{
				gaspi_printf("Wrong value in element %lu\n", i);
			}
		}
	}

	if(south_rank.rank_no != NO_NEIGHBOUR)
	{
		read_halo(south_rank, south_halo, SEG_V, typesize, queue_id);
		gaspi_wait(queue_id, GASPI_BLOCK);

		gaspi_pointer_t seg_ptr;
		gaspi_segment_ptr(vid, &seg_ptr);

		gaspi_printf("in south rank %d\n",south_rank.rank_no);
		for (int i = 0; i < south_rank.element_count; ++i)
		{
			double * v = (double *) (seg_ptr + south_halo.halo_offset + (south_halo.jump_offset * i));
			if (*v != (south_rank.rank_no + 1))
			{
				gaspi_printf("Wrong value in element %lu\n", i);
			}
		}
	}

	//~ FILE * out = get_file_handle();
	//~ if (out == NULL)
	//~ {
		//~ gaspi_printf("Get no file handle\n");
		//~ return;
	//~ }
	//~ print_map(out, SEG_V);
	//~ fclose(out);

	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
}

void test_notify_all()
{
	seg_type_t seg_type = SEG_U;
	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
	notify_all_neighbours(seg_type, 0);

	gaspi_return_t ret = GASPI_SUCCESS;
	gaspi_segment_id_t seg = SEG_ID(seg_type);

	FILE * out = get_file_handle("test-halo");

	for (uint64_t i = 0; i < neighbour_count; ++i)
	{
		gaspi_notification_id_t notify_id;
		gaspi_notification_t    notify_val;

		ret = blocking_waitsome(NOTIFY_ID_BEGIN, NOTIFY_ID_END, &notify_id, &notify_val, seg);
		check_errors(ret, "test_notify_all", "blocking_waitsome");

		if (check_notify(north_rank, notify_id, notify_val))
		{
			fprintf(out, "get notify from north rank %lu\n",north_rank.rank_no);
		}
		if (check_notify(south_rank, notify_id, notify_val))
		{
			fprintf(out, "get notify from south rank %lu\n",south_rank.rank_no);
			gaspi_printf("south rank\n");
		}
		if (check_notify(east_rank, notify_id, notify_val))
		{
			gaspi_printf("east rank\n");
			fprintf(out, "get notify from east rank %lu\n",east_rank.rank_no);
		}
		if (check_notify(west_rank, notify_id, notify_val))
		{
			fprintf(out, "get notify from west rank %lu\n",west_rank.rank_no);
		}
	}

	fclose(out);
	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
}

