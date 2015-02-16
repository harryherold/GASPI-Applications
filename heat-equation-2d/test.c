#include "test.h"

void
test_read_halo(uint64_t m, uint64_t n, size_t typesize)
{
	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

	if(west_rank.rank_no != NO_NEIGHBOUR)
	{
		read_halo(west_rank, west_halo, SEG_V, typesize);
		gaspi_wait(0, GASPI_BLOCK);

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
		read_halo(east_rank, east_halo, SEG_V, typesize);
		gaspi_wait(0, GASPI_BLOCK);

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
		read_halo(north_rank, north_halo, SEG_V, typesize);
		gaspi_wait(0, GASPI_BLOCK);

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
		read_halo(south_rank, south_halo, SEG_V, typesize);
		gaspi_wait(0, GASPI_BLOCK);

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

	FILE * out = get_file_handle();
	if (out == NULL)
	{
		gaspi_printf("Get no file handle\n");
		return;
	}
	print_map(out, SEG_V);
	fclose(out);

	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
}
