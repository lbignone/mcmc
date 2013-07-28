AM_CPPFLAGS = -I$(srcdir)/../src
bin_PROGRAMS = lines test parallel

lines_SOURCES = lines.c
lines_LDADD = ../src/libmcmc.la -lhdf5 -lhdf5_hl

test_sources = test.c
test_CPPFLAGS = $(MPI_CFLAGS)
test_LDADD = $(MPI_CLDFLAGS)

parallel_sources = parallel.c

parallel_CPPFLAGS = -I$(srcdir)/../src $(MPI_CFLAGS)
parallel_LDADD = ../src/libmcmc.la -lhdf5 -lhdf5_hl $(MPI_CLDFLAGS)
