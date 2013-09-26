AM_CPPFLAGS = -I$(srcdir)/../src
bin_PROGRAMS = lines test parallel continue

lines_SOURCES = lines.c
lines_CPPFLAGS = -I$(HOME)/local/include -I$(srcdir)/../src
lines_LDADD = -L$(HOME)/local/lib ../src/libmcmc.la -lhdf5 -lhdf5_hl -lgsl -lgslcblas

continue_SOURCES = continue.c
continue_CPPFLAGS = -I$(HOME)/local/include -I$(srcdir)/../src
continue_LDADD = -L$(HOME)/local/lib ../src/libmcmc.la -lhdf5 -lhdf5_hl -lgsl -lgslcblas

test_sources = test.c
test_CPPFLAGS = $(MPI_CFLAGS) -I$(HOME)/local/include/ -I$(srcdir)/../src
test_LDADD = $(MPI_CLDFLAGS)

parallel_sources = parallel.c

parallel_CPPFLAGS = -I$(srcdir)/../src $(MPI_CFLAGS) -I$(HOME)/local/include
parallel_LDADD = -L$(HOME)/local/lib ../src/libmcmc.la -lhdf5 -lhdf5_hl $(MPI_CLDFLAGS) -lgsl -lgslcblas
