Name = PNP
FFLAGS = -O3
OBJS =  	para.f90\
			data.f90\
			main_diff.f90\
			Transport_parameter.f90\
			main_F.f90\
			main.f90\
			initial.f90\
			diff_error.f90\
			F_results.f90\
			diff_results.f90\

$(Name): $(OBJS)
	gfortran $(FFLAGS) -fopenmp -ffree-line-length-none -o $@ $(OBJS)

clean:
	rm *.o
###### gfortran $(FFLAGS) -fopenmp -ffree-line-length-none -o $@ $(OBJS)
