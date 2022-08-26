#
#

FC = gfortran
# FC = f90
FFLAGS =

CWD = ./
MYL = ../lib/
MATH = dsp_math_MATRIX

OUT1 = re.exe
OUT2 = zin.exe

INC =
LIB =
# LIB =	-lmatmpp_sc -lMSL2

re:	$(CWD)re_header.f90\
	$(CWD)re_verification.f90\
	$(CWD)re_fbi_inhomo.f90\
	$(CWD)re_mag.f90\
	$(CWD)re_slv.f90\
	$(CWD)re_mi_couple.f90\
	$(CWD)re_main.f90

	$(FC) $(FFLAGS) -c $(CWD)re_header.f90
	$(FC) $(FFLAGS) -c $(CWD)re_verification.f90
	$(FC) $(FFLAGS) -c $(CWD)re_fbi_inhomo.f90
	$(FC) $(FFLAGS) -c $(CWD)re_mag.f90
	$(FC) $(FFLAGS) -c $(CWD)re_slv.f90
	$(FC) $(FFLAGS) -c $(CWD)re_mi_couple.f90
	$(FC) $(FFLAGS) -c $(CWD)re_main.f90

	$(FC) $(FFLAGS) $(CWD)re_header.o\
			$(CWD)re_verification.o\
			$(CWD)re_fbi_inhomo.o\
			$(CWD)re_mag.o\
			$(CWD)re_slv.o\
			$(CWD)re_mi_couple.o\
			$(CWD)re_main.o\
			-o $(OUT1) $(LIB)

zin: $(CWD)re_header.f90\
	$(CWD)re_verification.f90\
	$(CWD)re_fbi_inhomo.f90\
	$(CWD)re_mag.f90\
	$(CWD)re_zinsearch.f90\

	$(FC) $(FFLAGS) -c $(CWD)re_header.f90
	$(FC) $(FFLAGS) -c $(CWD)re_verification.f90
	$(FC) $(FFLAGS) -c $(CWD)re_fbi_inhomo.f90
	$(FC) $(FFLAGS) -c $(CWD)re_mag.f90
	$(FC) $(FFLAGS) -c $(CWD)re_zinsearch.f90

	$(FC) $(FFLAGS) $(CWD)re_header.o\
			$(CWD)re_verification.o\
			$(CWD)re_fbi_inhomo.o\
			$(CWD)re_mag.o\
			$(CWD)re_zinsearch.o\
			-o $(OUT2) $(LIB)


clean:
	rm -f ./*.o ./*.mod ./*.log ./*.exe ./*.err ./*.out

clear:
	rm -f ./*.o ./*.mod ./*.log ./*.err ./*.out
