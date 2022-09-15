SET PATH=c:\MinGW\bin;%PATH%
REM  no warnings from: -Waliasing -Wline-truncation -Wsurprising -Wunderflow
REM  -Wzerotrip not recognised
REM gfortran -fmax-errors=20 -fcheck=all -ffree-line-length-none   -Wunderflow -static -frecursive -c f_e_interpreter.f90
ECHO 'compiling f_electron'
gfortran -fmax-errors=20 -fcheck=all -static -c diagr.for
gfortran -fmax-errors=20  -fno-range-check -ffixed-line-length-none  -Wunderflow -static -c f_e_ESO.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_check.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_readFile.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_Block.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_calculate.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_data.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_Edipole.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_fit.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_fncrosserrors.f90
gfortran -fmax-errors=20 -fcheck=all -ffpe-trap=invalid,zero,overflow  -ffixed-line-length-none  -Wunderflow -static -c f_e_group.f90
gfortran -fmax-errors=20 -fcheck=all -ffpe-trap=invalid,zero,overflow  -ffixed-line-length-none  -Wunderflow -static -c f_e_input.f90
gfortran -fmax-errors=20 -fcheck=all -ffree-line-length-none   -Wunderflow -static -frecursive -c f_e_interpreter.f90
gfortran -fmax-errors=20 -fcheck=all -ffpe-trap=invalid,zero,overflow  -ffixed-line-length-none  -Wunderflow -static -c f_e_LF.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_magnetics.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_parameters.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static -c f_e_wigner.f90
gfortran -fmax-errors=20 -fcheck=all -ffixed-line-length-none  -Wunderflow -static f_electrons.f90 f_e_check.o f_e_block.o f_e_calculate.o f_e_data.o f_e_Edipole.o f_e_ESO.o f_e_fit.o f_e_fncrosserrors.o f_e_group.o f_e_input.o f_e_interpreter.o f_e_LF.o f_e_magnetics.o f_e_parameters.o f_e_readFile.o f_e_wigner.o diagr.o -o f_e

