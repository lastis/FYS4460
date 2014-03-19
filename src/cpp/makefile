#Basic makefile test

main: main3.cpp CPhys.a; g++ main3.cpp CPhys.a -fopenmp -O3
    
init: init_system2.cpp CPhys.a; g++ -o init.x init_system2.cpp CPhys.a

thermalize: thermalize.cpp CPhys.a; g++ -o thermalize.x thermalize.cpp CPhys.a -fopenmp -O3
    
clean: ;rm a.out init.x thermalize.x

fclean: ;rm a.out init.x ../../res/Measurements/* ../../res/States/*

