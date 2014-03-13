#Basic makefile test

main: main3.cpp CPhys.a; g++ main3.cpp CPhys.a -fopenmp
    
init: init_system2.cpp CPhys.a; g++ -o init.x init_system2.cpp CPhys.a
    
clean: ;rm a.out init.x