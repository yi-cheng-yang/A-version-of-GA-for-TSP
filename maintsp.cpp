#include "GATSP.hpp"

int main(int argc, char* argv[]){
    int run = std::atoi(argv[1]);
    int NE = std::atoi(argv[2]);
    int population = std::atoi(argv[3]);
    GTSP obj;
    obj.exe(run, NE, population);
}