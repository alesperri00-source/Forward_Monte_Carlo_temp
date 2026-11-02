//
// Created by Joris on 09/07/2018.
//
//
#ifndef Inverse_3D_new_Initialize_h
#define Inverse_3D_new_Initialize_h

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <string>
#include "global.h"
using namespace Eigen;

std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts(number_of_threads);
std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations(number_of_threads);

void initialize(std::vector<Vector3i> &polymer, int pol_length, int thread_num){
    if (init_from_file) { // configuration sampled from the MaxEnt model
        std::string filename = init_conf_path + "configuration_init_" + std::to_string(thread_num) + ".txt";
        std::ifstream configuration(filename);
        int coordinate;
        for(int i = 0; i < pol_length; i++){
            Vector3i monomer;
            for(int j = 0; j < 3; j++){
                configuration >> coordinate;
                monomer(j) = coordinate;
            }
            polymer.push_back(monomer);
        }
        if (!configuration) {
            std::cerr << "ERROR: could not read " << pol_length << " monomers from " << filename << std::endl;
            exit(1);
        }
        for (int i = 0; i < pol_length; i++) {
            if ((polymer[(i+1)%pol_length] - polymer[i]).squaredNorm() != 1) {
                std::cerr << "ERROR: monomers " << i << " and " << (i+1)%pol_length << " are not neighbours in " << filename << std::endl;
                exit(1);
            }
        }
    }
    else { // fixed compact conformation
        int x_i, y_i, z_i;
        Vector3i monomer;
        x_i=2; y_i=2; z_i=2;
        for (int i=0; i<pol_length; i++){
                if (i%4==0){x_i++;}
                else if (i%4==1){y_i++;}
                else if (i%4==2){x_i--;}
                else if (i%4==3){y_i--;}
            monomer << x_i,y_i,z_i; polymer.push_back(monomer);
        }
    }

    for (int i = 0; i < pol_length; i++){ // register occupied lattice sites and initial contacts
        std::vector<int> site = {polymer[i](0), polymer[i](1), polymer[i](2)};
        if (locations[thread_num].find(site)!=locations[thread_num].end()){
            for (auto elem : locations[thread_num].find(site)->second){
                contacts[thread_num][{elem,i}]=0;
            }
            locations[thread_num].find(site)->second.push_back(i);
        }
        else { locations[thread_num][site] = {i}; }
    }
}
#endif //Inverse_3D_new_Initialize.h
