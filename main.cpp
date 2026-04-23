// Created by Alessandro on 16/03/2026 based on code by Joris
//

// This code performs an iterative Monte Carlo procedure to obtain a Maximum Entropy
// model for the chromosome organization of C. crescentus

// Created by Joris on 09/07/2018.
// Modified by Capucine on 31/03/2025
// Modified by Alessandro on 16/03/2026 (add T)
// Modified by Alessandro on 23/04/2026 --> changed indexing of files

// This code performs an iterative Monte Carlo procedure to obtain a Maximum Entropy
// model for the chromosome organization of C. crescentus

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <thread>
#include <mutex>
#include "Initialize.h"
#include "Moves.h"
using namespace Eigen;

//const int number_of_threads = 36;
//const int number_of_threads = 48;


const double diameter = 6.4;
const int length_cylinder = 21;
const int midway = int(ceil(length_cylinder/2));
const int cap_length = 2;

const int pol_length = 1620;

// const long int mc_moves_start = 25000000;
const long int mc_moves_start = 25000000;
long int mc_moves;
//const int burn_in_time = 2000000;
const int burn_in_time = 2000000;
const int save_interval = 500000;



// create a vector of temperature

std::vector<double> linspace(float start, float end, size_t points)
{
  std::vector<float> res(points);
  float step = (end - start) / (points - 1);
  size_t i = 0;
  for (auto& e : res)
  {
    e = start + step * i++;
  }
  return res;
}

std::vector<double> betas = linspace(0.5, 3., 0.1);
// now I am running 1 thread per beta SO number of threads = number of betas
const int number_of_threads = static_cast<int>(betas.size());
std::vector<std::vector<Vector3i>> polymer(number_of_threads);

bool boundary_cond = 1; //enforces boundary conditions if 1
bool orient = 1; //orients the cell such that the origin is always in the left half

std::vector< std::vector<double>> Interaction_E(pol_length, std::vector<double>(pol_length,0));

std::uniform_real_distribution<double> unif(0.0,1.0);
std::uniform_int_distribution<int> unimove(0,2);
std::uniform_int_distribution<int> unisite(0,pol_length-1);
std::vector<std::vector<std::vector<double>>> total_contacts(number_of_threads, std::vector< std::vector<double>>(pol_length, std::vector<double>(pol_length, 0)));
std::vector<std::vector<double>> final_contacts(pol_length, std::vector<double>(pol_length, 0));

void move(std::vector<Vector3i> &polymer,int thread_num, int m, double beta){ //performs a single Monte Carlo step
    int action;
    int site;

    action = unimove(gen);
    site = unisite(gen);
    if (action==0){
        kink_move(polymer,site, thread_num,m, beta);
    }
    else if (action==1){
        crankshaft_move(polymer,site, thread_num, m, beta);
    }
    else if (action==2){
       loop_move(polymer,site, thread_num, m, beta);
    }
}

void run_burnin(int thread_num, int mc_moves, double beta) { //burns in the polymer configurations
    for (int m = 1; m < mc_moves; m++) {
        move(polymer[thread_num], thread_num, m, beta);
    }
}

const std::string base_path = "/home/alessandro/PhD_Alessandro/first_project/MaxEnt-Chromosome-Caulobacter-0.1/Forward_Monte_Carlo_with_T/";


void run(int thread_num, int mc_moves, double beta, int batch) {
    // save stuff every save_interval steps
    for (int m = 1; m < mc_moves; m++) {    //performs a forward polymer simulation
        move(polymer[thread_num], thread_num, m, beta);

        if (m % save_interval == 0) {
            std::string filename = base_path + "intermediate_confs/"
                                 + "conf_batch" + std::to_string(batch)
                                 + "_thread" + std::to_string(thread_num)
                                 + "_step" + std::to_string(m) + ".txt";
            std::ofstream out(filename);
            for (int i = 0; i < pol_length; i++) {
                for (int j = 0; j < 3; j++) {
                    out << polymer[thread_num][i][j] << '\n';
                }
            }
        }
    }
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Started!" << std::endl;

    mc_moves=mc_moves_start;
    // Load interaction energies once
    std::ifstream couplings("/home/alessandro/alessandro/PhD_Alessandro/first_project/MaxEnt-Chromosome-Caulobacter-0.1/Forward_Monte_Carlo/energies_ccrescentus_wt_rep1.txt");
    for (int i = 0; i < pol_length; i++) {
        for (int j = 0; j < pol_length; j++) {
            couplings >> Interaction_E[i][j];
        }
    }
    couplings.close();
    std::cout << "Loaded interaction energies." << std::endl;
    for (int i = 0; i < pol_length / 4; i++) {
        for (int j = i + 2; j < pol_length / 4; j++) {
            std::cout << "Interaction_E[4 * " << i << "][4 * " << j << "]: "
                    << Interaction_E[4 * i][4 * j]
                    << std::endl;
        }
    }
    const int batch_size = number_of_threads;
    //const int total_batches = 1000;
    const int total_batches = 100;
    int sample_counter = 0;
    std::mutex counter_mutex;


    for (int batch = 0; batch < total_batches; batch++) {
        std::cout << "Starting batch " << batch + 1 << " / " << total_batches << std::endl;

        // Initialize polymers
        for (int l = 0; l < batch_size; l++) {
            initialize(polymer[l], pol_length, l);
        }
        std::cout << "Initialized monomer positions (thread 0):\n";
        for (int i = 0; i < 10; ++i) {
            std::cout << polymer[0][i].transpose() << std::endl;
        }

        
        auto finish1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed1 = finish1 - start;
        std::cout << "Elapsed time: " << elapsed1.count() << " seconds\n";
        std::cout << "Finished Initializing "  << std::endl;

        // Burn-in
        std::vector<std::thread> threads(batch_size);
        for (int l = 0; l < batch_size; l++) {
            threads[l] = std::thread(run_burnin, l, burn_in_time, betas[l]);
        }
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start;
        std::cout << "Elapsed time: " << elapsed2.count() << " seconds\n";
        std::cout << "Finished Burn-in "<< std::endl;
        

        for (auto &&l : threads) l.join();

        // Forward simulation
        for (int l = 0; l < batch_size; l++) {
            threads[l] = std::thread(run, l, mc_moves, betas[l], batch);
        }

        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start;
        std::cout << "Elapsed time: " << elapsed3.count() << " seconds\n";
        std::cout << "Finished forward "<< std::endl;
        for (auto &&l : threads) l.join();

        // Output configurations
        auto write_start = std::chrono::high_resolution_clock::now();

        std::vector<std::thread> write_threads;
        // std::mutex counter_mutex;

        for (int thread_num = 0; thread_num < batch_size; thread_num++) {
                 write_threads.emplace_back([&, thread_num]() { // create a new thread and store it in the vector write_threads
                        if (polymer[thread_num][0][2] < 1000) {
                        std::string filename = "configuration_batch" + std::to_string(batch) + "_thread" + std::to_string(thread_num) + ".txt";    
                        std::ofstream out(base_path + "final_confs/final_configuration_" +
                                        filename);

                                for (int i = 0; i < pol_length; i++) {
                                        for (int j = 0; j < 3; j++) {
                                                out << polymer[thread_num][i][j] << '\n';
                                         }
                                 }

                                 // Safely increment counter
                                 // std::lock_guard<std::mutex> lock(counter_mutex);
                                 // sample_counter++;
                         }
                        else {
                            std::string filename1 = "configuration_batch" + std::to_string(batch) + "_thread" + std::to_string(thread_num) + ".txt";    
                            std::ofstream out(base_path + "rejected_confs/final_configuration_" +
                                            filename1);

                                for (int i = 0; i < pol_length; i++) {
                                        for (int j = 0; j < 3; j++) {
                                                out << polymer[thread_num][i][j] << '\n';
                                         }
                                 }

                                 // Safely increment counter
                                 // std::lock_guard<std::mutex> lock(counter_mutex);
                                 // sample_counter++;
                        }

                });
         }

        // Join all output threads
         for (auto& t : write_threads) t.join();

        // ⏱  Stop timing
         auto write_end = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double> write_elapsed = write_end - write_start;
         std::cout << "Finished writing files in " << write_elapsed.count() << " seconds\n";
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Completed all batches." << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;

}







// int main() {
//     auto start = std::chrono::high_resolution_clock::now();
//     std::cout << "started! " << std::endl;

//     for (int l = 0; l < number_of_threads; l++) {
//         initialize(polymer[l], pol_length, l);
//     }
//     std::cout << "Initialized " << std::endl;

//     // Read in converged interaction energies
//     std::ifstream couplings;
//     couplings.open("/home/capucine/Documents/test/Data/Initial_Energies/energies_ccrescentus_wt_rep1.txt");
//     //couplings.open ("/media/joris/raid_data/Chromosome_Maxent/Saved_energies/energies_smallcell.txt");
//     for (int i = 0; i < pol_length; i++) { //read in energies
//         for (int j = 0; j < pol_length; j++) {
//             couplings >> Interaction_E[i][j];
//         }
//     }
//     couplings.close();

//     // burn in configurations
//     std::vector<std::thread> threads(number_of_threads);
//     for (auto l = 0; l < number_of_threads; l++) {
//         threads[l] = std::thread(run_burnin, l, burn_in_time);
//     }
//     for (auto &&l : threads) {
//         l.join();
//     }

//     std::cout << "Done with burn in " << std::endl;

//     for (int i = 0; i < pol_length / 4; i++) {
//         for (int j = i + 2; j < pol_length / 4; j++) {
//             std::cout << "Interaction_E[4 * " << i << "][4 * " << j << "]: "
//                     << Interaction_E[4 * i][4 * j]
//                     << std::endl;
//         }
//     }


//     //run forward simulation
//     for (auto l = 0; l < number_of_threads; l++) {
//         threads[l] = std::thread(run, l, mc_moves);
//     }
//     for (auto &&l : threads) {
//         l.join();
//     }


//     std::cout << "Finished! "  << std::endl;

//     //output final configurations for each thread
//     int n_capture = 0;
//     while (n_capture <number_of_threads) {
//         for (int thread_num = 0; thread_num < number_of_threads; thread_num++) {
//             if (polymer[thread_num][0][2] < 20) {
//                 std::ofstream final_configuration;
//                 final_configuration.open(
//                         "/home/capucine/Documents/test/Data/Final_Configurations/configuration_init_" +
//                         std::to_string(n_capture) + ".txt");
//                 for (int i = 0; i < pol_length; i++) { //output contact frequencies
//                     for (int j = 0; j < 3; j++) {
//                         final_configuration << polymer[thread_num][i][j] << std::endl;
//                     }
//                 }
//                 final_configuration.close();
//                 n_capture++;
//             }
//         }
//     }




//     auto finish = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = finish - start;
//     std::cout << "Elapsed time: " <<  elapsed.count() << " seconds\n";
//     return 0;
// } 


