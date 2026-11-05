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
#include <cstdlib>
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
const long int mc_moves_start = 60000000;
long int mc_moves;
//const int burn_in_time = 2000000;
//long int burn_in_time = 200000000; // 0.9-2.0 run: too short for the z profile, which kept relaxing for ~6 batches
long int burn_in_time = 2000000000; // 2G moves per temperature, ~17 min at 30 s per 60M moves

// create a vector of temperature

std::vector<double> linspace(float start, float end, size_t points)
{
  std::vector<double> res(points);
  float step = (end - start) / (points - 1);
  size_t i = 0;
  for (auto& e : res)
  {
    e = start + step * i++;
  }
  return res;
}


std::vector<double> inverse_vector(const std::vector<double>& Ts)
{
    std::vector<double> betas(Ts.size());

    for (size_t i = 0; i < Ts.size(); i++)
    {
        betas[i] = 1.0 / Ts[i];
    }

    return betas;
}

// Temperature ladder, from lowest to highest T (THESE ARE TEMPERATURES).
// Below 0.9, Var(E) is up to 30 times larger than above 1.0, so evenly spaced temperatures would barely swap there.
// The 28 temperatures from 0.5 to 0.9 are spaced so that every neighbouring pair has the same predicted swap acceptance,
// about 0.6, from A = erfc(dbeta * sigma_E / 2) with sigma_E measured in the 0.5-1.5 benchmark.
// 0.9, 0.957895 and 1.01579 are temperatures of the 0.9-2.0 run (linspace(0.9, 2.0, 20)), where the two runs can be
// compared and joined.
// num_points must stay a literal: Initialize.h sizes contacts and locations with number_of_threads before Ts exists,
// which works only because number_of_threads is a compile-time constant.
const int num_points = 30;
//std::vector<double> Ts = linspace(0.9, 2.0, num_points); // 0.9-2.0 run, with num_points = 20
std::vector<double> Ts = {
    0.5000, 0.5183, 0.5353, 0.5512, 0.5656, 0.5798, 0.5942, 0.6079, 0.6219, 0.6361,
    0.6492, 0.6623, 0.6749, 0.6863, 0.6974, 0.7086, 0.7200, 0.7338, 0.7481, 0.7625,
    0.7773, 0.7924, 0.8063, 0.8207, 0.8356, 0.8519, 0.8703, 0.9000, 0.957895, 1.01579};
std::vector<double> betas = inverse_vector(Ts);

// now I am running 1 thread per beta SO number of threads = number of betas
const int number_of_threads = num_points;
std::vector<std::vector<Vector3i>> polymer(number_of_threads);
std::vector<long long> accepted_moves(number_of_threads, 0);
std::vector<long long> attempted_moves(number_of_threads, 0);
bool boundary_cond = 1; //enforces boundary conditions if 1
bool orient = 1; //orients the cell such that the origin is always in the left half
bool init_from_file = 1; //starts from the MaxEnt configurations init_conf_path/configuration_init_<thread>.txt if 1, from a fixed compact conformation if 0

// Replica exchange (parallel tempering): every exchange_interval moves, neighbouring temperatures try to swap their
// configurations, so that a configuration stuck at low T can travel up to high T, decorrelate there, and come back down.
bool replica_exchange = 1; //swaps configurations between neighbouring temperatures if 1
const long int exchange_interval = 100000; // MC moves per temperature between two rounds of swap attempts
std::mt19937_64 swap_generator;
std::vector<int> replica_at(number_of_threads); // replica_at[l] = label of the configuration currently at temperature Ts[l]
std::vector<int> trip_state(number_of_threads, 0); // per configuration: 0 never at lowest T, 1 left lowest T, 2 reached highest T since
std::vector<long long> swap_attempts(number_of_threads - 1, 0); // entry l refers to the pair (Ts[l], Ts[l+1])
std::vector<long long> swap_accepts(number_of_threads - 1, 0);
long long round_trips = 0; // lowest T -> highest T -> lowest T, summed over all configurations
long long exchange_round = 0;
// Flow: per temperature, the number of swap rounds in which the configuration there had last visited the lowest T
// (moving up) or the highest T (moving down). f_up = up / (up + down) goes from 1 at the lowest T to 0 at the highest;
// a steep drop between two temperatures marks a bottleneck.
std::vector<long long> flow_up(number_of_threads, 0);
std::vector<long long> flow_down(number_of_threads, 0);
// After every trace_every swap rounds: the energy and the configuration label at each temperature
std::ofstream trace_out;
const int trace_every = 1;

std::vector< std::vector<double>> Interaction_E(pol_length, std::vector<double>(pol_length,0));
std::uniform_real_distribution<double> unif(0.0,1.0);
std::uniform_int_distribution<int> unimove(0,2);
std::uniform_int_distribution<int> unisite(0,pol_length-1);
//std::vector<std::vector<std::vector<double>>> total_contacts(number_of_threads, std::vector< std::vector<double>>(pol_length, std::vector<double>(pol_length, 0)));
//std::vector<std::vector<double>> final_contacts(pol_length, std::vector<double>(pol_length, 0));
std::vector<std::mt19937_64> generators(number_of_threads);


struct EnergyStats {
    double mean;
    double variance;
};

EnergyStats compute_energy_stats(const std::vector<double>& energies)
{
    if (energies.empty()) {
        return {0.0, 0.0};
    }

    double mean = 0.0;
    double M2 = 0.0;
    size_t n = 0;

    for (double E : energies) {
        ++n;

        double delta = E - mean;
        mean += delta / static_cast<double>(n);

        double delta2 = E - mean;
        M2 += delta * delta2;
    }

    EnergyStats result;
    result.mean = mean;

    // <(E - <E>)^2>
    result.variance = M2 / static_cast<double>(n);

    return result;
}


void move(std::vector<Vector3i> &polymer,int thread_num, long int m, double beta){ //performs a single Monte Carlo step
    int action;
    int site;

    action = unimove(generators[thread_num]);
    site = unisite(generators[thread_num]);

    if (action==0){
        kink_move(polymer,site, thread_num, m, beta);
    }
    else if (action==1){
        crankshaft_move(polymer,site, thread_num, m, beta);
    }
    else if (action==2){
       loop_move(polymer,site, thread_num, m, beta);
    }
    attempted_moves[thread_num]++;
}

const std::string base_path = "/home/alessandro/alessandro/PhD_Alessandro/first_project/MaxEnt-Chromosome-Caulobacter-0.1/Forward_Monte_Carlo_with_T/";
const std::string init_path = "/home/alessandro/alessandro/PhD_Alessandro/first_project/MaxEnt-Chromosome-Caulobacter-0.1/Inverse_Monte_Carlo/";

const std::string init_conf_path = init_path + "sim_data/final_confs/";


void run(int thread_num, long int n_moves, double beta) { //performs n_moves Monte Carlo steps at fixed temperature
    for (long int m = 0; m < n_moves; m++) {
        move(polymer[thread_num], thread_num, m, beta);
    }
}

double instantaneous_energy(int thread_num)
// computes the energy by doing a scalar prodcuct bwtween epsilon and contacts
{
    double E = 0.0;

    // Each entry corresponds to one lattice position
    for (const auto& entry : locations[thread_num]) {

        const std::vector<int>& monomers = entry.second;

        // Every pair occupying the same lattice position is a contact
        for (size_t a = 0; a < monomers.size(); ++a) {
            for (size_t b = a + 1; b < monomers.size(); ++b) {

                int i = monomers[a];
                int j = monomers[b];

                // Follow the i < j convention
                if (i > j)
                    std::swap(i, j);

                E += Interaction_E[i][j];
            }
        }
    }

    return E;
}

void attempt_swaps() { //one round of replica exchange between neighbouring temperatures
    std::vector<double> E(number_of_threads);
    for (int l = 0; l < number_of_threads; ++l) {
        E[l] = instantaneous_energy(l);
    }

    // deterministic even/odd scheme: even rounds try the pairs (0,1),(2,3),..., odd rounds (1,2),(3,4),...
    for (int l = exchange_round % 2; l + 1 < number_of_threads; l += 2) {
        swap_attempts[l]++;
        double log_acc = (betas[l] - betas[l + 1]) * (E[l] - E[l + 1]);
        if (log_acc >= 0 || unif(swap_generator) < exp(log_acc)) {
            // O(1) swaps of the containers: the configurations change temperature, threads and generators do not
            std::swap(polymer[l], polymer[l + 1]);
            std::swap(locations[l], locations[l + 1]);
            std::swap(contacts[l], contacts[l + 1]);
            std::swap(replica_at[l], replica_at[l + 1]);
            std::swap(E[l], E[l + 1]); // E[l] stays the energy of the configuration now at Ts[l]
            swap_accepts[l]++;
        }
    }
    exchange_round++;

    int coldest = replica_at[0];
    int hottest = replica_at[number_of_threads - 1];
    if (trip_state[coldest] == 2) round_trips++;
    trip_state[coldest] = 1;
    if (trip_state[hottest] == 1) trip_state[hottest] = 2;

    for (int l = 0; l < number_of_threads; ++l) {
        int state = trip_state[replica_at[l]];
        if (state == 1) flow_up[l]++;
        else if (state == 2) flow_down[l]++;
    }

    if (trace_out.is_open() && exchange_round % trace_every == 0) {
        trace_out << exchange_round;
        for (int l = 0; l < number_of_threads; ++l) trace_out << '\t' << E[l];
        for (int l = 0; l < number_of_threads; ++l) trace_out << '\t' << replica_at[l];
        trace_out << '\n';
    }
}

void advance_all(long int n_moves) { //advances every temperature by n_moves MC moves, one thread per temperature
    long int done = 0;
    while (done < n_moves) {
        long int chunk = replica_exchange ? std::min(exchange_interval, n_moves - done) : n_moves;
        std::vector<std::thread> threads;
        for (int l = 0; l < number_of_threads; ++l) {
            threads.emplace_back(run, l, chunk, betas[l]);
        }
        for (auto& t : threads) {t.join();}
        done += chunk;
        if (replica_exchange) attempt_swaps();
    }
}

void print_swap_stats(std::ostream& out) {
    out << "# T_low\tT_high\tswap_attempts\tswap_acceptance\n";
    for (int l = 0; l + 1 < number_of_threads; ++l) {
        double rate = swap_attempts[l] > 0 ? static_cast<double>(swap_accepts[l]) / static_cast<double>(swap_attempts[l]) : 0.0;
        out << Ts[l] << '\t' << Ts[l + 1] << '\t' << swap_attempts[l] << '\t' << rate << '\n';
    }
    out << "# round trips (lowest T -> highest T -> lowest T): " << round_trips << '\n';
}

void print_flow_stats(std::ostream& out) {
    out << "# T\trounds_moving_up\trounds_moving_down\tf_up\n";
    for (int l = 0; l < number_of_threads; ++l) {
        long long n = flow_up[l] + flow_down[l];
        out << Ts[l] << '\t' << flow_up[l] << '\t' << flow_down[l] << '\t';
        if (n > 0) out << static_cast<double>(flow_up[l]) / static_cast<double>(n) << '\n';
        else out << "nan\n"; // no configuration there has visited the lowest T yet
    }
}

void open_trace(const std::string& filename) {
    trace_out.close();
    trace_out.open(base_path + filename);
    if (!trace_out.is_open()) {
        std::cerr << "ERROR: could not open " << base_path + filename << std::endl;
        exit(1);
    }
    trace_out << "# swap round, then the energy at each temperature (lowest to highest T), then the label of the "
                 "configuration at each temperature, after the swaps of that round; one round = "
              << exchange_interval << " moves per temperature\n";
}

int main() {
    if (static_cast<int>(Ts.size()) != number_of_threads) {
        std::cerr << "ERROR: the ladder has " << Ts.size() << " temperatures but num_points = " << num_points << std::endl;
        return 1;
    }
    for (int l = 0; l + 1 < number_of_threads; ++l) {
        if (!(Ts[l] < Ts[l + 1])) {
            std::cerr << "ERROR: temperatures must increase, but Ts[" << l << "] = " << Ts[l] << " and Ts[" << l + 1
                      << "] = " << Ts[l + 1] << std::endl;
            return 1;
        }
    }

    //system(("mkdir -p " + base_path + "intermediate_confs").c_str());
    system(("mkdir -p " + base_path + "final_confs").c_str());
    system(("mkdir -p " + base_path + "burnin_confs").c_str());
    system(("mkdir -p " + base_path + "rejected_confs").c_str());

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Started!" << std::endl;

    time_t time_now = time(0);

    for (int i = 0; i < number_of_threads; i++) {
        generators[i].seed(static_cast<unsigned long long>(time_now) + i);
        replica_at[i] = i;
    }
    swap_generator.seed(static_cast<unsigned long long>(time_now) + number_of_threads);
    mc_moves = mc_moves_start;

    // Load interaction energies once
    std::ifstream couplings("/home/alessandro/alessandro/PhD_Alessandro/first_project/MaxEnt-Chromosome-Caulobacter-0.1/Forward_Monte_Carlo/energies_ccrescentus_wt_rep1.txt");
    for (int i = 0; i < pol_length; i++) {
        for (int j = 0; j < pol_length; j++) {
            couplings >> Interaction_E[i][j];
            // Interaction_E[i][j] = 0; // to simulate random polymer
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
    const int total_batches = 800;

    std::vector<std::vector<double>> energy_samples(number_of_threads);
    for (int l = 0; l < number_of_threads; ++l) {
        energy_samples[l].reserve(total_batches);
    }

    //int sample_counter = 0;
    //std::mutex counter_mutex;
    std::cout << "Running with " << number_of_threads << " threads." << std::endl;

    // Initialize polymers
    for (int l = 0; l < batch_size; l++) 
    {
        initialize(polymer[l], pol_length, l);
    }

    std::cout << "Initialized monomer positions (thread 0):\n";
    for (int i = 0; i < 10; ++i) 
    {
        std::cout << polymer[0][i].transpose() << std::endl;
    }
    auto finish1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = finish1 - start;
    std::cout << "Elapsed time: " << elapsed1.count() << " seconds\n";
    std::cout << "Finished Initializing "  << std::endl;


    // Burn in
    std::vector<std::thread> write_threads_bi;

    // 1M moves = 10 swap rounds between saved checkpoints; which configuration is where at every round is in
    // burnin_energy_trace.txt
    const long int burnin_save_interval = 1000000;
    const int burnin_checkpoints = burn_in_time / burnin_save_interval;

    // Settings of this run, so that it can be reproduced
    {
        std::ofstream info(base_path + "run_info.txt");
        info << "seed\t" << time_now << "\t(generator of temperature l: seed + l; swap generator: seed + "
             << number_of_threads << ")\n";
        info << "temperatures";
        for (int l = 0; l < number_of_threads; ++l) info << '\t' << Ts[l];
        info << "\nreplica_exchange\t" << replica_exchange << "\nexchange_interval\t" << exchange_interval
             << "\nburn_in_time\t" << burn_in_time << "\nburnin_save_interval\t" << burnin_save_interval
             << "\nmoves_per_batch\t" << mc_moves << "\ntotal_batches\t" << total_batches
             << "\nboundary_cond\t" << boundary_cond << "\norient\t" << orient << "\ninit_from_file\t" << init_from_file
             << "\ninit_conf_path\t" << init_conf_path << '\n';
    }

    open_trace("burnin_energy_trace.txt");

    for (int checkpoint = 0; checkpoint < burnin_checkpoints; checkpoint++) {

        std::cout << "Burn-in checkpoint "<< checkpoint + 1 << " / " << burnin_checkpoints << std::endl;

        // Advance every temperature by another burnin_save_interval moves
        advance_all(burnin_save_interval);

        // Now polymer is not being modified, so it is safe to write
        std::vector<std::thread> write_threads_bi;

        for (int th_n = 0; th_n < batch_size; th_n++) {

            write_threads_bi.emplace_back([&, th_n, checkpoint]() {

                std::string filename =
                    "burnin_configuration_" + std::to_string(checkpoint) +"_thread" + std::to_string(th_n) +".txt";

                std::ofstream out(base_path + "burnin_confs/" + filename);

                if (!out.is_open()) {
                    std::cerr << "ERROR: could not open " << base_path + "burnin_confs/" + filename << '\n';
                    return;
                }

                for (int i = 0; i < pol_length; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        out << polymer[th_n][i][j] << '\n';
                    }
                }
            });
        }

        for (auto& t : write_threads_bi) {t.join();}
    }
    
    auto finish2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = finish2 - start;
    std::cout << "Elapsed time: " << elapsed2.count() << " seconds\n";
    std::cout << "Finished Burn-in "<< std::endl;     

    if (replica_exchange) {
        print_swap_stats(std::cout);
        print_flow_stats(std::cout);
        std::ofstream burnin_swap_out(base_path + "burnin_swap_statistics.txt");
        print_swap_stats(burnin_swap_out);
        std::ofstream burnin_flow_out(base_path + "burnin_flow_statistics.txt");
        print_flow_stats(burnin_flow_out);
    }

    std::fill(accepted_moves.begin(), accepted_moves.end(), 0);
    std::fill(attempted_moves.begin(), attempted_moves.end(), 0);
    std::fill(swap_attempts.begin(), swap_attempts.end(), 0);
    std::fill(swap_accepts.begin(), swap_accepts.end(), 0);
    std::fill(flow_up.begin(), flow_up.end(), 0);
    std::fill(flow_down.begin(), flow_down.end(), 0);
    round_trips = 0;

    open_trace("energy_trace.txt");

    // which configuration label each saved file (thread<l> = temperature Ts[l]) comes from
    std::ofstream replica_out(base_path + "replica_ids.txt");
    replica_out << "# batch\tlabel of the configuration at each temperature, from lowest to highest T\n";

    // Forward simulation

    for (int batch = 0; batch < total_batches; batch++) {
        std::cout << "Starting batch " << batch + 1 << " / " << total_batches << std::endl;
        advance_all(mc_moves);

        // One instantaneous energy sample per temperature
        for (int l = 0; l < batch_size; l++) {
            double E = instantaneous_energy(l);
            energy_samples[l].push_back(E);
        }        

        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start;
        std::cout << "Elapsed time: " << elapsed3.count() << " seconds\n";
        std::cout << "Finished forward "<< std::endl;        

        // Output configurations
        auto write_start = std::chrono::high_resolution_clock::now();

        std::vector<std::thread> write_threads;
        // std::mutex counter_mutex;

        for (int thread_num = 0; thread_num < batch_size; thread_num++) {
                 write_threads.emplace_back([&, thread_num]() { // create a new thread and store it in the vector write_threads
                        if (polymer[thread_num][0][2] < 1000) {
                        std::string filename = "configuration_batch" + std::to_string(batch) + "_thread" + std::to_string(thread_num) + ".txt";    
                        std::ofstream out(base_path + "final_confs/final_configuration_" + filename);
                        if (!out.is_open()) {
                        std::cerr << "ERROR: could not open file: " << base_path + "final_confs/final_configuration_" + filename << "\n";
                        return;
                        }
                
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

         for (int l = 0; l < number_of_threads; ++l) {
            double acceptance_rate = static_cast<double>(accepted_moves[l]) / static_cast<double>(attempted_moves[l]);
            std::cout << "T = " << Ts[l] << " acceptance = " << acceptance_rate << std::endl;
        }

        if (replica_exchange) {
            print_swap_stats(std::cout);
            replica_out << batch;
            for (int l = 0; l < number_of_threads; ++l) replica_out << '\t' << replica_at[l];
            replica_out << std::endl;
        }
    }

    if (replica_exchange) {
        std::ofstream swap_out(base_path + "swap_statistics.txt");
        print_swap_stats(swap_out);
        std::ofstream flow_out(base_path + "flow_statistics.txt");
        print_flow_stats(flow_out);
    }
    trace_out.close();

    std::ofstream energy_out(base_path + "energy_statistics.txt");

    energy_out << "# T\tbeta\tmean_E\tvariance_E\n";

    for (int l = 0; l < number_of_threads; ++l) {

        EnergyStats stats = compute_energy_stats(energy_samples[l]);

        energy_out
            << Ts[l] << '\t'
            << betas[l] << '\t'
            << stats.mean << '\t'
            << stats.variance << '\n';

        std::cout
            << "T = " << Ts[l]
            << "   <E> = " << stats.mean
            << "   Var(E) = " << stats.variance
            << std::endl;
    }

    energy_out.close();    

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Completed all batches." << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;

}




