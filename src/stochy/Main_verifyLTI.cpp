/*
 * Stochy main file
 *
 *  Created on: 1 Feb 2020
 *      Author: Abraham
 *      Goal  : Analyze the LTI system x_{k+1} = 0.8 x_k + w_k where w_k is a
 *              standard Gaussian with zero mean and covariance 0.05 I_d. Here,
 *              the dimension of the system is d. The safety objective is to
 *              remain within the safe set [-1,1]^d for 10 time steps.
 *
 *              The data generated from this code is used for comparison with
 *              SReachTools.
 */

#include "taskExec.h"
#include "time.h"
#include <iostream>
#include <nlopt.hpp>
void run_stochy_on_LTI(int d, float scaling_a, float cov, float grid_step_size, 
    float safety_set_bound_max, std::string additional_filename_string);
std::string create_additional_filename_string(int prob_case, int d, float cov, 
    float grid_step_size);

int main(int argc, char **argv) {

  std::cout << " _______  _______  _______  _______  __   __  __   __ "
            << std::endl;
  std::cout << "|       ||       ||       ||       ||  | |  ||  | |  |"
            << std::endl;
  std::cout << "|  _____||_     _||   _   ||       ||  |_|  ||  |_|  |"
            << std::endl;
  std::cout << "| |_____   |   |  |  | |  ||       ||       ||       |"
            << std::endl;
  std::cout << "|_____  |  |   |  |  |_|  ||      _||       ||_     _|"
            << std::endl;
  std::cout << " _____| |  |   |  |       ||     |_ |   _   |  |   |  "
            << std::endl;
  std::cout << "|_______|  |___|  |_______||_______||__| |__|  |___| "
            << std::endl;
  std::cout << std::endl;
  std::cout << " Welcome!  Copyright (C) 2018  natchi92 " << std::endl;
  std::cout << std::endl;

  // ------------------------- Case study 3 - Scaling in dimensions
  // -----------------------------------
  std::cout << "------------ Performing verification of LTI system -----------"
            << std::endl;
  std::cout << std::endl;

  // 1  - verification of 2D system with 0.05 covariance, 
  // 2  - verification of 2D system with changing variance and changing step size, 
  // 3X - scalability evaluation of stochy | X \in {1, 2, 3}
  int prob_case = 2;              
  float scaling_a = 0.8;
  float safety_set_bound_max = 1; // For Case 1 and 3 only

  // Variable definitions that are not going to be overwritten, but added here
  // for global scope

  // Case 1
  int case1_grid_step_size_index_max = 4;
  float case1_grid_step_size_vector [case1_grid_step_size_index_max] 
    = {0.05, 0.15, 0.2, 1};

  // Case 2
  float case2_safety_set_bound_max = 3;

  int case2_grid_step_size_index_max = 4;
  float case2_grid_step_size_vector [case2_grid_step_size_index_max] 
    = {1, 0.2, 0.15, 0.05};
  int cov_index_max = 4;
  float cov_vector [cov_index_max] = {0.05, 0.1, 0.2, 0.5};

  // Case 3
  // Define variables outside switch case for scope --- All of these will
  // be overwritten
  int d = 0;              
  float cov = 0;           
  float grid_step_size = 0;
  int grid_step_size_index = 0;
  int cov_index = 0;
  int d_max = 0;
  float t_est = 0;
  int case_int = 0;
  std::string additional_filename_string = "";

  switch(prob_case)
  {
      case 1:                   // Verification of a 2D LTI system
        d = 2;
        cov = 0.05;
        t_est = 0.5;

        // Get current dimension
        std::cout << "Case 1: Estimated time to complete is " 
            << t_est << "??? hour(s)."<< std::endl << std::endl;

        for(grid_step_size_index = 0; 
            grid_step_size_index < case1_grid_step_size_index_max; 
            grid_step_size_index ++)
        {
            std::cout << "Case 1: Verification of a 2D LTI system" << std::endl;
            grid_step_size = case1_grid_step_size_vector[grid_step_size_index];
            additional_filename_string = create_additional_filename_string(
                prob_case, d, cov, grid_step_size);
            run_stochy_on_LTI(d, scaling_a, cov, grid_step_size, 
                safety_set_bound_max, additional_filename_string);
        }
        break;
      case 2:
        d = 2;
        cov = 0.05;               // To be overwritten
        grid_step_size = 0;       // To be overwritten
        t_est = 1;
        
        std::cout << "Case 2: Estimated time to complete is " 
            << t_est << "??? hour(s)."<< std::endl << std::endl;
        // Get current dimension
        for(grid_step_size_index = 0; 
            grid_step_size_index < case2_grid_step_size_index_max; 
            grid_step_size_index ++)
        {
            for(cov_index = 0; cov_index < cov_index_max; cov_index ++)
            {
                std::cout << "Case 2: Verification of a 2D LTI system" 
                    << std::endl;
                grid_step_size = case2_grid_step_size_vector[
                  grid_step_size_index];
                cov = cov_vector[cov_index];
                additional_filename_string = create_additional_filename_string(
                    prob_case, d, cov, grid_step_size);
                run_stochy_on_LTI(d, scaling_a, cov, grid_step_size, 
                    case2_safety_set_bound_max, additional_filename_string);
            }
        }
        break;
      case 31:                   // Verification of a n-D LTI system
      case 32:                   // Verification of a n-D LTI system
      case 33:                   // Verification of a n-D LTI system
        switch(prob_case)
        {
            case 31:
                d_max = 10;
                grid_step_size = 1;
                t_est = 2;
                break;
            case 32:
                d_max = 3;
                grid_step_size = 0.2;
                t_est = 1;
                break;
            case 33:
                d_max = 3;
                grid_step_size = 0.15;
                t_est = 1;
                break;
            default:
              std::cout << "Invalid prob_case provided" << std::endl;
        }
        case_int = prob_case % 10;
        std::cout << "Case 3." << case_int <<": Estimated time to complete is " 
            << t_est << " hour(s)."<< std::endl << std::endl;
        // Switch dimension d between d_max, d_max - 1, ..., 2
        for(d=d_max; d >= 2 ; d = d-1)
        {
            std::cout << "Case 3." << case_int 
              << " Scalability of n-D LTI system" << std::endl;
            additional_filename_string = create_additional_filename_string(
                prob_case, d, cov, grid_step_size);
            run_stochy_on_LTI(d, scaling_a, cov, grid_step_size, 
                safety_set_bound_max, additional_filename_string);
        }
        break;
      default:
        std::cout << "Invalid prob_case provided" << std::endl;
  }
}


void run_stochy_on_LTI(int d, float scaling_a, float cov, float grid_step_size, 
    float safety_set_bound_max, std::string additional_filename_string)
{
  std::cout << "Dimension      = " << d << std::endl;
  std::cout << "Variance       = " << cov << std::endl;
  std::cout << "Scaling a      = " << scaling_a << std::endl;
  std::cout << "Safety bound   = " << safety_set_bound_max << std::endl;
  std::cout << "Grid step size = " << grid_step_size << std::endl << std::endl;
  
  // Define the boundary for each dimension
  arma::mat bound(1, 2);
  bound(0, 0) = -safety_set_bound_max;
  bound(0, 1) =  safety_set_bound_max;

  // Define grid size for each dimension
  arma::mat grid = grid_step_size * arma::ones<arma::mat>(1, 1);

  // Define relative tolerance for each dimension
  arma::mat reft = arma::ones<arma::mat>(1, 1);

  // Concatenate matrices to correct dimension
  // Boundary = [-1, 1]^d
  // Grid = [1]^d
  // Reft = [1]^d
  for (unsigned con = 0; con < d- 1; con++) {
    arma::mat tempBound(1, 2);
    tempBound(0, 0) = -safety_set_bound_max;
    tempBound(0, 1) =  safety_set_bound_max;
    bound = join_vert(bound, tempBound);
    grid = join_horiz(grid, grid_step_size * arma::ones<arma::mat>(1, 1));
    reft = join_horiz(reft, arma::ones<arma::mat>(1, 1));
  }
  // Define time horizon
  // Infinite time horizon:  T = -1
  // Finite time horizon: T = k (where k is an integer value)
  int T = 10;

  // Task definition (1 = simulator, 2 = faust^2, 3 = imdp)
  int lb  =imdp;

  // Property type
  // (1 = verify safety, 2= verify reach-avoid, 3 = safety synthesis, 4 = reach-avoid synthesis)
  int p = verify_safety;

  // task specification
  taskSpec_t cs3Spec(lb, T, p, bound, grid, reft);

  // Define model dynamics
  arma::mat Am3 = scaling_a * arma::eye(d,d);
  arma::mat Gm3 = cov * arma::eye(d,d);

  ssmodels_t model3(Am3, Gm3);

  std::vector<ssmodels_t> models3 = {model3};
  shs_t<arma::mat, int> cs3SHS(d, models3);

  // Combine
  inputSpec_t<arma::mat, int> cs3Input(cs3SHS, cs3Spec);

  // Perform task
  performTask(cs3Input, additional_filename_string);
  std::cout
      << "-----------------------------------------------------------------"
      << std::endl;
  std::cout << "Completed verification task." << std::endl;
  std::cout
      << "-----------------------------------------------------------------"
      << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

std::string create_additional_filename_string(int prob_case, int d, float cov, 
    float grid_step_size)
{
  return "Case_"+ std::to_string(prob_case) + "_dim_" + std::to_string(d) 
    + "_cov_" + std::to_string(int(cov * 1000)) + "_grid_step_size_" 
    + std::to_string(int(grid_step_size * 1000)) + "_";
}
