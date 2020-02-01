/*
 * Stochy main file
 *
 *  Created on: 1 Feb 2020
 *      Author: Abraham
 *      Goal  : Analyze the LTI system x_{k+1} = 0.8 x_k + w_k where w_k is a
 *              standard Gaussian with zero mean and covariance 0.2 I_d. Here,
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
void run_stochy_on_LTI(int d, float cov, float grid_step_size);

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

  int prob_case = 22;           // ( 1 - verification, 
                               //   21 - scalability with grid step size 1,
                               //   22 - scalability with grid step size 0.1)
  float cov = 0.05;
  int d = 0;
  int d_max = 0;
  float grid_step_size = 0;

  switch(prob_case)
  {
      case 1:                   // Verification of a 2D LTI system
        d = 2;
        
        // Get current dimension
        std::cout << "Case 1: Verification of a 2D LTI system" << std::endl;
        std::cout << "This experiment took about 4+60+1100 ~ 20 minutes" 
            << std::endl << std::endl;

        // Switch grid step size between 0.2, 0.1, 0.05
        for(grid_step_size=0.2; grid_step_size > 0.049 ; 
                grid_step_size = grid_step_size/2)
        {
            std::cout << "Case 1: Verification of a 2D LTI system" << std::endl;
            std::cout << "Dimension = " << d << std::endl;
            std::cout << "Variance = " << cov << std::endl;
            std::cout << "Grid step size " << grid_step_size << std::endl 
                << std::endl;
            run_stochy_on_LTI(d, cov, grid_step_size);
        }
        break;
      case 21:                   // Verification of a n-D LTI system
      case 22:                   // Verification of a n-D LTI system
        if(prob_case == 21)
        {
            grid_step_size = 1;
            d_max = 10;
            std::cout << "Case 2: Scalability of n-D LTI system" << std::endl;
            std::cout << "This experiment took about 2 hours" 
                << std::endl << std::endl;
        }
        else{
            grid_step_size = 0.15;
            d_max = 3;
            std::cout << "Case 2: Scalability of n-D LTI system" << std::endl;
            std::cout << "This experiment took about 1 hour" 
                << std::endl << std::endl;
        }
        // Switch dimension d between 2 to 10
        for(d=3; d >= 2 ; d = d-1)
        {
            std::cout << "Case 2: Scalability of n-D LTI system" << std::endl;
            std::cout << "Dimension = " << d << std::endl;
            std::cout << "Variance = " << cov << std::endl;
            std::cout << "Grid step size " << grid_step_size << std::endl 
                << std::endl;
            run_stochy_on_LTI(d, cov, grid_step_size);
        }
        break;
      default:
        std::cout << "Invalid prob_case provided" << std::endl;
  }
}


void run_stochy_on_LTI(int d, float cov, float grid_step_size)
{
  // Define the boundary for each dimesnion
  arma::mat bound(1, 2);
  bound(0, 0) = -1;
  bound(0, 1) = 1;

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
    tempBound(0, 0) = -1;
    tempBound(0, 1) = 1;
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
  arma::mat Am3 = 0.8 * arma::eye(d,d);
  arma::mat Gm3 = cov * arma::eye(d,d);

  ssmodels_t model3(Am3, Gm3);

  std::vector<ssmodels_t> models3 = {model3};
  shs_t<arma::mat, int> cs3SHS(d, models3);

  // Combine
  inputSpec_t<arma::mat, int> cs3Input(cs3SHS, cs3Spec);

  // Perform task
  performTask(cs3Input);
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
