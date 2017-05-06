/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

//std::default_random_engine generator;

std::random_device rd_engine;
std::mt19937 generator(rd_engine());


std::normal_distribution<double> distribution(0,1);

/*
* Calculates the bivariate normal pdf of a point given a mean and std and assuming zero correlation
*/
double bivariate_normal(double x, double y, double mu_x, double mu_y, double sig_x, double sig_y) {
  return exp(-((x-mu_x)*(x-mu_x)/(2*sig_x*sig_x) + (y-mu_y)*(y-mu_y)/(2*sig_y*sig_y))) / (2.0*3.14159*sig_x*sig_y);
}




void ParticleFilter::init(double x, double y, double theta, double std_var[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */


	 num_particles = 200;
	 weights.resize(num_particles);
   std::default_random_engine gen;
   std::normal_distribution<double> N_x_init(0, std_var[0]);
   std::normal_distribution<double> N_y_init(0, std_var[1]);
   std::normal_distribution<double> N_theta_init(0, std_var[2]);



	 for (int i = 0; i < num_particles; i++) {
		 		Particle p_i; // Define i-th particle of type stuct Particle
				// Fill out different fields for particle
				p_i.id = i;
				p_i.x = x + N_x_init(gen);
				p_i.y = y + + N_y_init(gen);
				p_i.theta = theta + N_theta_init(gen);
        p_i.weight = 1.0f;
        weights[i] = 1.0f;
				particles.push_back(p_i);

			}
      is_initialized = true;



}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */

   std::default_random_engine gen;
   std::normal_distribution<double> N_x_init(0, std_pos[0]);
   std::normal_distribution<double> N_y_init(0, std_pos[1]);
   std::normal_distribution<double> N_theta_init(0, std_pos[2]);


	 for (int i = 0; i < num_particles; i++) {

		 double th = particles[i].theta;
		 double v = velocity;
		 double dth = yaw_rate;


	 	if (fabs(dth)>0.001){
				particles[i].x +=   v/dth*(sin(th+dth*delta_t) - sin(th))+ N_x_init(gen);
				particles[i].y +=  -v/dth*(cos(th+dth*delta_t) - cos(th))+ N_y_init(gen);
				particles[i].theta += dth*delta_t+ N_theta_init(gen);

			} else {
				particles[i].x += v*cos(th)*delta_t+ N_x_init(gen);
				particles[i].y += v*sin(th)*delta_t+ N_y_init(gen);
				particles[i].theta += dth*delta_t+ N_theta_init(gen);
		}

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


	double dist_po;
	for (int i=0;i<observations.size();i++){ // loop over observations
		double min_dist = 1000000;
		int closest_id = -1;
		for (int j=0;j<predicted.size();j++){

			dist_po = dist(observations[i].x, observations[i].y,
										predicted[j].x, predicted[j].y);
			if (dist_po<min_dist){
				min_dist = dist_po;
				closest_id = predicted[j].id;
			}
		}
		observations[i].id  = closest_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	//std::cout<<"Obs_size: "<<observations.size()<<"Map_size: "<<map_landmarks.landmark_list.size()<<"Range:"<<sensor_range<<std::endl;
	weights.clear();


	for (int i_p=0;i_p<particles.size();i_p++){


    // Convert observations to ground frame
		std::vector<LandmarkObs> observations_grnd;
    for (int i=0;i<observations.size();i++){
      LandmarkObs obsv_i;
      obsv_i.x = 0;
			obsv_i.x += observations[i].x * cos(particles[i_p].theta);
			obsv_i.x += -observations[i].y * sin(particles[i_p].theta);
			obsv_i.x += particles[i_p].x;


      obsv_i.y = 0;
      obsv_i.y += observations[i].x * sin(particles[i_p].theta);
			obsv_i.y += observations[i].y * cos(particles[i_p].theta);
			obsv_i.y += particles[i_p].y;

			obsv_i.id = -1; // Temporary ID.
      observations_grnd.push_back(obsv_i);
			}



    std::vector<LandmarkObs> predicted_meas;
		// Compute predicted measurements
		for (int i=0;i<map_landmarks.landmark_list.size();i++){
			double distance_particle_obs;
			distance_particle_obs = dist(particles[i_p].x,particles[i_p].y,
																	map_landmarks.landmark_list[i].x_f,
																	map_landmarks.landmark_list[i].y_f);
		 	if (distance_particle_obs <= sensor_range){
				LandmarkObs predicted_i;
				predicted_i.id = map_landmarks.landmark_list[i].id_i;
				predicted_i.x = map_landmarks.landmark_list[i].x_f;
				predicted_i.y = map_landmarks.landmark_list[i].y_f;

			 	predicted_meas.push_back(predicted_i);
		 	}
		 }





		 dataAssociation(predicted_meas, observations_grnd);


		 double prob = 1.0;
     double prob_i;



		 for (int i =0; i<predicted_meas.size(); i++){
			 int ind_min = -1;
			 double dist_min = 1000000;


			 for (int j =0; j<observations_grnd.size(); j++){
         //Use measurement closest to predicted in case of multiple measurements
         // assignend to the same observation
				 if (predicted_meas[i].id == observations_grnd[j].id ){
					 double check_dist = dist(predicted_meas[i].x,
                                  predicted_meas[i].y,
						 											observations_grnd[j].x,
                                  observations_grnd[j].y);
					if (check_dist<dist_min){
						ind_min = j;
						dist_min = check_dist;
					}
				 }
			 }
			 if (ind_min!=-1){
			 		prob_i = bivariate_normal(predicted_meas[i].x, predicted_meas[i].y,
																observations_grnd[ind_min].x, observations_grnd[ind_min].y,
																 std_landmark[0], std_landmark[1]);

         prob = prob*prob_i;
			 }


		 }
		 weights.push_back(prob);
		 particles[i_p].weight = prob;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	  std::random_device rd_wts;
	  std::mt19937 generator_wts(rd_wts());

    //std::default_random_engine generator_wts;

    // Creates a discrete distribution for weight.
    std::discrete_distribution<int> distribution_wts(weights.begin(), weights.end());
    std::vector<Particle> resampled_particles;
    // Resample
    for(int i=0;i<num_particles;i++){
				Particle particles_i = particles[distribution_wts(generator_wts)];
        resampled_particles.push_back(particles_i);
    }
    particles = resampled_particles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
