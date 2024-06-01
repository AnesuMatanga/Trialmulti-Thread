// System Headers
#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>

// Project Headers
#include "nbody.h"

//#define GRAPHICS
#ifdef GRAPHICS
	#include <SFML/Window.hpp>
	#include <SFML/Graphics.hpp>
#endif

// Number of particles
// #define SMALL
 #define LARGE
// #define MASSIVE

#if defined(SMALL)
	const int N = 1000;
#elif defined(LARGE)
	const int N = 5000;
#elif defined(MASSIVE)
	const int N = 10000;
#endif

// Constants
const double min2 = 2.0;
const double G = 1 * 10e-10;
const double dt = 0.01;
const int NO_STEPS = 500;

// Size of Window/Output image
const int width = 1920;
const int height = 1080;

// Bodies
body bodies[N];
//body *bodies = new body[N];

// Update Nbody Simulation
void update() {
	// Acceleration
	vec2 acc[N] = {};

	// Initiallise threads number
	int num_threads = omp_get_max_threads();

	// Thread local, local storage for local acc
	vec2 local_acc[num_threads][N];

	// Initialise local_acc arrays to zero
	#pragma omp parallel
	{
		int ID = omp_get_thread_num();
		for(int i = 0; i < N; ++i){
			// Initiallise
			local_acc[ID][i] = vec2(0, 0);
		}
	}

	// Clear Acceleration
	// Using pragma for to parallelise for
	/*#pragma omp parallel for
	for(int i = 0; i < N; ++i) {
		acc[i] = vec2(0,0);
	}*/

	// For each body
	// Schedule(dynamic) To avoid race conditions, load imbalance, will locallise acc and for loops
	#pragma omp parallel
	{
		int ID = omp_get_thread_num();

		#pragma omp for schedule(dynamic)
		for(int i = 0; i < N; ++i) {
			//vec2 local_acc[N] = {};
			/*for (int j = 0; j < N; ++j){
				local_acc[j] = vec2(0, 0);
			}*/

			// For each following body
			for(int j = i+1; j < N; ++j) {
				// No point computing if i and j are the same
				if(i != j){
					// Difference in position
					vec2 dx = bodies[i].pos - bodies[j].pos;

					// Normalised difference in position
					vec2 u = normalise(dx);

					// Calculate distance squared
					double d2 = length2(dx);

					// If greater than minimum distance
					if(d2 > min2) {
						// Smoothing factor for particles close to minimum distance
						double x = smoothstep(min2, 2 * min2, d2);

						// Force between bodies
						double f = -G*bodies[i].mass*bodies[j].mass / d2;

						// Add to acceleration
						//local_acc[ID][i] += (u * f / bodies[i].mass) * x;
						//local_acc[ID][j] -= (u * f / bodies[j].mass) * x;
						// Add to local acceleration
                    	vec2 ai = (u * f / bodies[i].mass) * x;
                    	vec2 aj = (u * f / bodies[j].mass) * x;

                    	
                    	local_acc[ID][i] += ai;
                    	
                    	local_acc[ID][j] -= aj;

					}
				}
			}
		}
	}

	// Combine local accelerations
	#pragma omp parallel for
	//acc[i] += local_acc[i];
	for (int j = 0; j < N; ++j){
		for (int t = 0; t < num_threads; ++t){
			acc[j] += local_acc[t][j];
		}
	}
		

	// For each body
	// Create threads to execute the following block of code & split loop iterations between threads
	#pragma omp parallel for 
		for(int i = 0; i < N; ++i) {
			// Update Position
			bodies[i].pos += bodies[i].vel * dt;

			// Update Velocity
			bodies[i].vel += acc[i] * dt;
		}
}

// Initialise NBody Simulation
void initialise() {
	// Create a central heavy body (sun)
	bodies[0] = body(width/2, height/2, 0, 0, 1e15, 5);

	// For each other body
	for(int i = 1; i < N; ++i) {
		// Pick a random radius, angle and calculate velocity
		double r = (uniform() + 0.1) * height/2;
		double theta = uniform() * 2 * M_PI;
		double v = sqrt(G * (bodies[0].mass + bodies[i].mass) / r);

		// Create orbiting body
		bodies[i] = body(width/2 + r * cos(theta), height/2 + r * sin(theta), -sin(theta) * v, cos(theta)*v, 1e9, 2);
	}
}

#ifdef GRAPHICS
	// Main Function - Graphical Display
	int main() {
		// Create Window
		sf::ContextSettings settings;
		settings.antialiasingLevel = 1;
		sf::RenderWindow window(sf::VideoMode(width, height), "NBody Simulator", sf::Style::Default, settings);

		// Initialise NBody Simulation
		initialise();
		int i = 0;
		// run the program as long as the window is open
		while (window.isOpen()) {
			// check all the window's events that were triggered since the last iteration of the loop
			sf::Event event;
			while (window.pollEvent(event)) {
				// "close requested" event: we close the window
				if (event.type == sf::Event::Closed) {
					window.close();
				}
			}

			if(i < NO_STEPS) {
				// Update NBody Simluation
				update();
				i++;
			}

			// Clear the window with black color
			window.clear(sf::Color::Black);

			// Render Objects
			for(int i = 0; i < N; ++i) {
				// Create Circle
				sf::CircleShape shape(bodies[i].radius);
				shape.setFillColor(sf::Color(255, 0, 0));
				shape.setPosition(bodies[i].pos.x, bodies[i].pos.y);
				
				// Draw Object
				window.draw(shape);
			}

			// Display Window
			window.display();
		}
	}
#else
	// Main Function - Benchmark
	int main() {

		// Set the number of threads to the maximum available according to hardware
		omp_set_num_threads(omp_get_max_threads());
		// omp_set_num_threads(6);

		// Initialise NBody Simulation
		initialise();

		// Get start time
		std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

		// Run Simulation
		for(int i = 0; i < NO_STEPS; i++) {
			// Update NBody Simluation
			update();
		}

		// Get end time
		std::chrono::system_clock::time_point end = std::chrono::system_clock::now();

		// Generate output image
		unsigned char *image = new unsigned char[width * height * 3];
		memset(image, 0, width * height * 3);

		// For each body
		for(int i = 0; i < N; ++i) {
			// Get Position
			vec2 p = bodies[i].pos;

			// Check particle is within bounds
			if(p.x >= 0 && p.x < width && p.y >= 0 && p.y < height) {
				// Add a red dot at body
				image[((((int)p.y * width) + (int)p.x) * 3)] = 255;
			}
		}

		// Write position data to file
		char data_file[200];
		sprintf(data_file, "output%i.dat", N);
		write_data(data_file, bodies, N);

		// Write image to file
		char image_file[200];
		sprintf(image_file, "output%i.png", N);
		write_image(image_file, bodies, N, width, height);

		// Check Results
		char reference_file[200];
		sprintf(reference_file, "reference%i.dat", N);
		calculate_maximum_difference(reference_file, bodies, N);

		// Time Taken
		std::cout << "Time Taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1000000.0 << std::endl;
	}
#endif