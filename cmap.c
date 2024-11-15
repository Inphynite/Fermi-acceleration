#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>

// C code for the calculation of (final or every step) phases and velocities of each particle after (just 1 or n) collisions
// Data written on a txt file, to be read and plotted in a Python script



// ~1h runtime (~90 min for sawtooth) for 40k particles @ 600k collisions (MinGW, AMD Ryzen 7 7840HS)

// Function that generates a random number under a normal distribution (around 0)
static double gauss(){
  	double x = (double) rand() / RAND_MAX;
	double y = (double) rand() / RAND_MAX;
	double z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
  	return z;
}

// Triangle wave function
double tri(double x){
    double y = (2 / M_PI) * asin(sin(x));
    return y;
}

// Sawtooth wave function
double saw(double x){
    double y = 2 * ((x / (2 * M_PI)) - floor(0.5 + (x / (2 * M_PI))));
    return y;
}

// Function that, given the initial velocity and phase, returns the velocity and phase of the particle after n collisions (sine wave function)
double *phase_velocity(double u0, double psi0, int n){
    double u = u0; // Initial velocity
    double psi = psi0; // Initial phase
    double M = 100.; // M parameter value (M = l/2*pi*a) with (l: wall distance) and (a: oscillating wall amplitude)
    static double res[2]; // 2 element array containing final u and psi
    
    for (int i = 0; i < n; i++){
        
        u = fabs(u + sin(psi));
        psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI)); 
    }  

    // Conditions to avoid writing NaN values in the file (division by 0)
    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        // Conditions to center stability zones in the graph
        if (psi < M_PI){
            res[1] = psi + M_PI;
        }
        else {
            res[1] = psi - M_PI;
        }
        
        return res;
    }
}

// Function that, given the initial velocity and phase, returns the velocity and phase of the particle after 1 collision (sine wave function)
double *phase_velocity_single(double u0, double psi0){
    double u = u0; // Initial velocity
    double psi = psi0; // Initial phase
    double M = 180.; // M parameter value (M = l/2*pi*a) with (l: wall distance) and (a: oscillating wall amplitude)
    static double res[2]; // 2 element array containing final u and psi
            
    u = fabs(u + sin(psi));
    psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI)); 
     

    // Conditions to avoid writing NaN values in the file (division by 0)
    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        
        res[1] = psi;
        
        return res;
    }
}

// Function that, given the initial velocity and phase, returns the velocity and phase of the particle after 1 collision (sine wave function) with noise velocity term
double *phase_velocity_stoch_single(double u0, double psi0, double u_st){
    double u = u0; // Initial velocity
    double psi = psi0; // Initial phase
    double M = 10.; // M parameter value (M = l/2*pi*a) with (l: wall distance) and (a: oscillating wall amplitude)
    static double res[2]; // 2 element array containing final u and psi
            
    u = fabs(u + sin(psi) + u_st);
    psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI)); 
     

    // Conditions to avoid writing NaN values in the file (division by 0)
    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        // Conditions to center stability zones in the graph
        res[1] = psi;
        
        return res;
    }
}

// Function that, given the initial velocity and phase, returns the velocity and phase of the particle after 1 collision (sawtooth wave function)
double *phase_velocity_saw_single(double u0, double psi0){
    double u = u0; // Initial velocity
    double psi = psi0; // Initial phase
    double M = 10.; // M parameter value (M = l/2*pi*a) with (l: wall distance) and (a: oscillating wall amplitude)
    static double res[2]; // 2 element array containing final u and psi
            
    u = fabs(u + saw(psi));
    psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI)); 
     

    // Conditions to avoid writing NaN values in the file (division by 0)
    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        // Conditions to center stability zones in the graph
        res[1] = psi;
        
        return res;
    }
}

// Phase velocity function with additional stochastic velocity term (sine wave function)
double *phase_velocity_stoch(double u0, double psi0, double u_st, int n){
    double u = u0;
    double psi = psi0;
    double M = 100.; 
    static double res[2]; 
    
    for (int i = 0; i < n; i++){
        
        u = fabs(u + sin(psi) + u_st);
        psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI));
    }  

    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        // Conditions to center stability zones in the graph
        if (psi < M_PI){
            res[1] = psi + M_PI;
        }
        else {
            res[1] = psi - M_PI;
        }
        
        return res;
    }
}

// Phase velocity function for a sawtooth wave function
double *phase_velocity_saw(double u0, double psi0, int n){
    double u = u0;
    double psi = psi0;
    double M = 100.; 
    static double res[2]; 
    
    for (int i = 0; i < n; i++){
        
        u = fabs(u + saw(psi));
        psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI));
    }  
    // Conditions to avoid writing NaN values in the file (division by 0) - causes issues when plotting
    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        // Conditions to center stability zones in the graph
        if (psi < M_PI){
            res[1] = psi + M_PI;
        }
        else {
            res[1] = psi - M_PI;
        }
        
        return res;
    }
}

// Phase velocity function with sawtooth wave function and additional noise velocity term
double *phase_velocity_saw_stoch(double u0, double psi0, double u_st, int n){
    double u = u0;
    double psi = psi0;
    double M = 100.; 
    static double res[2]; 
    
    for (int i = 0; i < n; i++){
        
        u = fabs(u + saw(psi) + u_st);
        psi = fmod(psi + ((2 * M_PI * M) / u), (2 * M_PI));
    }  

    if (isnan(u) || isnan(psi) || isinf(u)){
        res[0] = 0;
        res[1] = 0;
        return res;
    }
    else {
        res[0] = u;
        // Conditions to center stability zones in the graph
        if (psi < M_PI){
            res[1] = psi + M_PI;
        }
        else {
            res[1] = psi - M_PI;
        }
        
        return res;
    }
}

int main(){
    srand(time(NULL)); // Random number generation seed
    int n = 10000000; // Number of collisions

    // Adjustable variables, only useful for final maps and multiple trajectories
    int npsi = 50; // Number of initial phase points (i.e. Number of particles)
    int nv = 50; // Number of initial velocity points (gaussian)

    double sigma = pow(10, -3); // Gaussian sigma (standard deviation)

    // Array containing npsi initial phase points in [0, 2pi], therefore representing npsi particles
    double *psi_arr = calloc(npsi, sizeof(double));
    for (int k = 0; k < npsi; k++){
        *(psi_arr + k) = (k * (2 * M_PI / npsi));
    }
    // Array containing nv initial velocities (for nv different particles with the same initial phase)
    double *v_arr = calloc(nv, sizeof(double));
    for (int i = 0; i < nv; i++){
        *(v_arr + i) = gauss() * sigma + 1; // + 1 sets the average value as 1
    }
    
    // File to be written on
    FILE *data;
    data = fopen("data.txt", "w");
    fprintf(data, "%s,%s\n", "vit", "ph");
    
    
    // 2-element list containing initial conditions (adjustable)
    double *p0 = phase_velocity_single(1, M_PI / 2);
    // Main for-loop (for 1 single trajectory) (adjustable for each function to be called)
    for (int l = 0; l < n; l++){
        fprintf(data, "%lf,%lf\n", *p0, *(p0 + 1));
        //p0 = phase_velocity_stoch_single(*p0, *(p0 + 1), gauss() * sigma); // Function call for a map with a noise velocity term
        //p0 = phase_velocity_saw_single(*p0, *(p0 + 1));
        p0 = phase_velocity_single(*p0, *(p0 + 1));
            
    }
    
    fclose(data);
    /*
    // For velocity distribution (writing and displaying only final phase-velocities of  several particles)
    for (int j = 0; j < nv; j++){
        //double *p0 = phase_velocity_saw(*(v_arr + j), M_PI / 2, n); 
        //double *p0 = phase_velocity_saw(1, *(psi_arr + j), n);
        
        fprintf(data, "%lf,%lf\n", *p0, *(p0 + 1)); 
    }
    fclose(data);
    */
    free(psi_arr);
    free(v_arr);
    

}