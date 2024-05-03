#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

struct Atom {
    char element;
    double x, y, z;
};

// Function to read atomic coordinates from a file
std::vector<Atom> readCoordinates(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Atom> atoms;
    char element;
    double x, y, z;
    while (file >> element >> x >> y >> z) {
        Atom atom;
        atom.element = element;
        atom.x = x;
        atom.y = y;
        atom.z = z;
        atoms.push_back(atom);
    }

    file.close();
    return atoms;
}

// Function to solve Poisson equation using Successive Over-Relaxation (SOR) method in 3D
void solvePoisson3D_SOR(std::vector<std::vector<std::vector<double> > >& phi, const std::vector<Atom>& atoms, double dx, double dy, double dz, int max_iter, double tolerance, double omega) {
    int Nx = phi.size();
    int Ny = phi[0].size();
    int Nz = phi[0][0].size();

    // // Set up charge density (consider default charge for other elements)
    // for (const auto& atom : atoms) {
    //     int i = std::floor(atom.x / dx);  // Ensure index within grid boundaries (using floor)
    //     int j = std::floor(atom.y / dy);
    //     int k = std::floor(atom.z / dz);
    //     if (i >= 0 && i < Nx && j >= 0 && j < Ny && k >= 0 && k < Nz) {
    //         phi[i][j][k] = atom.element == 'O' ? -8.0 : 1.0; // Adjust charge for Oxygen
    //     }
    // }

    // Set boundary conditions (Dirichlet: potential = 0 at the boundaries)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            phi[i][j][0] = phi[i][j][Nz - 1] = 0.0; // Set potential to 0 at z boundaries
            phi[i][0][j] = phi[i][Ny - 1][j] = 0.0; // Set potential to 0 at y boundaries
        }
    }
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            phi[0][j][k] = phi[Nx - 1][j][k] = 0.0; // Set potential to 0 at x boundaries
        }
    }

    // Timing the iterations
    std::vector<std::chrono::microseconds> iteration_times;

    // Open the output file
    std::ofstream outputFile("poisson_SOR_3D.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Main iterative loop
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_error = 0.0;

        // Start timing the iteration
        auto start_time = std::chrono::steady_clock::now();

        // Iterate over interior points
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int k = 1; k < Nz - 1; ++k) {
                    double old_phi = phi[i][j][k];
                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega * (
                        (phi[i - 1][j][k] + phi[i + 1][j][k]) / (dx * dx) +
                        (phi[i][j - 1][k] + phi[i][j + 1][k]) / (dy * dy) +
                        (phi[i][j][k - 1] + phi[i][j][k + 1]) / (dz * dz)
                    ) / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
                    double error = std::abs(phi[i][j][k] - old_phi);
                    if (error > max_error) {
                        max_error = error;
                    }
                }
            }
        }

        // End timing the iteration
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        iteration_times.push_back(duration);

        // Check convergence
        if (max_error < tolerance) {
            outputFile << "Converged after " << iter + 1 << " iterations." << std::endl;

            // Output time for each iteration
            for (int i = 0; i < iteration_times.size(); ++i) {
                outputFile << "Iteration " << i + 1 << " time: " << iteration_times[i].count() << " microseconds" << std::endl;
            }

            outputFile.close();
            return;
        }
    }

    outputFile << "Warning: Did not converge after " << max_iter << " iterations." << std::endl;
    outputFile.close();
}

int main() {
    // Parameters
    const int Nx = 100;          // Number of grid points in x-direction
    const int Ny = 100;          // Number of grid points in y-direction
    const int Nz = 100;          // Number of grid points in z-direction
    const double Lx = 10.0;       // Length of the domain in x-direction
    const double Ly = 1-.0;       // Length of the domain in y-direction
    const double Lz = 10.0;       // Length of the domain in z-direction
    const double dx = Lx / (Nx - 1); // Grid spacing in x-direction
    const double dy = Ly / (Ny - 1); // Grid spacing in y-direction
    const double dz = Lz / (Nz - 1); // Grid spacing in z-direction
    const int max_iter = 10000;   // Maximum number of iterations
    const double tolerance = 1e-8; // Convergence criterion
    const double omega = 1.5;     // Relaxation parameter (optimal for convergence)

    // Read atomic coordinates from file
    std::vector<Atom> atoms = readCoordinates("H2O.txt");

    // Initialize arrays
    std::vector<std::vector<std::vector<double> > > phi(Nx, std::vector<std::vector<double> >(Ny, std::vector<double>(Nz, 0.0))); // Electrostatic potential

    // Solve Poisson equation using SOR method
    solvePoisson3D_SOR(phi, atoms, dx, dy, dz, max_iter, tolerance, omega);

    return 0;
}
