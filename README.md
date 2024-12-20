# Optimization Algorithms Comparison

This MATLAB project implements and compares four optimization algorithms on three different optimization problems.

## Algorithms Implemented

1. Particle Swarm Optimization (PSO)
2. Second-Order Vibrating PSO (SecVibratPSO)
3. Simulated Annealing (SA)
4. Genetic Algorithm (GA)

## Optimization Problems

1. Ackley Function (continuous optimization)
2. G06 Problem (constrained optimization)
3. 0-1 Knapsack Problem (combinatorial optimization)

## Project Structure

```
.
├── algorithms/
│   ├── PSO.m
│   ├── SecVibratPSO.m
│   ├── APSO.m
│   ├── SA.m
│   └── GA.m
├── problems/
│   ├── Ackley.m
│   ├── G06.m
│   └── Knapsack.m
├── utils/
│   ├── plotConvergence.m
│   ├── plotStatistics.m
│   └── displayResults.m
├── results/
├── main.m
└── README.md
```

## How to Run

1. Make sure all files are in their respective directories
2. Open MATLAB and navigate to the project directory
3. Run the main script:
   ```matlab
   main
   ```

## Results

The program will generate:
- Convergence plots for each problem (saved as PNG files)
- Statistical comparison plots (saved as PNG files)
- A text file with detailed results (optimization_summary.txt)

All results will be saved in the `results/` directory.

## Parameters

You can modify the following parameters in `main.m`:
- `max_iter`: Maximum number of iterations
- `pop_size`: Population size for PSO, SOPSO, and GA
- `runs`: Number of independent runs for statistical analysis

## Contact

For any questions or issues, please open an issue on the repository.
