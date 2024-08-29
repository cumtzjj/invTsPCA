import numpy as np
import multiprocessing as mp

class Particle:
    def __init__(self, bounds):
        self.position = np.array([np.random.uniform(low, high) for low, high in bounds])
        self.velocity = np.array([np.random.uniform(-abs(high-low), abs(high-low)) for low, high in bounds])
        self.best_position = np.copy(self.position)
        self.best_value = float('inf')

def rosenbrock_function(x):
    return np.sum(100.0 * (x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)

def evaluate_particle(particle, objective_function):
    current_value = objective_function(particle.position)
    if current_value < particle.best_value:
        particle.best_value = current_value
        particle.best_position = np.copy(particle.position)
    return particle.position, particle.best_position, particle.best_value

def update_particle(particle, global_best_position, inertia, cognitive, social, bounds):
    r1, r2 = np.random.rand(len(particle.position)), np.random.rand(len(particle.position))
    cognitive_velocity = cognitive * r1 * (particle.best_position - particle.position)
    social_velocity = social * r2 * (global_best_position - particle.position)
    particle.velocity = inertia * particle.velocity + cognitive_velocity + social_velocity
    particle.position += particle.velocity
    
    # Apply boundary constraints
    for i in range(len(particle.position)):
        if particle.position[i] < bounds[i][0]:
            particle.position[i] = bounds[i][0]
            particle.velocity[i] = 0
        elif particle.position[i] > bounds[i][1]:
            particle.position[i] = bounds[i][1]
            particle.velocity[i] = 0

def pso(objective_function, bounds, num_particles, max_iter, inertia_start=0.9, inertia_end=0.4, cognitive=1.5, social=1.5, patience=100, min_best_increase=0.01):
    particles = [Particle(bounds) for _ in range(num_particles)]
    global_best_position = np.copy(particles[0].position)
    global_best_value = float('inf')
    
    best_value_history = []

    pool = mp.Pool(mp.cpu_count())

    for iteration in range(max_iter):
        # Linearly decreasing inertia
        inertia = inertia_start - (inertia_start - inertia_end) * (iteration / max_iter)
        
        results = pool.starmap(evaluate_particle, [(particle, objective_function) for particle in particles])

        for i, (position, best_position, best_value) in enumerate(results):
            particles[i].position = position
            particles[i].best_position = best_position
            particles[i].best_value = best_value

            if best_value < global_best_value:
                global_best_value = best_value
                global_best_position = best_position

        # Update the best value history
        best_value_history.append(global_best_value)

        # Check for early stopping
        if len(best_value_history) > patience:
            if abs(best_value_history[-patience] - global_best_value) < min_best_increase:
                print(f"No improvement after {patience} iterations. Stopping early.")
                break

        for particle in particles:
            update_particle(particle, global_best_position, inertia, cognitive, social, bounds)

        print(f"Iteration {iteration+1}/{max_iter}, Best Value: {global_best_value}")

    pool.close()
    pool.join()

    return global_best_position, best_value_history

# Example usage
if __name__ == "__main__":
    # Define the bounds for a 10-dimensional search space
    #bounds = [(-5, 5) for _ in range(10)]  # 10-dimensional space with boundaries [-5, 5]
    lb=[-5,-5,-5,-5]
    ub=[5,5,5,5]
    bounds=[(i,j) for i,j in zip(lb,ub)]
    print(bounds)
    num_particles = 30
    max_iter = 1000

    best_position, best_value_history = pso(rosenbrock_function, bounds, num_particles, max_iter, min_best_increase=0.01)
    print(f"Best Position: {best_position}")
    print(f"Best Value: {best_value_history}")
