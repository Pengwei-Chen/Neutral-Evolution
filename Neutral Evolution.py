import json
import random
import numpy
import matplotlib.pyplot as plt
import math


def calculate_distance(location_1, location_2):
    x1, y1 = location_1
    x2, y2 = location_2
    return math.sqrt(math.pow(x1 - x2, 2) + math.pow(y1 - y2, 2))

# Individual
class Individual:
    def __init__(self, age, child_bearing_age, death_rate, gender, genotype, lifespan, location,
                    map_size, mated, mating_radius, mobility, mutation_rate, reproductive_cycle, resource_pressure):
        self.age = age
        self.child_bearing_age = child_bearing_age
        self.genotype = genotype
        self.death_rate = death_rate
        self.gender = gender
        self.lifespan = lifespan
        self.location = location
        self.map_size = map_size
        self.mated = mated
        self.mating_radius = mating_radius
        self.mobility = mobility
        self.mutation_rate = mutation_rate
        self.reproductive_cycle = reproductive_cycle
        self.resource_pressure = resource_pressure

    def move(self):
        x = self.location[0] + random.randint(-self.mobility, self.mobility)
        y = self.location[1] + random.randint(-self.mobility, self.mobility)
        while calculate_distance((x, y), self.location) > self.mobility:
            x = self.location[0] + random.randint(-self.mobility, self.mobility)
            y = self.location[1] + random.randint(-self.mobility, self.mobility)
        # Moving out of the map would be considered as emigration
        if x < 1 or x > self.map_size[0] or y < 1 or y > self.map_size[1]:
            return None
        self.location = (x, y)
        return self
    
    def grow(self):
        self.age += 1
        if self.age == self.child_bearing_age:
            self.mated = self.child_bearing_age
        if self.age > self.child_bearing_age:
            self.mated += 1 
        return self
    
    def mutate(self):
        for i in self.genotype:
            genotype = list(self.genotype)
            if random.random() < self.mutation_rate:
                genotype[i] = not genotype[i]
            self.genotype = tuple(genotype)
        return self

    def die(self):
        if self.age >= self.lifespan or random.random() < self.death_rate or random.random() < self.resource_pressure:
            return None
        return self

# To generate true age accoring to age structure
def generate_age(life_span, death_rate):
    age = 0
    for i in range(life_span):
        if random.random() < 1 - death_rate:
            age +=  1
        else:
            return age
    return generate_age(life_span, death_rate)

# Initialize individuals in the population
def initialize(config):
    population = []
    for i in range(config["Initial population size"]):
        age = generate_age(config["Life span"], config["Death rate"])
        location_shift = (random.randint(-config["Initial distribution radius"], config["Initial distribution radius"]),
                            random.randint(-config["Initial distribution radius"], config["Initial distribution radius"]))
        while calculate_distance(location_shift, (0, 0)) > config["Initial distribution radius"]:
            location_shift = (random.randint(-config["Initial distribution radius"], config["Initial distribution radius"]),
                    random.randint(-config["Initial distribution radius"], config["Initial distribution radius"]))
        population.append(Individual(age = age,
                            child_bearing_age = config["Child-bearing age"],
                            death_rate = config["Death rate"],
                            gender = random.randint(0, 1),
                            genotype = (random.random() < config["Initial gene A percent"], random.random() < config["Initial gene A percent"]),
                            lifespan = config["Life span"],
                            location = (int(config["Map X"] / 2) + location_shift[0], int(config["Map Y"] / 2) + location_shift[1]),
                            map_size = (config["Map X"], config["Map Y"]),
                            mated = age % config["Child-bearing age"],
                            mating_radius = config["Mating radius"],
                            mobility = config["Mobility"],
                            mutation_rate = config["Mutation rate"],
                            reproductive_cycle = config["Reproductive cycle"],
                            resource_pressure = 0))
    return population

# Update the population
def develop(population, config):
    # Move, grow, and mutate
    for i in range(len(population)):
        population[i] = population[i].move()
    population = list(filter(None, population))
    for i in range(len(population)):
        population[i] = population[i].grow().mutate()
    # Mating
    new_individuals = []
    for i in range(len(population) - 1):
        if population[i].mated >= population[i].reproductive_cycle:
            distance = []
            for j in range(i, len(population)):
                if population[i].gender != population[j].gender and population[j].mated >= population[j].reproductive_cycle:
                    distance.append((j, calculate_distance(population[i].location, population[j].location)))
            if distance != []:
                distance = sorted(distance, key = lambda x : x[1])
                population[i].mated = 0
                population[distance[0][0]].mated = 0
                new_individuals.append(Individual(age = 0,
                            child_bearing_age = config["Child-bearing age"],
                            death_rate = config["Death rate"],
                            gender = random.randint(0, 1),
                            genotype = (population[i].genotype[random.randint(0, 1)], population[distance[0][0]].genotype[random.randint(0, 1)]),
                            lifespan = config["Life span"],
                            location = (int((population[i].location[0] + population[distance[0][0]].location[0]) / 2), 
                                            int((population[i].location[1] + population[distance[0][0]].location[1]) / 2)),
                            map_size = (config["Map X"], config["Map Y"]),
                            mated = 0,
                            mating_radius = config["Mating radius"],
                            mobility = config["Mobility"],
                            mutation_rate = config["Mutation rate"],
                            reproductive_cycle = config["Reproductive cycle"],
                            resource_pressure = 0))
    
    # Resource competition
    distribution = numpy.zeros((config["Map X"], config["Map Y"]), dtype = int)
    for individual in population:
        for i in range(-3, 4):
            for j in range(-3, 4):
                if individual.location[0] - 1 + i >= 0 and individual.location[0] - 1 + i < config["Map X"] \
                        and individual.location[1] - 1 + j >= 0 and individual.location[1] - 1 + j < config["Map Y"]:
                    distribution[individual.location[0] - 1 + i, individual.location[1] - 1 + j] += 1
    for i in range(len(population)):
        population[i].resource_pressure = distribution[population[i].location[0] - 1, population[i].location[1] - 1] * 0.01 - 0.01
    
    # Resource shortage outside the initial range
    middle_point = (int(config["Map X"] / 2), int(config["Map Y"] / 2))
    for i in range(len(population)):
        distance = int(calculate_distance(population[i].location, middle_point))
        if distance > config["Initial distribution radius"]:
            population[i].resource_pressure += distance / config["Initial distribution radius"] * 0.3

    # Die because of aging, resource competition, or random events
    for i in range(len(population)):
        population[i] = population[i].die()
    population = list(filter(None, population))

    # Add newborns
    population = population + new_individuals
    return population

# Generate map for visualization
def generate_map(population, map_size):
    map = numpy.zeros((map_size[0], map_size[1] + 1), dtype = int)
    for i in range(4):
        map[i, map_size[1]] = i
    for individual in population:
        genotype = individual.genotype[0] + individual.genotype[1] + 1
        map[individual.location[0] - 1, individual.location[1] - 1] = genotype
    return map

# Generate statistics for visualization
def generate_statistics(population):
    statistics = [0, 0]
    for individual in population:
        for gene in individual.genotype:
            if gene:
                statistics[0] += 1
            else:
                statistics[1] += 1
    return statistics

# Main
# Read configuration
config_file = open("Config.json", "r")
config = json.loads(config_file.read())
config_file.close()

# Initialize population and develop
population = initialize(config)
map = generate_map(population, map_size = (config["Map X"], config["Map Y"]))
a, b = generate_statistics(population)
gene_a = [a / (a + b)]
gene_b = [b / (a + b)]
timeline = [0]
plt.ion()
figure = plt.figure()
figure.canvas.manager.set_window_title("Neutral Evolution")
ax1 = plt.subplot(1, 2, 1)
plt.title("Map")
ax2 = plt.subplot(1, 2, 2)
plt.title("Statistics")
plt.xlabel("Time")
plt.ylabel("Percentage")
ax2.set_ylim(ymin = 0, ymax = 1)
image = ax1.imshow(map)
for time_point in range(1, config["Time scale"] + 1):
    population = develop(population, config)
    if len(population) == 0:
        print("The population extincts.")
        # Visualize
        map = generate_map(population, map_size = (config["Map X"], config["Map Y"]))
        image.set_data(map)
        ax2.plot(timeline, gene_a, color = "red")
        ax2.plot(timeline, gene_b, color = "blue")
        figure.canvas.draw()
        figure.canvas.flush_events()
        input()
        exit()
    a, b = generate_statistics(population)
    gene_a.append(a / (a + b))
    gene_b.append(b / (a + b))
    timeline.append(time_point)
    if time_point % 10 == 0 or a == 0 or b == 0:
        # Visualize
        map = generate_map(population, map_size = (config["Map X"], config["Map Y"]))
        image.set_data(map)
        ax2.plot(timeline, gene_a, color = "red")
        ax2.plot(timeline, gene_b, color = "blue")
        figure.canvas.draw()
        figure.canvas.flush_events()
        if a == 0 or b == 0:
            input()
            exit()
input()