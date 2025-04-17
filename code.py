# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation
import time
import pickle

class Cell():
    def __init__(self, row, column, size, is_cancer, args = None):
        self.is_cancer = is_cancer
        self.is_dying = False
        self.is_invading = False
        self.y_position = column
        self.x_position = row
        self.cell_cycle_stage = "G0"
        self.stage_duration = None
        self.stage_counter = 0
        self.oxygen_level = .3
        self.nutrient_level = .3
        self.chemo_level = 0.0
        self.neighbors = self.initialize_neighbors(row, column, size)
        self.num_mutation = 0
        self.damage = 0.0
        if is_cancer:
            number = int(np.random.normal(loc=5000, scale=200))
            self.lifespan = max(0, number)
            if args == None:
                self.mutation_probability = 0.02
                self.invasion_prob = 0.2
                self.growth_acceleration_factor = 2
                self.resistance_apoptosis = 0.2
                self.recovery_rate = 0.3
                self.age = np.random.randint(0, int(self.lifespan * .5))
        else:
            number = int(np.random.normal(loc=1000, scale=20))
            self.lifespan = max(0, number)
            if args == None:
                self.mutation_probability = 0.000001
                self.invasion_prob = 0.0 
                self.growth_acceleration_factor = 1
                self.resistance_apoptosis = 0.0 
                self.recovery_rate = 0.4
                self.age = np.random.randint(0,self.lifespan)
                
        if args != None:
            self.num_mutation = args[0]
            self.mutation_probability = args[1]
            self.invasion_prob = args[2]
            self.growth_acceleration_factor = args[3]
            self.resistance_apoptosis = args[4]
            self.recovery_rate = args[5]
            self.age = 0

    
    def update_cell(self, environment_nutrient, environment_oxygen, environment_chemo, is_empty_neighbor):
        if self.is_dying == False:
            self.update_internal_concentrations(environment_nutrient, environment_oxygen, environment_chemo)
            
            if self.oxygen_level < 0.05 or self.nutrient_level < 0.05:
                if np.random.random() > self.resistance_apoptosis:
                    self.is_dying = True
                    self.stage_counter = 0
                    
            if self.cell_cycle_stage == "S" or self.cell_cycle_stage == "M":
                damage_multiplier = 1.5
            elif self.cell_cycle_stage == "G1" or self.cell_cycle_stage == "G2":
                damage_multiplier = 1 
            else: 
                damage_multiplier = .5
            self.damage += self.chemo_level * damage_multiplier - self.recovery_rate 
            
            if self.damage < 0:
                self.damage = 0.0
            
            if self.damage > 1:
                self.damage = 1
            
            if self.age == self.lifespan:
                self.is_dying = True
                self.stage_counter = 0
            
            if np.random.random() <= self.mutation_probability:
                self.mutate_cell()
            
            if np.random.random() < self.damage - self.resistance_apoptosis:
                self.is_dying = True
                self.stage_counter = 0
                
            self.age += 1
            
            #11 hours in the G1 phase, 8 hours in the S phase, 4 hours in the G2 phase, and roughly 1 hour in the M (mitosis) phase
            if is_empty_neighbor and self.cell_cycle_stage == "G0":
                self.stage_duration = np.zeros(4)
                self.stage_duration[0] = int(np.random.normal(loc= 11 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[1] = int(np.random.normal(loc= 8 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[2] = int(np.random.normal(loc= 4 / self.growth_acceleration_factor, scale=2))
                self.stage_duration[3] = 1 
                self.cell_cycle_stage = "G1"
                self.stage_counter = 0
                is_ready_to_split = False
            
            elif self.cell_cycle_stage == "G1":
                if self.stage_counter == self.stage_duration[0]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "S"   
                else:
                    if self.oxygen_level > 0.1 and self.nutrient_level > 0.1:
                        self.stage_counter += 1
                is_ready_to_split = False
                
            elif self.cell_cycle_stage == "S":
                if self.stage_counter == self.stage_duration[1]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "G2"
                else:
                    self.stage_counter += 1 
                is_ready_to_split = False
                
            elif self.cell_cycle_stage == "G2":
                if self.stage_counter == self.stage_duration[2]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "M"
                else:
                    if self.oxygen_level > 0.1 and self.nutrient_level > 0.1:
                        self.stage_counter += 1
                is_ready_to_split = False
            elif self.cell_cycle_stage == "M":
                if self.stage_counter == self.stage_duration[3]:
                    self.stage_counter = 0
                    self.cell_cycle_stage = "G0"
                    if is_empty_neighbor or self.is_invading:
                        is_ready_to_split = True
                        self.is_invading = False
                    else: 
                        is_ready_to_split = False
                else:
                    self.stage_counter += 1
                    is_ready_to_split = False
                        
            elif np.random.random() < self.invasion_prob:
                self.stage_duration = np.zeros(4)
                self.stage_duration[0] = int(np.random.normal(loc= 11 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[1] = int(np.random.normal(loc= 8 / self.growth_acceleration_factor, scale=3))
                self.stage_duration[2] = int(np.random.normal(loc= 4 / self.growth_acceleration_factor, scale=2))
                self.stage_duration[3] = 1 
                self.cell_cycle_stage = "G1"
                self.is_invading = True
                self.stage_counter = 0
                is_ready_to_split = False
            
            else:
                is_ready_to_split = False
            
            return is_ready_to_split
        
        else:
            self.stage_counter += 1
            is_ready_to_split = False
            return is_ready_to_split
    
    def get_cell_attributes(self):
        return self.num_mutation, self.mutation_probability, self.invasion_prob, self.growth_acceleration_factor, self.resistance_apoptosis, self.recovery_rate, self.is_cancer
            
    def initialize_neighbors(self, row, column, size):
        neighbors = []
        neighbor_offsets = [
            (-1, -1), (-1, 0), (-1, 1),
            ( 0, -1),          ( 0, 1),
            ( 1, -1), ( 1, 0), ( 1, 1)]
        
        for dr, dc in neighbor_offsets:
            neighbor_row = row + dr
            neighbor_col = column + dc
            if 0 <= neighbor_row < size and 0 <= neighbor_col < size:
                neighbors.append((neighbor_row, neighbor_col))
        return neighbors
    
    def update_internal_concentrations(self, environment_nutrient, environment_oxygen, environment_chemo):
        if environment_oxygen > self.oxygen_level:
            self.oxygen_level += (environment_oxygen - self.oxygen_level) * environment_oxygen
        else:
            self.oxygen_level += (environment_oxygen - self.oxygen_level) * self.oxygen_level
        
        if self.oxygen_level > 1:
            self.oxygen_level = 1
        
        if environment_nutrient > self.nutrient_level:
            self.nutrient_level += (environment_nutrient - self.nutrient_level) * environment_nutrient
        else:
            self.nutrient_level += (environment_nutrient - self.nutrient_level) * self.nutrient_level
        
        if self.nutrient_level > 1:
            self.nutrient_level = 1
        
        if environment_chemo > self.chemo_level:
            self.chemo_level += (environment_chemo - self.chemo_level) * environment_chemo
        else:
            self.chemo_level += (environment_chemo - self.chemo_level) * self.chemo_level
        
        if self.chemo_level > 1:
            self.chemo_level = 1
            
        if self.chemo_level < 0.05:
            self.chemo_level = 0
    
    def mutate_cell(self):
        random_number = np.random.random()
        
        if random_number < .2:
            if self.mutation_probability < .99:
                self.mutation_probability += .01
        if random_number >= .2 and random_number < .4:
            if self.invasion_prob < .98:
                self.invasion_prob += .02
        if random_number >= .4 and random_number < .6:
            self.growth_acceleration_factor += .05
        if random_number >= .6 and random_number < .8:
            if self.resistance_apoptosis < .98:
                self.resistance_apoptosis += .02
        if random_number >= .8 and random_number < 1:
            if self.recovery_rate < .8:
                self.recovery_rate += .01
            
        self.damage += .05
    
        if self.damage > 1:
            self.damage = 1
        
        self.num_mutation += 1
            
        
class Extracellular_Matrix():
    def __init__(self, size):
        # [:,:,0] is Nutrients, [:,:,1] is oxygen, and [:,:,3] is chemo 
        self.size = size                              
        self.ecm = np.zeros((size, size, 3), dtype=float)
        self.vasculature = np.zeros((size, size))
        self.midline = size // 2
        
        self.vasculature[:, self.midline - 1] = 1
        
        nutrient_step_size = 0.02
        self.ecm[:,:self.midline, 0] = np.arange(1 - (self.midline) * nutrient_step_size, 1, nutrient_step_size)
        self.ecm[:, self.midline + 1:, 0] = np.arange(1 - nutrient_step_size,1 - (self.midline + 1) * nutrient_step_size, -nutrient_step_size)
        self.ecm[:, self.midline, 0] = 1 
         
        oxygen_step_size = 0.01
        self.ecm[:,:self.midline, 1] = np.arange(1 - (self.midline) * oxygen_step_size, 1, oxygen_step_size)
        self.ecm[:, self.midline + 1:, 1] = np.arange(1 - oxygen_step_size,1 - (self.midline + 1) * oxygen_step_size, -oxygen_step_size)
        self.ecm[:, self.midline, 1] = 1
    
    def set_chemo_concentration(self, max_concentration):
        chemo_step_size = 0.04
        self.ecm[:, :self.midline, 2] = np.arange(max_concentration - (self.midline) * chemo_step_size, max_concentration, chemo_step_size)
        self.ecm[:, self.midline + 1:, 2] = np.arange(max_concentration - chemo_step_size,max_concentration - (self.midline + 1) * chemo_step_size, - chemo_step_size)
        self.ecm[:, self.midline, 2] = max_concentration
        self.ecm[:,:,2] = np.where(self.ecm[:,:,2] > 0.005, self.ecm[:,:,2], 0)
    
    def reduce_chemo_concentration(self, halflife):
        self.ecm[:,:,2] = self.ecm[:,:,2] * (.5 ** (1/halflife))
        self.ecm[:,:,2] = np.where(self.ecm[:,:,2] > 0.005, self.ecm[:,:,2], 0)
    
        
class Tissue():
    def __init__(self, size):
        self.size = size
        self.cell_matrix = np.full((size, size), None, dtype= object)
        for row in range(size):
            for column in range(size):
                rand_chance = np.random.uniform()
                if rand_chance < .99:
                    self.cell_matrix[row, column] = Cell(row, column, size, is_cancer = False)
                else:
                    self.cell_matrix[row, column] = None
        
        self.extracellular_matrix = Extracellular_Matrix(size)
        
        a = np.indices((size,size))
        self.all_indices = list(zip(a[0].ravel(), a[1].ravel()))
    
    def simulate_step(self, chemo_concentration = None):
        if chemo_concentration != None:
            self.extracellular_matrix.set_chemo_concentration(chemo_concentration)
        else:
            self.extracellular_matrix.reduce_chemo_concentration(18)
            
        next_step_cells = np.full((self.size, self.size), None, dtype= object)
        sample_indices = np.random.choice(len(self.all_indices), len(self.all_indices), replace=False)
        for cell_index in sample_indices:
            cell_location = self.all_indices[cell_index]
            cell_row, cell_column = cell_location
        
            cell = copy.deepcopy(self.cell_matrix[cell_row, cell_column])
            
            if cell is None:
                continue
            
            cell_neighbors = cell.neighbors
            print
            is_empty_neighbor = False
            neighbor_status = []
            for neighbor in cell_neighbors:
                if self.cell_matrix[neighbor[0], neighbor[1]] == None:
                    is_empty_neighbor = True
                    neighbor_status.append(True)
                else:
                    neighbor_status.append(False)
                    
            environment_nutrient = self.extracellular_matrix.ecm[cell_row, cell_column, 0]
            environment_oxygen = self.extracellular_matrix.ecm[cell_row, cell_column, 1]
            environment_chemo = self.extracellular_matrix.ecm[cell_row, cell_column, 2]
            
            if cell.is_dying and cell.stage_counter == 12:
                next_step_cells[cell_row, cell_column] = None
            else:
                is_dividing = cell.update_cell(environment_nutrient, environment_oxygen, environment_chemo, is_empty_neighbor)
                
                if is_dividing and is_empty_neighbor:
                    replacement_index = np.random.choice(np.where(np.array(neighbor_status))[0])
                    replacement_location = cell_neighbors[replacement_index]
                    cell_atributes = cell.get_cell_attributes()
                    next_step_cells[replacement_location[0], replacement_location[1]] = Cell(replacement_location[0], replacement_location[1], self.size, cell_atributes[6], cell_atributes[:6])
                
                elif is_dividing and is_empty_neighbor == False:
                    replacement_index = np.random.choice(range(len(cell_neighbors)))
                    replacement_location = cell_neighbors[replacement_index]
                    cell_atributes = cell.get_cell_attributes()
                    next_step_cells[replacement_location[0], replacement_location[1]] = Cell(replacement_location[0], replacement_location[1], self.size, cell_atributes[6], cell_atributes[:6])
                
                next_step_cells[cell_row, cell_column] = cell
        
        self.cell_matrix = next_step_cells
    
    def generate_plot_array(self):
        plotting_array = np.full((self.size, self.size), 0, dtype=int)
        for row in range(self.size):
            for column in range(self.size):
                if self.cell_matrix[row,column] != None:
                    # Then check if the cell is cancer
                    if self.cell_matrix[row,column].is_cancer:
                        plotting_array[row, column] = 2
                    # Otherwise, mark normal cells
                    else:
                        plotting_array[row, column] = 1
        return plotting_array
    
    def place_cancer_cell(self, row, column):
        self.cell_matrix[row, column] = Cell(row, column, self.size, True)
    
    def get_metric(self, metric = "Number of Cancer"):
        metric_array = np.full((self.size, self.size), 0, dtype=float)
        if metric == "Number of Cancer":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        # Then check if the cell is cancer
                        if self.cell_matrix[row,column].is_cancer:
                            metric_array[row, column] = 1
                        # Otherwise, mark normal cells
                        else:
                            metric_array[row, column] = 0
                            
        if metric == "Number of Mutations":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].num_mutation
        
        if metric == "Apoptosis Resistance":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].resistance_apoptosis
        
        if metric == "Recovery Rate":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].recovery_rate
        
        if metric == "Mutation Probability":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].mutation_probability
        
        if metric == "Invasion Probability":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].invasion_prob
        if metric == "Damage":
            for row in range(self.size):
                for column in range(self.size):
                    if self.cell_matrix[row,column] != None:
                        metric_array[row, column] = self.cell_matrix[row, column].damage
        return metric_array


file_name = "Test"
chemo_admin_interval_days = 28
chemo_admin_interval = chemo_admin_interval_days * 24
chemo_max_concentration = 0.6 
treatment_delay_days = 60
treatment_delay = treatment_delay_days * 24


test_duration_days = 250
num_steps = test_duration_days * 24
size = 41
tissue = Tissue(size)
tissue.place_cancer_cell(20, 30)



#%%
cmap = ListedColormap(["Black", mcolors.CSS4_COLORS['salmon'],mcolors.CSS4_COLORS['saddlebrown']])

plotting_array = tissue.generate_plot_array()
frames = []
metrics = {"Number of Cancer": [],
           "Number of Mutations": [], 
           "Apoptosis Resistance": [],
           "Recovery Rate": [],
           "Mutation Probability": [],
           "Invasion Probability": [], 
           "Damage": []}

def get_all_metrics():
    for metric in metrics.keys():
        metrics[metric].append(tissue.get_metric(metric))
        

# plt.figure()
# img = plt.imshow(plotting_array, cmap=cmap)
# plt.pause(.1)

    
# plt.pause(5)
start_time = time.time()
for step in range(num_steps):
    if step >= treatment_delay and (step - treatment_delay) % chemo_admin_interval == 0:
        print("Delivered")
        tissue.simulate_step(chemo_max_concentration)
    else:
        tissue.simulate_step()
    
    if step % 12 == 0:
        print(step)
        plotting_array = tissue.generate_plot_array()
        frames.append(plotting_array)
        # if step % 96 == 0:
        #     img.set_data(plotting_array)
        #     plt.draw()
        #     plt.title(f" Day {step/24:.2f}")
        #     plt.pause(.1)
        get_all_metrics()
        if np.sum(tissue.get_metric()) == 0:
            print("Cancer Eradicated")
            break

        
end_time = time.time()
print(f"Time to simulate {test_duration_days} days: {end_time - start_time}")

with open(f'{file_name}.pkl', 'wb') as f:
    pickle.dump(metrics, f)
#%%
fig, ax = plt.subplots() 
img_anim = ax.imshow(frames[0], cmap=cmap)
title = ax.set_title("Day 0.00") 

def update(frame_idx):
    img_anim.set_data(frames[frame_idx])
    
    days = frame_idx * 12 / 24  
    title.set_text(f"Day {days:.2f}") 
    
    return [img_anim, title]

ani = FuncAnimation(fig, update, frames=len(frames), interval=500)

ani.save(f"{file_name}.mp4", fps=15, dpi=150, writer='ffmpeg')




#test_ecm = Extracellular_Matrix(41)
# nutrient = test_ecm.ecm[:,:,0]
# oxygen = test_ecm.ecm[:,:,1]
# chemo = test_ecm.ecm[:,:,2]
# print(chemo)
# test_ecm.set_chemo_concentration(.5)
# chemo = test_ecm.ecm[:,:,2]
# print(chemo)
# test_ecm.reduce_chemo_concentration(24)
# chemo_2 = test_ecm.ecm[:,:,2]

# test_cell = Cell(1,2, 4, False)
# test_cell_2 = Cell(1,2, 4, False, [10,0.1, 0.1, 1.5, 0.3, 0.03])
# test_cell_3 = Cell(1,2,4, True)

# test_tissue = Tissue(11)

        
    
    
