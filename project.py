import tkinter as tk
import random
import pandas as pd
from time import strftime
import numpy as np
import matplotlib.pyplot as plt
window = tk.Tk() # you may also see it named as "root" in other sources

window.title("FLIGHT SCRAPER") # self explanatory!
window.minsize(width=600,height=400)
window.maxsize(width=600,height=600)

window.resizable(width="false", height="false") # change to false if you want to prevent resizing
frame_header = tk.Frame(window, borderwidth=2, pady=2)
center_frame = tk.Frame(window, borderwidth=2, pady=5)
bottom_frame = tk.Frame(window, borderwidth=2, pady=5)
frame_header.grid(row=0, column=0)
center_frame.grid(row=1, column=0)
bottom_frame.grid(row=2, column=0)

# label header to be placed in the frame_header
header = tk.Label(frame_header, text = "SPECIES EVOLUTION", bg='green', fg='white', height='3', width='50', font=("Helvetica 16 bold"))
# inside the grid of frame_header, place it in the position 0,0
header.grid(row=0, column=0)
# WIDGETS
number1Label = tk.Label (center_frame,text="Enter Number of Generations")
number1Label.pack()
number1Box = tk.Entry(center_frame)
number1Box.pack()

number2Label = tk.Label (center_frame,text="Enter the Gene Length")
number2Label.pack()
number2Box = tk.Entry(center_frame)
number2Box.pack()


def LCS_function(X, Y): # X and Y are the strings to compare
      
    m = len(X) 
    n = len(Y) 
  
    L = [[0 for i in range(n+1)] for j in range(2)]  
    bi = bool
      
    for i in range(m): 
        bi = i&1
  
        for j in range(n+1): 
            if (i == 0 or j == 0): 
                L[bi][j] = 0
  
            elif (X[i] == Y[j - 1]): 
                L[bi][j] = L[1 - bi][j - 1] + 1
  
            else: 
                L[bi][j] = max(L[1 - bi][j],  
                               L[bi][j - 1]) 
  
    return L[bi][n] #returns the length of longest common sequence






class Genetic_Class:
    def __init__(self, refrence_gene_number, reference_genes, population_number, gene_length, max_generation):  #reference_genes is a nd array of reference genes. If you want to have random reference genes, just pass an integer.
        
        if refrence_gene_number <2 or refrence_gene_number>10:
            raise Exception('refrence_gene_number is less than 2 or more than 10')
        if population_number <10 or population_number>100:
            raise Exception('population_number is less than 10 or more than 100')
        
        self.reference_gene_number = 5
        self.population_number = 100
        self.gene_length = gene_length
        self.max_generation = max_generation
        self.generation = 0
        #print(self.population_number)
        
        char = ['A', 'C', 'G', 'T']
        self.population = np.empty([2* population_number], dtype='<U'+str(gene_length))
        for i in range(population_number):
            temp = ''
            for j in range(gene_length):
                temp += np.random.choice(char)
            self.population[i] = temp
        #print(self.population)
        
        if type(reference_genes) == int:
            self.reference_genes = np.empty([refrence_gene_number], dtype='<U'+str(gene_length))
            for i in range(refrence_gene_number):
                temp = ''
                for j in range(gene_length):
                    temp += np.random.choice(char)
                self.reference_genes[i] = temp
        else:
            self.reference_genes = reference_genes
        del char
        #print(self.population_number)
    
    
    
    def cross_over(self, A, B):                 # 30% of each gene goes to cross over
        temp_start = np.random.randint(len(A))
        temp_end = temp_start + int(0.3 * self.gene_length)
        if temp_end > self.gene_length:
            temp = temp_end - self.gene_length
            temp_end = temp_start
            temp_start = temp
        child_one = A[0:temp_start] + B[temp_start:temp_end] + A[temp_end:]
        child_two = B[0:temp_start] + A[temp_start:temp_end] + B[temp_end:]
        return child_one, child_two
    
    
    
    def mutation(self, A):                      # 50% of each gene gets mutated
        result = np.empty([len(A)], dtype=str)
        for j in range(len(A)):
            result[j] = A[j]
        for i in range(int(self.gene_length * 0.5)):
            result[np.random.randint(self.gene_length-1)] = np.random.choice(['A', 'C', 'G', 'T'])
        return ''.join(result)
    

    
    def squared_average(self,A):
        temp = np.sqrt(np.sum(np.power(A,2)))/A.shape[0]
        return temp
        


    def next_generation(self,generation_number):
        if generation_number == 'max':
            generation_number = self.max_generation - self.generation
        if self.max_generation - self.generation < generation_number:
            raise Exception('generation_number parameter of next_generation function is more than allowed.')
        scores = np.empty([2*self.population_number, self.reference_gene_number], dtype=int)
        squared_average_score = np.empty([2*self.population_number], dtype=float)
        for rep in range(generation_number):
            self.generation += 1
            perm = np.random.permutation(np.array(range(self.population_number)))
            for i in range(0, int(0.7*perm.shape[0]), 2):      # 0.7 of population go to cross over
                ch1, ch2 = self.cross_over(self.population[perm[i]], self.population[perm[i+1]])
                self.population[i+self.population_number] = ch1
                self.population[i+self.population_number+1] = ch2
                
            for i in range(int(0.7*perm.shape[0]), perm.shape[0]):  # 0.3 of population go to mutation
                self.population[i+self.population_number] = self.mutation(self.population[perm[i]])
                
            if self.generation == 0:
                for i in range(self.population_number):
                    for j in range(self.reference_gene_number):
                        scores[i,j] = LCS_function(str(self.population[i]), str(self.reference_genes[j]))
                    squared_average_score[i] = self.squared_average(scores[i])

            for i in range(self.population_number, 2*self.population_number):
                for j in range(self.reference_gene_number):
                    scores[i,j] = LCS_function(str(self.population[i]), str(self.reference_genes[j]))
                squared_average_score[i] = self.squared_average(scores[i])
                
            scored_population = np.empty([2*self.population_number, 3], dtype=object)
            scored_population.T[0] = squared_average_score
            scored_population.T[1] = self.population
            scored_population.T[2] = np.array(range(2*self.population_number))
            scored_population = scored_population[np.argsort(scored_population[:,0])]
            scored_population = np.flip(scored_population,0)
            self.population = scored_population[:,1]
            order = scored_population[:,2]                    #indexes of sorted population
            temp_scores = np.empty([2*self.population_number, self.reference_gene_number], dtype=int)
            for i in range(2*self.population_number):
                temp_scores[i] = scores[order[i]]
            scores = temp_scores
            squared_average_score = scored_population[:,0]
            del scored_population
            del temp_scores
        population_score = scores[0:self.population_number]
        
        return np.array((population_score/self.gene_length)*100, dtype=int)         # reurns the percentage of similarity of each gene to each reference gene.












def addNo():
    a = int(number1Box.get())
    b = int(number2Box.get())
    GC = Genetic_Class(refrence_gene_number=5, reference_genes=-1, population_number=100, gene_length=a, max_generation=50)   #Initialization of genetic algorithm
    #results1 = GC.next_generation(6)     # run 7 generations and return the result after 7 generations
    # result is a matrix of size (population_number*refrence_gene_number), that contains the percentage of similarity of each gene of population to each reference gene.
    #results2 = GC.next_generation(3)    # run 4 more generations and return the result after 4 more generations (returns the result of 7 + 4 = 11th generation)
    
    results3 = GC.next_generation(a)    # run remaining generations and return the result. (generation 12 to max generation)
    
    count = [0,0,0,0,0]
    
    for i in range(results3.shape[0]):
        
        res = results3[i]
        largest = max(res)
        indices = [idx for idx, val in enumerate(res) if val == largest]
        bucket = random.choice(indices)
            
            
        count[bucket] = count[bucket] + 1
            #print(bucket)
    #print(results3[1][2])
    #print(count)
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    references = ["Swimming","Gliding","Flying","Jumping","Running"]
    
    ax.bar(references,count)
    plt.show()    

    
but = tk.Button(bottom_frame,text="Generate Data", command=addNo)
but.pack()

window.mainloop()

