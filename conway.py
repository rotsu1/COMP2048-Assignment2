# -*- coding: utf-8 -*-
"""
The Game of Life (GoL) module named in honour of John Conway

This module defines the classes required for the GoL simulation.

Created on Tue Jan 15 12:21:17 2019

@author: shakes
"""
import numpy as np
from scipy import signal
import rle

class GameOfLife:
    '''
    Object for computing Conway's Game of Life (GoL) cellular machine/automata
    '''
    def __init__(self, N=256, finite=False, fastMode=False):
        self.grid = np.zeros((N, N), np.int64)
        self.neighborhood = np.ones((3,3), np.int64) # 8 connected kernel
        self.neighborhood[1,1] = 0 #do not count centre pixel
        self.finite = finite
        self.fastMode = fastMode
        self.aliveValue = 1
        self.deadValue = 0
        self.length = N
        self.row = N
        self.column = N
        self.pad = 0
        
    def getStates(self):
        '''
        Returns the current states of the cells
        '''
        return self.grid
    
    def getGrid(self):
        '''
        Same as getStates()
        '''
        return self.getStates()
               
    def evolve(self):
        '''
        Given the current states of the cells, apply the GoL rules:
        - Any live cell with fewer than two live neighbors dies, as if by underpopulation.
        - Any live cell with two or three live neighbors lives on to the next generation.
        - Any live cell with more than three live neighbors dies, as if by overpopulation.
        - Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction
        '''
        #get weighted sum of neighbors
        #PART A & E CODE HERE

        #implement the GoL rules by thresholding the weights
        #PART A CODE HERE

        def live(number: int, isLive: int):
            if isLive:
                # next generation
                if number == 2 or number == 3:
                    return self.aliveValue
                # underpopulation or overpopulation
                else: 
                    return self.deadValue
            elif not isLive:
                # reproduction
                if number == 3:
                    return self.aliveValue
                else:
                    return self.deadValue
        
        new_grid = np.zeros((self.row, self.column), np.int64)
        if not self.fastMode:
            for row_index, row in enumerate(self.grid):
                for column_index, column in enumerate(row):
                    if self.pad > 0:
                        if row_index < self.pad or column_index < self.pad or (self.row - row_index) <= self.pad or (self.column - column_index) <= self.pad:
                            pass
                        else:
                            total = self.grid[row_index - 1, column_index - 1] + self.grid[row_index - 1, column_index] + self.grid[row_index - 1, column_index + 1] + self.grid[row_index, column_index - 1] + self.grid[row_index, column_index + 1] + self.grid[row_index + 1, column_index - 1] + self.grid[row_index + 1, column_index] + self.grid[row_index + 1, column_index + 1]
                            new_grid[row_index, column_index] = live(total, column)
                    else:
                        if row_index == 0 and column_index == 0:
                            total = self.grid[0, 1] + self.grid[1, 0] + self.grid[1, 1]
                            new_grid[0, 0] = live(total, column)
                        elif row_index == 0 and column_index == self.length - 1:
                            total = self.grid[0, self.length - 2] + self.grid[1, self.length - 2] + self.grid[1, self.length - 1]
                            new_grid[0, self.length - 1] = live(total, column)
                        elif row_index == self.length - 1 and column_index == 0:
                            total = self.grid[self.length - 2, 0] + self.grid[self.length - 2, 1] + self.grid[self.length - 1, 1]
                            new_grid[self.length - 1, 0] = live(total, column)
                        elif row_index == self.length - 1 and column_index == self.length - 1:
                            total = self.grid[self.length - 2, self.length - 1] + self.grid[self.length - 2, self.length - 2], self.grid[self.length - 1, self.length - 2]
                            new_grid[self.length - 1, self.length - 1] = live(total, column)
                        elif row_index == 0:
                            total = self.grid[0, column_index - 1] + self.grid[1, column_index - 1] + self.grid[1, column_index] + self.grid[1, column_index + 1] + self.grid[0, column_index + 1]
                            new_grid[0, column_index] = live(total, column)
                        elif row_index == self.length - 1:
                            total = self.grid[self.length - 1, column_index - 1] + self.grid[self.length - 2, column_index - 1] + self.grid[self.length - 2, column_index] + self.grid[self.length - 2, column_index + 1] + self.grid[self.length - 1, column_index + 1]
                            new_grid[self.length - 1, column_index] = live(total, column)
                        elif column_index == 0:
                            total = self.grid[row_index - 1, 0] + self.grid[row_index - 1, 1] + self.grid[row_index, 1] + self.grid[row_index + 1, 1] + self.grid[row_index + 1, 0]
                            new_grid[row_index, 0] = live(total, column)
                        elif column_index == self.length - 1:
                            total = self.grid[row_index - 1, self.length - 1] + self.grid[row_index - 1, self.length - 2] + self.grid[row_index, self.length - 2] + self.grid[row_index + 1, self.length - 2] + self.grid[row_index + 1, self.length - 1]
                            new_grid[row_index, self.length - 1] = live(total, column)
                        else:
                            total = self.grid[row_index - 1, column_index - 1] + self.grid[row_index - 1, column_index] + self.grid[row_index - 1, column_index + 1] + self.grid[row_index, column_index - 1] + self.grid[row_index, column_index + 1] + self.grid[row_index + 1, column_index - 1] + self.grid[row_index + 1, column_index] + self.grid[row_index + 1, column_index + 1]
                            new_grid[row_index, column_index] = live(total, column)
            #update the grid
            self.grid = new_grid
        else:
            # set mode=same to return the same size as self.grid
            new_grid = signal.convolve(self.grid, self.neighborhood, mode='same')
            new_grid[((new_grid != 2) & (new_grid != 3)) & (self.grid == 1)] = self.deadValue
            new_grid[(new_grid != 3) & (self.grid == 0)] = self.deadValue
            new_grid[((new_grid == 2) | (new_grid == 3)) & (self.grid == 1)] = self.aliveValue
            new_grid[(new_grid == 3) & (self.grid == 0)] = self.aliveValue
            self.grid = new_grid
    
    def insertBlinker(self, index=(0,0)):
        '''
        Insert a blinker oscillator construct at the index position
        '''
        self.grid[index[0], index[1]+1] = self.aliveValue
        self.grid[index[0]+1, index[1]+1] = self.aliveValue
        self.grid[index[0]+2, index[1]+1] = self.aliveValue
        
    def insertGlider(self, index=(0,0)):
        '''
        Insert a glider construct at the index position
        '''
        self.grid[index[0], index[1]+1] = self.aliveValue
        self.grid[index[0]+1, index[1]+2] = self.aliveValue
        self.grid[index[0]+2, index[1]] = self.aliveValue
        self.grid[index[0]+2, index[1]+1] = self.aliveValue
        self.grid[index[0]+2, index[1]+2] = self.aliveValue
        
    def insertGliderGun(self, index=(0,0)):
        '''
        Insert a glider construct at the index position
        '''
        self.grid[index[0]+1, index[1]+25] = self.aliveValue
        
        self.grid[index[0]+2, index[1]+23] = self.aliveValue
        self.grid[index[0]+2, index[1]+25] = self.aliveValue
        
        self.grid[index[0]+3, index[1]+13] = self.aliveValue
        self.grid[index[0]+3, index[1]+14] = self.aliveValue
        self.grid[index[0]+3, index[1]+21] = self.aliveValue
        self.grid[index[0]+3, index[1]+22] = self.aliveValue
        self.grid[index[0]+3, index[1]+35] = self.aliveValue
        self.grid[index[0]+3, index[1]+36] = self.aliveValue
        
        self.grid[index[0]+4, index[1]+12] = self.aliveValue
        self.grid[index[0]+4, index[1]+16] = self.aliveValue
        self.grid[index[0]+4, index[1]+21] = self.aliveValue
        self.grid[index[0]+4, index[1]+22] = self.aliveValue
        self.grid[index[0]+4, index[1]+35] = self.aliveValue
        self.grid[index[0]+4, index[1]+36] = self.aliveValue
        
        self.grid[index[0]+5, index[1]+1] = self.aliveValue
        self.grid[index[0]+5, index[1]+2] = self.aliveValue
        self.grid[index[0]+5, index[1]+11] = self.aliveValue
        self.grid[index[0]+5, index[1]+17] = self.aliveValue
        self.grid[index[0]+5, index[1]+21] = self.aliveValue
        self.grid[index[0]+5, index[1]+22] = self.aliveValue
        
        self.grid[index[0]+6, index[1]+1] = self.aliveValue
        self.grid[index[0]+6, index[1]+2] = self.aliveValue
        self.grid[index[0]+6, index[1]+11] = self.aliveValue
        self.grid[index[0]+6, index[1]+15] = self.aliveValue
        self.grid[index[0]+6, index[1]+17] = self.aliveValue
        self.grid[index[0]+6, index[1]+17] = self.aliveValue
        # missing
        self.grid[index[0]+6, index[1]+18] = self.aliveValue
        self.grid[index[0]+6, index[1]+23] = self.aliveValue
        self.grid[index[0]+6, index[1]+25] = self.aliveValue
        
        self.grid[index[0]+7, index[1]+11] = self.aliveValue
        self.grid[index[0]+7, index[1]+17] = self.aliveValue
        self.grid[index[0]+7, index[1]+25] = self.aliveValue
        
        self.grid[index[0]+8, index[1]+12] = self.aliveValue
        self.grid[index[0]+8, index[1]+16] = self.aliveValue
        
        self.grid[index[0]+9, index[1]+13] = self.aliveValue
        self.grid[index[0]+9, index[1]+14] = self.aliveValue
        
    def insertFromPlainText(self, txtString, pad=0):
        '''
        Assumes txtString contains the entire pattern as a human readable pattern without comments
        '''
        self.pad = pad

        with open(txtString, 'r') as file:
            ignored_lines = 0
            whole_lines = 0
            length_list = []
            if pad > 0:
                for line in file:
                    whole_lines += 1
                    if line.startswith('!'):
                        ignored_lines += 1
                    elif not line.startswith('!'):
                        length_list.append(len(line))
                self.row = whole_lines - ignored_lines + pad * 2
                self.column = max(length_list) + pad * 2
                # pad 4 sides
                self.grid = np.zeros((self.row, self.column), np.int64)

        with open(txtString, 'r') as file:
            for row, line in enumerate(file):
                row -= ignored_lines
                for column, character in enumerate(line):
                    if character == 'O':
                        self.grid[row+pad, column+pad] = self.aliveValue

    def insertFromRLE(self, rleString, pad=0):
        '''
        Given string loaded from RLE file, populate the game grid
        '''

        self.pad = pad
        
        def parser(line):
            '''
            Changes into .cells format
            '''
            number = ""
            string = ""
            list = []
            for char in line:
                if char.isdigit():
                    number += char
                elif char == "b":
                    if not number:
                        string += '.'
                    else:
                        i = 0
                        while i < int(number):
                            i += 1
                            string += '.'
                    number = ""
                elif char == 'o':
                    if not number:
                        string += 'O'
                    else:
                        i = 0
                        while i < int(number):
                            i += 1
                            string += 'O'
                    number = ""
            return string
        
        with open(rleString, 'r') as file:
            ignored_lines = 0
            for line in file:
                if line.startswith('#'):
                    ignored_lines += 1
                elif line.startswith('x'):
                    ignored_lines += 1
                    elements = line.split(',')
                    self.row = int(elements[1].split('=')[1].strip()) + pad * 2
                    self.column = int(elements[0].split('=')[1].strip()) + pad * 2
            self.grid = np.zeros((self.row, self.column), np.int64)

        with open (rleString, 'r') as file:    
            for _ in range(ignored_lines):
                next(file)
            joined_lines = " ".join(line.strip() for line in file)
            separated_lines = joined_lines.split("$")
            terminate = False
            for row, line in enumerate(separated_lines):
                parsed_line = parser(line)
                for column, character in enumerate(parsed_line):
                    if character == 'O':
                        self.grid[row+pad, column+pad] = self.aliveValue
                    elif character == '!':
                        terminate = True
                        break
                if terminate == True:
                    break

# Turing complete is when you can simultate any other turing machines.
# Turing machine consists read/write head, memory(tape), states, and set of rules(transition function).
# The link shows that there are patterns that represent reading and writing, states, set of rules (GoL rules), and (assuming the grid is infinite) the grid serves as the tape.
# Therefore, game of life can run a turing machine so it is turing complete.