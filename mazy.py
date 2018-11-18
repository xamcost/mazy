#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 19:40:06 2018

@author: costalongam
"""

import random
import numpy as np
import matplotlib.pyplot as plt

class Cell(object):
    """ 
    A cell object as an element of a maze
    
    Parameters
    ----------
    
    i: int
        row position in maze
    j: int
        column position in maze
    """
    def __init__(self, i, j):
        self.i, self.j = i, j
        self.walls = {'top':True, 'bottom':True, 'left':True, 'right':True}
        self.toDeadEnd = False
    
    def break_wall(self, other):
        """ 
        Break the wall between cell self and cell other
        
        Parameters
        ----------
        
        other: Cell object
            row position in maze
        """
        if self.i == other.i:
            if self.j < other.j:
                self.walls['right'] = False
                other.walls['left'] = False
            elif self.j > other.j:
                self.walls['left'] = False
                other.walls['right'] = False
            else:
                raise ValueError('Can break a wall only between two neighboring cells')
        elif self.j == other.j:
            if self.i < other.i:
                self.walls['bottom'] = False
                other.walls['top'] = False
            elif self.i > other.i:
                self.walls['top'] = False
                other.walls['bottom'] = False
            else:
                raise ValueError('Can break a wall only between two neighboring cells')
        else:
            raise ValueError('Can break a wall only between two neighboring cells')
            
    def __repr__(self):
        return 'cell at ({0}, {1})'.format(self.i, self.j)


class Maze(object):
    """ 
    Creates a p-rows and q-columns maze, starting from cell located at i0, j0
    
    Parameters
    ----------
    
    p: int
        number of rows of the maze
    q: int
        number of columns of the maze
    i0: int
        row-coordinate of the starting point for maze construction
        Default: 0
    j0: int
        column-coordinate of the starting point for maze construction
        Default: 0
    """
    def __init__(self, p, q, i0=0, j0=0):
        self.p, self.q = p, q
        self.i0, self.j0 = i0, j0
        self.maze_map = [[Cell(i, j) for j in range(q)] for i in range(p)]
        self.make()
        
    def make(self):
        """ Build the maze using deep first search algorithm """
        cell_stack = [self.maze_map[self.i0][self.j0]]
        nv = 1
        N = self.p * self.q
        while nv < N:
            neighbours = self.get_unvisited_neighbours(cell_stack[-1])
            if not neighbours:
                cell_stack.pop()
                continue
            cell_stack.append(random.choice(neighbours))
            Cell.break_wall(cell_stack[-2], cell_stack[-1])
            nv += 1
            
    def get_unvisited_neighbours(self, cell):
        """ 
        Returns the unvisited neighbours of a cell
        
        Parameters
        ----------
        
        cell: Cell object
            cell of interest
            
        Returns
        ----------
        
        list of cell object corresponding to unvisited neighbours
        """
        delta = [(-1,0), (1,0), (0,1), (0,-1)]
        neighbours = []
        for di, dj in delta:
            i2, j2 = cell.i + di, cell.j + dj
            if (0 <= i2 < self.p) and (0 <= j2 < self.q):
                neighbour = self.maze_map[i2][j2]
                if all(neighbour.walls.values()):
                    neighbours.append(neighbour)
        return neighbours
    
    def get_image(self):
        """ Returns an image of the maze """
        im = np.ones((10*self.p + 1, 10*self.q + 1))
        for i in range(self.p):
            for j in range(self.q):
                if self.maze_map[i][j].walls['top']:
                    im[10*i, 10*j:(10*(j + 1) + 1)] = 0
                if self.maze_map[i][j].walls['left']:
                    im[10*i:(10*(i + 1) + 1), 10*j] = 0
        im[10*self.p, :] = 0
        im[:, 10*self.q] = 0
        return im
    
    def plot(self):
        """ Dsplay the maze in a matplotlib figure """
        fig, ax = plt.subplots()
        ax.imshow(self.get_image(), cmap='gray')
        ax.tick_params(axis='both', bottom=False, top=False, labelbottom =False, 
                       left=False, right=False, labeltop =False, 
                       labelright =False, labelleft =False)
        fig.show()
        
    def __str__(self):
        maze_rows = ['-' * self.q*2]
        for i in range(self.p):
            maze_row = ['|']
            for j in range(self.q):
                if self.maze_map[i][j].walls['right']:
                    maze_row.append(' |')
                else:
                    maze_row.append('  ')
            maze_rows.append(''.join(maze_row))
            maze_row = ['|']
            for j in range(self.q):
                if self.maze_map[i][j].walls['bottom']:
                    maze_row.append('-+')
                else:
                    maze_row.append(' +')
            maze_rows.append(''.join(maze_row))
        return '\n'.join(maze_rows)
    
    def __repr__(self):
        return '{0} x {1} maze'.format(self.p, self.q)


if __name__ == "__main__":
    m = Maze(4, 5)
    print(m)
    m.plot()