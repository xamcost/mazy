#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 19:40:06 2018

@author: costalongam
"""

import random
import numpy as np
from fractions import Fraction
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
        self.prim_visited = False
        self.solve_visited = False
        self.type = None
    
    def break_wall(self, other):
        """ 
        Break the wall between cell self and cell other
        
        Parameters
        ----------
        
        other: Cell object
            neighbouring cell
        """
        if self.i == other.i:
            if self.j == other.j - 1:
                self.walls['right'] = False
                other.walls['left'] = False
            elif self.j == other.j + 1:
                self.walls['left'] = False
                other.walls['right'] = False
            else:
                raise ValueError('Can break a wall only between two neighboring cells')
        elif self.j == other.j:
            if self.i == other.i - 1:
                self.walls['bottom'] = False
                other.walls['top'] = False
            elif self.i == other.i + 1:
                self.walls['top'] = False
                other.walls['bottom'] = False
            else:
                raise ValueError('Can break a wall only between two neighboring cells')
        else:
            raise ValueError('Can break a wall only between two neighboring cells')
            
    def isDeadEnd(self):
        """ Return True if the cell has 3 walls """
        comp = 0
        for w in self.walls.values():
            if w:
                comp += 1
        return comp == 3
    
    def get_common_wall(self, other):
        """ 
        Return the common wall between cell self and cell other, with respect 
        to self
        
        Parameters
        ----------
        
        other: Cell object
            neighbouring cell
            
        Returns
        ---------- 
        
        str of corresponding wall
        """
        if self.i == other.i:
            if self.j < other.j:
                return 'right'
            elif self.j > other.j:
                return 'left'
        elif self.j == other.j:
            if self.i < other.i:
                return 'bottom'
            elif self.i > other.i:
                return 'top'
        else:
            raise ValueError('{0} and {1} are not neighbouring cells'.format(self, other))
            
    def __repr__(self):
        return 'cell at ({0}, {1})'.format(self.i, self.j)
    
    def __eq__(self, other):
        return self.i == other.i and self.j == other.j
    
    def __hash__(self):
        return hash((self.i, self.j))


class Maze(object):
    """ 
    Creates a p-rows and q-columns maze, starting from cell located at i0, j0
    
    Parameters
    ----------
    
    p: int
        number of rows of the maze
    q: int
        number of columns of the maze
    method: str, optional
        algorithm for maze generation. Possible values: dfs (depth first search), 
        prim
        Default: dfs
    braid_factor: float, optional
        Between 0 and 1. If 0, perfect maze (only one solution). If 1, maze 
        has no dead ends.
        Default: 0
    i0: int, optional
        row-coordinate of the starting point for maze construction. 
        Default: 0
    j0: int
        column-coordinate of the starting point for maze construction. 
        Default: 0
    """
    def __init__(self, p, q, method='dfs', braid_factor=0, i0=0, j0=0):
        self.p, self.q = p, q
        self.i0, self.j0 = i0, j0
        self.method = method
        self.braid_factor = braid_factor
        self.maze_map = [[Cell(i, j) for j in range(q)] for i in range(p)]
        self.make()
        
    def make(self):
        """ Build the maze using algorithm method. If braid_factor is 0, mazes are 
        perfect and have only one solution. """
        if self.method == 'dfs':
            cell_stack = [self.maze_map[self.i0][self.j0]]
            nv = 1
            N = self.p * self.q
            while nv < N:
                neighbours = self.get_neighbours(cell_stack[-1], kind='unvisited')
                if not neighbours:
                    cell_stack.pop()
                    continue
                cell_stack.append(random.choice(neighbours))
                Cell.break_wall(cell_stack[-2], cell_stack[-1])
                nv += 1
        elif self.method == 'prim':
            current_cell = self.maze_map[self.i0][self.j0]
            current_cell.prim_visited = True
            cell_stack = self.get_neighbours(current_cell)
            next_cell = random.choice(cell_stack)
            Cell.break_wall(current_cell, next_cell)
            next_cell.prim_visited = True
            cell_stack = list(set(cell_stack).union(self.get_neighbours(next_cell, kind='unvisited')))
            cell_stack.remove(next_cell)
            while cell_stack:
                next_cell = random.choice(cell_stack)
                next_cell.prim_visited = True
                valid_neighbours = [c for c in self.get_neighbours(next_cell) if c.prim_visited]
                if valid_neighbours:
                    other_cell = random.choice(valid_neighbours)
                    Cell.break_wall(next_cell, other_cell)
                    cell_stack = list(set(cell_stack).union(self.get_neighbours(next_cell, kind='unvisited')))
                cell_stack.remove(next_cell)
        else:
            raise ValueError('{0} is an unknow/unsupported method for maze generation'.format(self.method))
        self.break_dead_ends()  
            
    def solve(self, start=None, end=None, method='dead_end_filler'):
        """ 
        Solves the maze, starting at start and ending at end, using method
        
        Parameters
        ----------
        
        start: tuple of int
            starting cell coordinates. Randomly chosen if not provided.
        end: tuple of int
            ending cell coordinates. Randomly chosen if not provided.
        method: str
            method used to solve the maze
        """
        if start is None and not hasattr(self, 'start'):
            start = (random.randint(0, self.p - 1), random.randint(0, self.q - 1))
        if end is None and not hasattr(self, 'end'):
            end = (random.randint(0, self.p - 1), random.randint(0, self.q - 1))
            while end == start:
                end = (random.randint(0, self.p - 1), random.randint(0, self.q - 1))
        self.start = start
        self.end = end
        self.maze_map[self.start[0]][self.start[1]].type = 'start'
        self.maze_map[self.end[0]][self.end[1]].type = 'end'
        
        if method == 'dead_end_filler':
            new_dead_ends = self.fill_dead_end(tag='direct_dead_end')
            restore = new_dead_ends
            while new_dead_ends:
                for c, w in new_dead_ends:
                    c.walls[w] = True
                new_dead_ends = self.fill_dead_end(tag='indirect_dead_end')
                restore = restore + new_dead_ends
            for c, w in restore:
                c.walls[w] = False
        else:
            raise ValueError('{0} is not a known/implemented method for maze solving'.format(method))
            
    def get_dead_ends(self):
        """ Returns the list of all dead ends in the maze, except starting and 
        ending cells """
        dead_ends = []
        for i in range(self.p):
            for j in range(self.q):
                c = self.maze_map[i][j]
                if c.isDeadEnd():
                    dead_ends.append(c)
        return dead_ends
            
    def fill_dead_end(self, tag='direct_dead_end'):
        """ 
        Fills 'direct' dead ends in a maze 
        
        Parameters
        ----------
        
        tag: str
            tag to assign to cell found in dead ends.
            
        Returns
        ----------
        
        list of tuples of type (cell, wall) indicating walls to create to 
        isolate dead ends.
        """
        dead_ends = []
        for i in range(self.p):
            for j in range(self.q):
                c = self.maze_map[i][j]
                if c.isDeadEnd() and not c.solve_visited:
                    if (i, j) != self.start:
                        if (i, j) != self.end:
                            dead_ends.append(c)
                            c.type = tag
        walls_to_create = []
        for d in dead_ends:
            current_cell = d
            next_cell = self.get_neighbours(d, kind='accessible')
            while len(next_cell) == 1:
                current_cell.type = tag
                current_cell.solve_visited = True
                old_cell = current_cell
                current_cell = next_cell[0]
                if next_cell[0].type == 'start' or next_cell[0].type == 'end':
                    break
                next_cell = self.get_neighbours(current_cell, kind='accessible')
                next_cell = [c for c in next_cell if c != old_cell]
            walls_to_create.append((old_cell, old_cell.get_common_wall(current_cell)))
            walls_to_create.append((current_cell, current_cell.get_common_wall(old_cell)))
        return walls_to_create
    
    def break_dead_ends(self):
        """ Break dead ends in the maze. FUNCTION UNDER DEVELOPMENT: when braid_factor=1, 
        the maze might still have a tiny amount of dead ends... So far, this function 
        iterates over the dead ends, and choose randomly in a weighted list of boolean 
        if it has to break the dead end or not"""
        if self.braid_factor == 0:
            return
        else:
            # Creating a list to weight the number of dead end, according to braid_factor
            frac = Fraction.from_float(self.braid_factor)
            frac = frac.limit_denominator(10)
            rate = frac.denominator*[False]
            for ind in range(frac.numerator):
                rate[ind] = True
#            dead_ends = random.shuffle(self.get_dead_ends())
            dead_ends = self.get_dead_ends()
            for c in dead_ends:
                if c.isDeadEnd() and random.choice(rate):
                    n = self.get_neighbours(c, 'accessible')[0]
                    if n.i == c.i:
                        i2 = c.i
                        j2 = 2*c.j - n.j
                    elif n.j == c.j:
                        j2 = c.j
                        i2 = 2*c.i - n.i
                    try:
                        Cell.break_wall(c, self.maze_map[i2][j2])
                    except (ValueError, IndexError):
                        pass
                
            
    def get_neighbours(self, cell, kind='all'):
        """ 
        Returns the unvisited neighbours of a cell
        
        Parameters
        ----------
        
        cell: Cell object
            cell of interest
        kind: str
            type of neighbours to return. Possible values: all, unvisited (neighbours 
            that have all their walls), visited (neighbours that have at least one 
            broken wall), accessible (neighbours that has no common wall with cell)
            Default: all
            
        Returns
        ----------
        
        list of cell object corresponding to unvisited neighbours
        """
        delta = [(-1,0), (1,0), (0,1), (0,-1)]
        neighbours = []
        if kind == 'accessible':
            pair = {'top':(-1,0), 'bottom':(1,0), 'left':(0,-1), 'right':(0,1)}
            for k, v in cell.walls.items():
                if not v:
                    neighbours.append(self.maze_map[cell.i + pair[k][0]][cell.j + pair[k][1]])
            return neighbours
        for di, dj in delta:
            i2, j2 = cell.i + di, cell.j + dj
            if (0 <= i2 < self.p) and (0 <= j2 < self.q):
                neighbour = self.maze_map[i2][j2]
                if kind == 'all':
                    neighbours.append(neighbour)
                elif kind == 'unvisited':
                    if all(neighbour.walls.values()):
                        neighbours.append(neighbour)
                elif kind == 'visited':
                    if not all(neighbour.walls.values()):
                        neighbours.append(neighbour)
                elif kind == 'accessible':
                    pass
                else:
                    raise ValueError('Unknown kind of neighbour')
        return neighbours
    
    def get_image(self):
        """ Returns an image of the maze """
        im = np.ones((10*self.p + 1, 10*self.q + 1, 3))
        for i in range(self.p):
            for j in range(self.q):
                if self.maze_map[i][j].walls['top']:
                    im[10*i, 10*j:(10*(j + 1) + 1), :] = 0
                if self.maze_map[i][j].walls['left']:
                    im[10*i:(10*(i + 1) + 1), 10*j, :] = 0
                if self.maze_map[i][j].type == 'direct_dead_end':
                    im[(10*i + 1):10*(i + 1), (10*j + 1):10*(j + 1), 1:] = 0
                if self.maze_map[i][j].type == 'indirect_dead_end':
                    im[(10*i + 1):10*(i + 1), (10*j + 1):10*(j + 1), :] = 0.5
        im[10*self.p, :, :] = 0
        im[:, 10*self.q, :] = 0
        if hasattr(self, 'start'):
            istart = self.start[0]
            jstart = self.start[1]
            im[(10*istart + 1):10*(istart + 1), (10*jstart + 1):10*(jstart + 1), :2] = 0
        if hasattr(self, 'end'):
            iend = self.end[0]
            jend = self.end[1]
            im[(10*iend + 1):10*(iend + 1), (10*jend + 1):10*(jend + 1), ::2] = 0
        return im
    
    def plot(self):
        """ Dsplay the maze in a matplotlib figure """
        fig, ax = plt.subplots()
        ax.imshow(self.get_image())
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
    m = Maze(30, 30, method='dfs', braid_factor=0.5)
    print(m)
    m.solve(start=(0, 0), end=(m.p - 1, m.q - 1))
#    m.solve()
    m.plot()
