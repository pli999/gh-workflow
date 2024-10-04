from itertools import combinations, permutations
from typing import List
import math
from tqdm import tqdm
import numpy as np


def expected_length(field: List[str], K: int) -> float:

    # find all checkpoints
    cp = list()
    m = len(field[0])
    n = len(field)
    for i in range(n):
        for j in range(m):
            if field[i][j] == '*':
                cp.append([i,j])

    # get all combinations of k checkpoints from all cp
    # e.g. [[c1,c2],[c1,c3][c2,c3]]
    total_sp = 0

    memo = [[[[-1 for _ in range(m)] for _ in range(n)] for _ in range(m)] for _ in range(n)]

    # calculate shortestPath between any two checkpoints
    cp_size = len(cp)
    shortest = [[math.inf for _ in range(cp_size)] for _ in range(cp_size)]


    print(f"Calculating shortest path of {cp_size} checkpoints...")
    filename = f"sp_{n}_{m}_{cp_size}.npy"    
    
    try:
        memo = np.load(filename)
        print(f"File {filename}. loaded")
    except OSError:

        for i in tqdm(range(cp_size)):
            for j in range(i+1,cp_size):
                shortest[i][j] = shortest[j][i] = shortestPath(cp[i][0], cp[i][1], cp[j][0], cp[j][1],n,m,field, memo)
        np.save(filename, memo)
        print(f"File {filename}. saved")

    subset = list(combinations(cp, K))
    subset_size = len(subset)
    print(f"Calculating expected length of {cp_size}C{K}={subset_size} combinations...")
    for x in tqdm(range(subset_size)):
        total_sp += sp_set3(subset[x],n,m,field,memo)

    return total_sp/len(subset)
        


def sp_set(cpSet,n,m,field,memo):
    permu = list(permutations(cpSet))
    min_total = math.inf
    for i in range(len(permu)//2):  # only need to calculate half of the permutation as A-B-C and C-B-A have same distance
        sp = 0
        for j in range(K-1):
            sp += shortestPath(permu[i][j][0], permu[i][j][1], permu[i][j+1][0], permu[i][j+1][1],n,m,field, memo)
        min_total = min(min_total, sp)
    
    return min_total
    

def sp_set2(cpSet,n,m,field,memo):

    l_cpSet = len(cpSet)

    if l_cpSet == 2:
        return shortestPath(cpSet[0][0], cpSet[0][1], cpSet[1][0], cpSet[1][1],n,m,field, memo)

    min_cycle = 0
    min_sequence = []

    # shortest cycle of the first 3 points (a triangle)
    for i in range(3):
        min_cycle += shortestPath(cpSet[i][0], cpSet[i][1], cpSet[(i+1)%3][0], cpSet[(i+1)%3][1],n,m,field, memo)

    min_sequence = [0,1,2]

    # insert remaining points to the cycles
    for k in range(3, l_cpSet):

        insert_after = -1
        insert_min = math.inf

        for i in range(k):
            new_cycle = min_cycle + \
                        shortestPath(cpSet[k][0], cpSet[k][1], cpSet[min_sequence[i]][0], cpSet[min_sequence[i]][1],n,m,field,memo) + \
                        shortestPath(cpSet[k][0], cpSet[k][1], cpSet[min_sequence[(i+1)%k]][0], cpSet[min_sequence[(i+1)%k]][1],n,m,field,memo) - \
                        shortestPath(cpSet[min_sequence[i]][0], cpSet[min_sequence[i]][1], cpSet[min_sequence[(i+1)%k]][0], cpSet[min_sequence[(i+1)%k]][1],n,m,field,memo)
            if new_cycle < insert_min:
                insert_min = new_cycle
                insert_after = i
        
        assert(insert_after != -1)

        min_cycle = insert_min
        min_sequence = min_sequence[0:insert_after+1] + [k] + min_sequence[insert_after+1:]

    # cut cycle at the longest edge
    assert(l_cpSet == len(min_sequence))
    longest_edge = 0
    break_after = 0

    for i in range(len(min_sequence)):

        edge = shortestPath(cpSet[min_sequence[i]][0], cpSet[min_sequence[i]][1], cpSet[min_sequence[(i+1)%l_cpSet]][0], cpSet[min_sequence[(i+1)%l_cpSet]][1],n,m,field,memo)
        if edge > longest_edge:
            longest_edge = edge
            break_after = i
    
    return min_cycle - longest_edge
        
def sp_set3(cpSet,n,m,field,memo):

    l_cpSet = len(cpSet)

    min_sequence = [0,1]
    min_sp = shortestPath(cpSet[0][0], cpSet[0][1], cpSet[1][0], cpSet[1][1],n,m,field, memo)

    if l_cpSet == 2:
        return min_sp
    
    # insert remaining points to the cycles
    for k in range(2, l_cpSet):

        insert_before = -1
        insert_min = math.inf

        for i in range(k+1):
            new_sp = min_sp
            if i>0:
                new_sp += shortestPath(cpSet[k][0], cpSet[k][1], cpSet[min_sequence[i-1]][0], cpSet[min_sequence[i-1]][1],n,m,field,memo)
            if i<k:
                new_sp += shortestPath(cpSet[k][0], cpSet[k][1], cpSet[min_sequence[i]][0], cpSet[min_sequence[i]][1],n,m,field,memo)
            if i>0 and i<k:
                new_sp -= shortestPath(cpSet[min_sequence[i-1]][0], cpSet[min_sequence[i-1]][1], cpSet[min_sequence[i]][0], cpSet[min_sequence[i]][1],n,m,field,memo)
            if new_sp < insert_min:
                insert_min = new_sp
                insert_before = i
        
        assert(insert_before != -1)
        
        min_sp = insert_min
        min_sequence = min_sequence[0:insert_before] + [k] + min_sequence[insert_before:]

    return min_sp


def shortestPath(s_i,s_j,e_i,e_j,n,m,field,memo):
   
    if memo[s_i][s_j][e_i][e_j] != -1:
        return memo[s_i][s_j][e_i][e_j]

    visited = [[False]*m for _ in range(n)]
    adj = [[0,-1],[-1,0],[0,1],[1,0]]
     
    memo[s_i][s_j][e_i][e_j] = shortestPathHelper(s_i,s_j,e_i,e_j,n,m,field,adj,visited,memo)
    memo[e_i][e_j][s_i][s_j] = memo[s_i][s_j][e_i][e_j]  # distance should be symmetric

    return memo[s_i][s_j][e_i][e_j]


def shortestPathHelper(s_i,s_j,e_i,e_j,n,m,field,adj,visited,memo):

    min_path = math.inf

    #depth first search of the field
    visited[s_i][s_j] = True
    if s_i==e_i and s_j==e_j:
        min_path = 0
    else:
        # try to move to all adjacent cells
        for x in range(4):
            new_i = s_i+adj[x][0]
            new_j = s_j+adj[x][1]

            if(isSafe(new_i, new_j,n,m,field) and not visited[new_i][new_j]):
                min_path=min(min_path, 1+shortestPathHelper(new_i,new_j,e_i,e_j,n,m,field,adj,visited,memo))

    visited[s_i][s_j] = False

    memo[s_i][s_j][e_i][e_j] = min(min_path,  memo[s_i][s_j][e_i][e_j])

    return min_path



# check if i,j position is valid
def isSafe(i,j,n,m,field):
    return i>=0 and i<n and j>=0 and j<m and field[i][j] != '#'





field = [\
"*#..#",\
".#*#.",\
"*...*"]
K = 2
print(3.833333, expected_length(field,K))




field = [
 "*#..#",
 ".#*#.",
 "*...*"]
K = 4

print("Case 2\r", field)
print(8.0, expected_length(field,K))


field = [
 "#.#**",
 "....#",
 "#*#**",
 "**#*#",
 "#..##",
 "*#..#",
 ".#.#.",
 "....*"]
K = 3

print("Case 3\r", field)
print(10.825000000000024, expected_length(field,K))



field = [    	
 "###################",
 "#*###############*#",
 "#.....#######.....#",
 "#*###*.#.*.#.*###*#",
 "#*####*.*#*.*####*#",
 "#*#####*###*#####*#",
 "###################"]
K = 9

print("Case 4\r", field)
print(30.272233648704244, expected_length(field,K))



field = [    	
 "**##*.**#..#.*...*#...*#..#.##..#..#.#*...#.##*##.",
 ".#..###..#..#.#.##..#.#.*#.*..#..#.#*..##.#*...*..",
 "..#.....###.#*.##..#.#.#*..#.#..#....#..#...#*####",
 ".#.##*#.*#..#*#*.#.#...*.#.*#.#.##.#*.##.#.#..*...",
 "..*.*#*.###.#..#.#..##.##.*#..#.....#.....#..#.#.#",
 ".#.##.#..##..*#..#.#...#*##*#*..#.#.#.#.##.##.#.#*",
 "..##....#..#.#*#...*.##...#.#.####...#.#*.....#...",
 ".#.*#.##.*#*.#*.#.#.#..#.#..#.#*#.###..##.##.#.##*",
 ".*.#*..*.#.#...#.*##.#.**.#.*...**..*#..#.#.#*.#..",
 ".#*.#*##....##.#.#*..*.###.#.##.##.#.#.#....#.#*.#",
 "*.#..#*#.#*#*....#.#.#..*#**...##.#.#.**#*##.*.#..",
 ".#*.##..##..##.#.#..#.#.###.###...#...#*#..##*#.#.",
 "#..#*.#..*.###..#.#...#.###.#.#*#.#.#**##.#...*.#*",
 "..#..#.#.##.#..#.**.##*#.#**.**..#.#..#...#.##*#..",
 ".#*#.#.*..#.*#...#.#...#...#.##.#..*#*.##*....###.",
 ".*.#.#.#.#*#..*##.**.##*##..#.*#.#*###..*.#.##.#..",
 ".#......#...#.#.*#.#.#..#..#.#*#....#*.#*#.*#..*.#",
 "#..####..#*#...#*.#..#.###...#.#.#.###*#..##*##.#.",
 ".#.*..#.#...#.#..#.##...#..#.#.#.#.###..##..*.*.*.",
 ".#.#.#.#..##.*..#.*.#.##.#..##*...#.#..#.#.##.#.##",
 ".#..#*.#.#..#.##..##..#.*..#.*#.#...##....#...###.",
 ".#.#.#.#*.#.#..#.#..#..#.#.*#...#.##...#.##.##.*..",
 ".#...#.#.##.#.#..*#.*#..###..#.#.#*###.##...#*.##.",
 ".#.##.*.......*.#.*#.#.#*###..*...*..#.*.##.#.#..#",
 "...###*####*#.#..##*...#..#..##.#.#.#..##*#*.*.*#.",
 "#.#.#....*#..#.#.#.#.##..#*.#...#..#.#*#...#.##.*.",
 "..*.#*##.#.#*#.###...#..##.#.#.#*###*#.*#.#.*###.#",
 "##*##..##...#.....##.#.#.**#..#*.....##.#..#*.#.*.",
 ".....#.*.##..##.##*.*#...#.#.#.##.#*#.**..#..#.#.#",
 "##.#.#*##.#.#.*.*.#.#*#.#.#....*...#*##*##.#....#.",
 "*.**#**....*..##.#*.*.**..##.###.##.....##...##.**",
 "#.####.##*#*##..#.*#*#.##*...#.##..#.##....#*..##.",
 "....#...##.#...#*.#..##.##.#*..*.#....##.#.*##...#",
 "#.#..*##*..#.#..#..#..#*....#.##..##.#*##.##.*##..",
 "..#.#*.*.##.#.#*#.#*##.###.##...#............#*.#.",
 "#.#.##.#....*....*..##..*#.#.#.###.#.#.#.###..#..#",
 ".#**..#*#.#*#*#.#.#...*##....##.#*..#..#*..*#..#..",
 "...#*#.....#..#.#..#*#.*##.#..#.#.##..#.*#*#.#...#",
 ".#*.###.#.#.#.#.*#*##.##..#.#*..#...#.#.#..#*.*#..",
 "#*.#.#.#..#..#..#....*#.*##..##.#.#..#...##.#.#..#",
 "*.#..#..#...#..##.#*#..#.#*#.#.#.###..#.#*...#.#..",
 "#...#.#...#.#.#..#.*.#*.....**.*..#*##.#*.##....##",
 "#*#....#*#..#.*.###*#..#*##.##.#.#...#.*.##.##.##.",
 "..##*##*..#*#.#..#*.*##*.##.#...#.#.#.#.#..*#.##..",
 "#...#*##.#*#**.##.*#.*.##..*.#*#**....#**##...*.*#",
 "*#.##......*#.##.#.#.##**.#.#.#.#.#.##..#...#*#*#*",
 "*....##.#.#..#.....#..##.#....*....#.#.##.#.#.##**",
 "#.##*#...#..#.#.##..#..##.##.##.##........##.#*#.#",
 "..#...#.#*#*..*#..*#.*#.#......##.#.#.#*#..#..****",
 ".###.#..#...#.#..#..#.#...#.#.#...**.#..*#*.*##*#."]
K = 150

print("Case 5\r", field)
print(1309.4951033725558, expected_length(field,K))
