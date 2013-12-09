# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <stdlib.h>
# include <algorithm>
# include <cmath>
# include <map>
# include <vector>
# include <time.h>
# include <climits>
# include <cstdio>
# include <ctime>

using namespace std;

//FUNCTIONS TO SOLVE PUZZLE
int* generate_random_array(int grid_size);
void solveMatrix(int* array, int grid_size);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattanDistances(int* array, int grid_size);
void swap(int* array, int loc1, int loc2);
bool dfs(int steps, int currManhat, int* array, int grid_size, int bound);
int ida_star(int* array, int grid_size);
int manhattan(int start, int end, int* array, int n);
bool checkLeft(int j);
bool checkUp(int i);
bool checkDown(int i, int n);
bool checkRight(int j, int n);
void printPath(int previous);
bool isSolution(int* array, int grid_size);

//OTHER FUNCTIONS
int** convertArraytoMatrix(int* array, int n);
int convertCoordinatestoIndex(int i, int j, int n);
void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int counter;
int loc[2] = {0, 0};
int nextBound;
map<char, char> path;
map<unsigned long long, int> visited;

int main ( int argc, char *argv[] ) {

  int N;

  //grabs value of N from command line
  if (argc < 2) {
    cout << "Please specify N" << endl;
  } else {
   N = atoi( argv[1]);
  }
  
  clock_t start = clock();
  int GRID = N*N;

  cout << "\n---------------------------------" << endl;

  int* numbers = new int[GRID];
  numbers = generate_random_array(GRID);
  int** matrix = convertArraytoMatrix(numbers, N);
  printMatrix(matrix, N, N);
  solveMatrix(numbers, GRID);
  int** solved = convertArraytoMatrix(numbers, N);
  cout << endl;
  printMatrix(solved, N, N);

  delete numbers;
  freeMatrix(matrix, N);
  freeMatrix(solved, N);
  clock_t end = clock();

  cout << "---------------------------------" << endl;

  cout << setprecision(15) << "EXECUTION TIME: " << double(end-start)/CLOCKS_PER_SEC << " seconds\n" << endl;

  return 0;
}

//creates a random puzzle
int* generate_random_array(int grid_size) {

  int *nums = new int[grid_size];
  //creates a list of N*N numbers
  for (int k=0; k<grid_size; k++) {
    nums[k] = k;
  }
  srand(time(0));
  //then shuffles them up to be random
  random_shuffle(nums, nums+grid_size);
  return nums;
}

//for readability purposes, converts an array to matrix
int** convertArraytoMatrix(int* array, int n) {
  int index = 0;
  int** matrix = new int*[n];
  for (int i=0; i<n; i++) {
    matrix[i] = new int[n];
    for (int j=0; j<n; j++) {
      matrix[i][j] = array[index];
      index++;
    }
  }
  return matrix;
}

//checks whether a given arrangement of the puzzle is solvable
bool isSolvable(int* array, int grid_size) {
  std::map<int, int> num_loc_pairs;

  //creates a map with all keys being the numbers in correct order
  //and the values being the current random arrangement of numbers
  for (int i=0; i<grid_size; i++) {
    num_loc_pairs.insert(std::pair<int, int>(i, array[i]));
  }
  std::map<int,int>::iterator it;

  //locates all of the cycles in the map
  for (it=num_loc_pairs.begin(); it!= num_loc_pairs.end(); it++) {
    if (!(std::find(searchedElements.begin(), searchedElements.end(), it->first)!= searchedElements.end()) ) {
      findCycle(num_loc_pairs, it->first, it->first, sqrt(grid_size));
    }
  }

  //calculates the manhattan distance between where the current blank spot
  //is and where it should be (bottom right corner)
  int manhattan = 2*(sqrt(grid_size) -1) - loc[0] - loc[1];

  //if the sum of the manhattan distance and the number of cycles
  //is even, then the puzzle is solvable!
  if ((manhattan+counter) % 2 == 0) {
    return true;
  } else {
    return false;
  }

}

//converts matrix coordinate pairs to an array index
int convertCoordinatestoIndex(int i, int j, int n) {
  return i*n+j;
}

//recursively finds all cycles and their lengths in the grid
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n) {
  std::map<int, int>::iterator it;
  it = theMap.find(curr_key);
  searchedElements.push_back(it->second);

  //if blank spot is located, record its position
  if (it->second == ((n*n)-1)) {
    loc[0] = it->first / n;
    loc[1] = it->first % n;
  }

  //base case: if a cycle is found, exit recursion
  if (initial_key == it->second) {
    return;
  //otherwise record the length of the cycle so far and keep searching
  } else {
    counter++;
    findCycle(theMap, initial_key, it->second, n);
  }
}

//if a solvable arrangement is found, solve it!
//otherwise print out appropriate message
void solveMatrix(int* array, int grid_size) {
  bool solvable = isSolvable(array, grid_size);
  if (solvable) {
    int num_steps = ida_star(array, grid_size);
    printPath(num_steps), printf("\n");

  } else {
    cout <<"\nPuzzle is not solvable :/" << endl;
  }
}

//iterative deepening: keep increasing the bound until
//a solution is found
int ida_star(int* array, int grid_size) {
  int bound = manhattanDistances(array, grid_size);
  cout << "INITIAL BOUND: " << bound << endl;
  while (true) {

    //reset path, visited, and the nextBound everytime the bound is increased
    nextBound = INT_MAX; // next limit
    path.clear();
    visited.clear();
    if (dfs(0, manhattanDistances(array, grid_size), array, grid_size, bound)) {
      cout << "\nFINAL NUMBER OF STEPS: " << bound << endl;
      return bound;
    }
    cout << "NEW BOUND: " << bound << endl;
    bound = nextBound; // nextBound > bound
  }
}

//DEPTH-FIRST SEARCH: searches all possibilites within a given bound 
bool dfs(int steps, int currManhat, int* array, int grid_size, int bound) {

  //if more steps have been taken than the bound, increase it
  if (steps + currManhat > bound) {
    nextBound = min(nextBound, steps + currManhat);
    return false;
  }

  //if the solution is found, we are finished
  if (isSolution(array, grid_size)) {
    return true;
  }

  //store the state of the grid that is currently being explored
  unsigned long long state = 0;
  for (int i = 0; i < grid_size; i++) { // transform the grid into an unsigned long long, to store the state
    state <<= int(sqrt(grid_size)); // move to the left n size bits
    state += array[i]; // add it to the state
  }

  if (visited.count(state) && (visited[state] <= steps)) { // we want to prevent cycling on the grid
    return false; 
  }

  visited[state] = steps;

  //locate the blank spot of the current grid
  int i, new_i, new_j;
  for (i = 0; i < grid_size; i++) {
    if (array[i] == grid_size-1) {
      break;
    }
  }

  int j = i % int(sqrt(grid_size));
  i /= sqrt(grid_size);

  int index = convertCoordinatestoIndex(i, j, sqrt(grid_size));

  //check if the blank spot can move to the right
  if(checkRight(j+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i, j+1, sqrt(grid_size));
    //calculate if the new board will be further or closer away from goal
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    //perform the swap & record the step just taken
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1);
    path[steps + 1] = 'R';
    //recursively call DFS on this new board
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) {
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1); // restore
  }
  //check if the blank spot can move to upwards
  if (checkUp(i-1)) {
    int new_index = convertCoordinatestoIndex(i-1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j);
    path[steps + 1] = 'U';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) {
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j); // restore
  }
  //check if the blank spot can move to the left
  if (checkLeft(j-1)) {
    int new_index = convertCoordinatestoIndex(i, j-1, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1));
    path[steps + 1] = 'L';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) {
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); // restore
  }
  //check if the blank spot can move downwards
  if (checkDown(i+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i+1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j);
    path[steps + 1] = 'D';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) {
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j); // restore
  }

  return false;
}

//print out the path taken, working backwards in the order the path was saved
void printPath(int previous) {
  if (previous == 0) {
    return;
  }
  printPath(previous - 1);
  printf("%c", path[previous]);
}

bool checkLeft(int j) {
  return (j >= 0);
}

bool checkRight(int j, int n) {
  return (j < n);
}

bool checkUp(int i) {
  return (i >= 0);
}

bool checkDown(int i, int n) {
  return (i < n);
}

//check if solution has been found, ie all numbers in correct order
bool isSolution(int* array, int grid_size) {
  for (int i=0; i<grid_size; i++) {
    if (array[i] != grid_size-1 && array[i] != i) {
      return false;
    }
  }
  return true;
}

void swap(int* array, int loc1, int loc2) {
  int temp = array[loc1];
  array[loc1] = array[loc2];
  array[loc2] = temp;
}

//rather than calculating all of the manhattan distances again,
//calculates if the new swap performed is closer or further away
//from solution
int manhattan(int start, int end, int* array, int n) {
    int goal_i = array[end] / n; //where the number at possible swap should end up (i)
    int goal_j = array[end] % n; //where the number at possible swap should end up (j)
    int initial_i = start / n; //current i location of space
    int initial_j = start % n; //current j location of space
    int possible_i = end / n; //the position of i possible swap
    int possible_j = end % n; //the position of j possible swap
    int dist_blank_from_goal = (abs(initial_i - goal_i) + abs(initial_j - goal_j));
    int dist_possible_from_goal = (abs(possible_i - goal_i) + abs(possible_j - goal_j));
    return dist_blank_from_goal - dist_possible_from_goal;
}

//calculate the cumulative manhattan distances of the whole grid
int manhattanDistances(int* array, int grid_size) {
  int distance = 0;
  int* startPosition = new int[2];
  int* goalPosition = new int[2];
  for (int i=0; i<grid_size; i++) {
    startPosition[0] = i / sqrt(grid_size);
    startPosition[1] = i % int(sqrt(grid_size));
    goalPosition[0] = array[i] / sqrt(grid_size);
    goalPosition[1] = array[i] % int(sqrt(grid_size));
    distance += (abs(startPosition[0]-goalPosition[0])+abs(startPosition[1]-goalPosition[1]));
  }
  return distance;
}

void freeMatrix(int** matrix, int row) {

  for (int i=0; i<row; i++) {
    delete [] matrix[i];
  }
  delete matrix;
}

void printArray(int* array, int size) {

  for (int i=0; i<size; i++) {
    cout << array[i] << " ";
  }
  cout << endl;
}

void printMatrix(int** matrix, int row, int col) {

  for(int i=0; i<row; i++) {
    for(int j=0; j<col; j++) {
      cout << matrix[i][j] << " ";
    }
  cout << endl;
  }
}
