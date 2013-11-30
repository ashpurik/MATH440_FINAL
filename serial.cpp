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

int* generate_random_array(int grid_size);
int** convertArraytoMatrix(int* array, int n);
void solveMatrix(int* array, int grid_size);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattanDistances(int* array, int grid_size);
void swap(int* array, int loc1, int loc2);

bool dfs(int steps, int currManhat, int* array, int grid_size, int bound);
int ida_star(int* array, int grid_size);
int manhattan(int start, int end, int* array, int n);
int convertCoordinatestoIndex(int i, int j, int n);
bool checkMove(int i, int j, int n);
bool checkLeft(int j);
bool checkUp(int i);
bool checkDown(int i, int n);
bool checkRight(int j, int n);

void printPath(int previous);
bool isSolution(int* array, int grid_size);

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
  
  clock_t start = clock();
  int N = 3;
  int GRID = N*N;

  cout << "\n---------------------------------" << endl;

  int* numbers = new int[GRID];
  int stuff[9] = {8, 0, 4, 7, 2, 6, 3, 5, 1};
  //numbers = generate_random_array(GRID);
   for (int i=0; i<GRID; i++) {
     numbers[i] = stuff[i];
   }
  //int* numbers = {2, 4, 6, 7, 3, 5, 8, 1, 0};
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

int* generate_random_array(int grid_size) {

  int *nums = new int[grid_size];
  for (int k=0; k<grid_size; k++) {
    nums[k] = k;
  }
  srand(time(0));
  random_shuffle(nums, nums+grid_size);
  return nums;
}

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

bool isSolvable(int* array, int grid_size) {
  std::map<int, int> num_loc_pairs;
  for (int i=0; i<grid_size; i++) {
    num_loc_pairs.insert(std::pair<int, int>(i, array[i]));
  }
  std::map<int,int>::iterator it;
  for (it=num_loc_pairs.begin(); it!= num_loc_pairs.end(); it++) {
    if (!(std::find(searchedElements.begin(), searchedElements.end(), it->first)!= searchedElements.end()) ) {
      findCycle(num_loc_pairs, it->first, it->first, sqrt(grid_size));
    }
  }

  int manhattan = 2*(sqrt(grid_size) -1) - loc[0] - loc[1];

  if ((manhattan+counter) % 2 == 0) {
    return true;
  } else {
    return false;
  }

}

int convertCoordinatestoIndex(int i, int j, int n) {
  return i*n+j;
}

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

void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n) {
  std::map<int, int>::iterator it;
  it = theMap.find(curr_key);
  searchedElements.push_back(it->second);

  if (it->second == ((n*n)-1)) {
    loc[0] = it->first / n;
    loc[1] = it->first % n;
  }

  if (initial_key == it->second) {
    return;
  } else {
    counter++;
    findCycle(theMap, initial_key, it->second, n);
  }
}

void solveMatrix(int* array, int grid_size) {
  bool solvable = isSolvable(array, grid_size);
  if (solvable) {
    int num_steps = ida_star(array, grid_size);
    printPath(num_steps), printf("\n");

  } else {
    cout <<"\nPuzzle is not solvable :/" << endl;
  }
}

int ida_star(int* array, int grid_size) {
  int bound = manhattanDistances(array, grid_size);
  cout << "INITIAL BOUND: " << bound << endl;
  while (true) {
    nextBound = INT_MAX; // next limit
    path.clear();
    visited.clear();
    if (dfs(0, manhattanDistances(array, grid_size), array, grid_size, bound)) {
      cout << "\nFINAL NUMBER OF STEPS: " << bound << endl;
      cout << "VISITED SIZE:" << visited.size() << endl;
      return bound;
    }
    //cout << "VISITED SIZE:" << visited.size() << endl;
    cout << "NEW BOUND: " << bound << endl;
    bound = nextBound; // nextBound > bound
  }
}

bool dfs(int steps, int currManhat, int* array, int grid_size, int bound) {

  if (steps + currManhat > bound) {
    nextBound = min(nextBound, steps + currManhat);
    return false;
  }

  if (isSolution(array, grid_size)) {
    return true;
  }

  unsigned long long state = 0;
  for (int i = 0; i < grid_size; i++) { // transform the grid into an unsigned long long, to store the state
    state <<= int(sqrt(grid_size)); // move to the left n size bits
    state += array[i]; // add it to the state
  }

  if (visited.count(state) && (visited[state] <= steps)) { // we want to prevent cycling on the grid
    //cout << "cycling" << endl;
    return false; 
  }

  visited[state] = steps;

  int i, new_i, new_j;
  for (i = 0; i < grid_size; i++) {
    if (array[i] == grid_size-1) {
      break;
    }
  }

  int j = i % int(sqrt(grid_size));
  i /= sqrt(grid_size);

  int index = convertCoordinatestoIndex(i, j, sqrt(grid_size));

  if(checkRight(j+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i, j+1, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1);
    path[steps + 1] = 'R';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) { // if ok, no need to restore, just go ahead
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1); // restore
  }
  if (checkUp(i-1)) {
    int new_index = convertCoordinatestoIndex(i-1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j);
    path[steps + 1] = 'U';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) { // if ok, no need to restore, just go ahead
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j); // restore
  }
  if (checkLeft(j-1)) {
    int new_index = convertCoordinatestoIndex(i, j-1, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); //swap
    path[steps + 1] = 'L';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) { // if ok, no need to restore, just go ahead
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); // restore
  }
  if (checkDown(i+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i+1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j);
    path[steps + 1] = 'D';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound)) { // if ok, no need to restore, just go ahead
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j); // restore
  }

  return false;
}

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

bool checkMove(int i, int j, int n) {
  return (checkLeft(j) && checkRight(j, n) && checkUp(i) && checkDown(i, n));
}

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
