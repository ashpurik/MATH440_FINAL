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

using namespace std;

int* generate_random_array(int grid_size);
int** convertArraytoMatrix(int* array, int n);
void solveMatrix(int* array, int grid_size);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattanDistances(int* array, int grid_size);
void swap(int* array, int loc1, int loc2);
bool isSolution(int* array, int grid_size);
bool checkLeft(int position, int n);
bool checkUp(int position, int n);
int findSmallest(int* array);

inline int h2(int i1, int j1, int i2, int j2, int* array);
int ida_star(int* array, int grid_size);
bool dfs(int g, int h, int* array, int grid_size);
inline bool valid(int r, int c);
//inline void swap(int i, int j, int new_i, int new_j, int* array);

void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int grid = pow(3,2);
int counter, last(grid-1);
int loc[2] = {0, 0};
int iterate = 0;
map<int, int> predecessors;
map<unsigned long long, int> visited;
int dr[] = { 0,-1, 0, 1}; // E,N,W,S
int dc[] = { 1, 0,-1, 0}; // R,U,L,D
int bound, nextBound;

int main ( int argc, char *argv[] ) {
  
  int* numbers = generate_random_array(grid);
  int** matrix = convertArraytoMatrix(numbers, sqrt(grid));
  printMatrix(matrix, sqrt(grid), sqrt(grid));
  cout << endl;
  //printArray(numbers, grid);
  //cout << endl;
  solveMatrix(numbers, grid);  
  //int** solved = convertArraytoMatrix(numbers, sqrt(grid));
  //printMatrix(solved, sqrt(grid), sqrt(grid));

  delete numbers;
  freeMatrix(matrix, sqrt(grid));
  //freeMatrix(solved, sqrt(grid));
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
      //cout << "( ";
      findCycle(num_loc_pairs, it->first, it->first, sqrt(grid_size));
      //cout << ")";
    }
  }

  //cout << "\n\nCounter = " << counter;
  int manhattan = 2*(sqrt(grid_size) -1) - loc[0] - loc[1];
  //cout << "\n\nManhattan Distance = " << manhattan;
  //cout << "\n\nNumber of Interchanges = " << manhattan+counter << endl;

  if ((manhattan+counter) % 2 == 0) {
    return true;
  } else {
    return false;
  }

}

void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n) {
  std::map<int, int>::iterator it;
  it = theMap.find(curr_key);
  searchedElements.push_back(it->second);

  if (it->second == last) {
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
    int t = ida_star(array, grid_size);
    cout << "\nBound: " << t << endl;
  } else {
    cout <<"\nPuzzle is not solvable :/" << endl;
  }
}

int ida_star(int* array, int grid_size) {
  bound = manhattanDistances(array, grid_size);
  while(true) {
    //cout << "in while" << endl;
    nextBound = INT_MAX;
    predecessors.clear();
    visited.clear();
    if (dfs(0, manhattanDistances(array, grid_size), array, grid_size)) {
      return bound;
    }
    if (nextBound == INT_MAX) {
      cout << "\nNextBound is infinity" << endl;
      return -1;
    }
    bound = nextBound;
    if (bound > 45) {
      cout << "more than 45" << endl;
      return -1;
    }
  }
}

bool dfs(int g, int h, int* array, int grid_size) {

  if (iterate == 3) {
    return false;
  }

  iterate++;
  cout << "\nDFS, Iteration: " << iterate << endl;

  if (g + h > bound) {
    cout << "in first if" << endl;
    nextBound = min(nextBound, g + h);
    return false;
  }

  if (isSolution(array, grid_size)) {
    cout << "solution found" << endl;
    return true;
  }
  
  unsigned long long state = 0;
  for (int i=0; i<grid_size; i++) {
    state <<= 4;
    state += array[i];
  }

  if(visited.count(state) && visited[state] <= g) {
    return false;
  }
  visited[state] = g;

  int i;
  for (i=0; i<grid_size; i++) {
    if (array[i] == grid_size-1) {
      break;
    }
  }

  int j = i % int(sqrt(grid_size));
  i /= sqrt(grid_size);

  cout << "Blank: (" << i << ", " << j << ")" << endl;

  int new_i, new_j;
  for (int d = 0; d < sqrt(grid_size)+1; d++) {
    new_i = i + dr[d]; new_j = j = dc[d];
    cout << "New possible Spot: (" << new_i << ", " << new_j << ")" << endl;
    if (valid(new_i, new_j)) {
      cout << "New Valid Spot: (" << new_i << ", " << new_j << ")" << endl;
      int dh = h2(i, j, new_i, new_j, array);
      //swap(i, j, new_i, new_j, array);
      swap(array, i*sqrt(grid_size)+j, new_i*sqrt(grid_size)+new_j);

      predecessors[g+1] = d;

       cout << "swap occured" << endl;
       int** test = convertArraytoMatrix(array, sqrt(grid_size));
       printMatrix(test, sqrt(grid_size), sqrt(grid_size));

      if (dfs(g+1, h+dh, array, grid_size)) {
        return true;
      }
      cout << "after recursive call swap occured " << endl;
      swap(array, i*sqrt(grid_size)+j, new_i*sqrt(grid_size)+new_j);
      //swap(i, j, new_i, new_j, array);
    }
  }
  return false;
}

int findSmallest(int* array) {
  int smallest = INT_MAX;
  int index = 0;
  for (int i=0; i<2; i++) {
    if (array[i] < smallest) {
      smallest = array[i];
      index = i;
    }
  }
  return index;
}

inline int h2(int i1, int j1, int i2, int j2, int* array) { // heuristic: sum of manhattan distances (compute delta)
  int tgt_i = array[i2 * 4 + j2] / 4, tgt_j = array[i2 * 4 + j2] % 4;
  return -(abs(i2 - tgt_i) + abs(j2 - tgt_j)) + (abs(i1 - tgt_i) + abs(j1 - tgt_j));
}

inline bool valid(int r, int c) {
  return 0 <= r && r < 4 && 0 <= c && c < 4;
}

bool checkLeft(int position, int n) {
  //cout << "Left check: " << position << ". Previous: " << prev << endl;
  if (position % n == 0) {
    return false;
  } else {
    return true;
  }
}

bool checkUp(int position, int n) {
  //cout << "Up check: " << position << ". Previous: " << prev << endl;
  if (position / n == 0) {
    return false;
  } else {
    return true;
  }
}

bool isSolution(int* array, int grid_size) {
  for (int i=0; i<grid_size; i++) {
    if (array[i] != grid_size && array[i] != i) {
      return false;
    }
  }
  return true;
}

void swap(int* array, int loc1, int loc2) {
  cout << "loc1 = " << loc1 << ". loc2 = " << loc2 << endl;
  int temp = array[loc1];
  array[loc1] = array[loc2];
  //cout << "IN SWAP: Loc1 now contains " << array[loc1] << endl;
  array[loc2] = temp;
  //cout << "IN SWAP: Loc2 now contains " << array[loc2] << endl;
}


//inline void swap(int i, int j, int new_i, int new_j, int* array) {
//  int temp = array[i * 4 + j];
//  array[i * 4 + j] = array[new_i * 4 + new_j];
//  array[new_i * 4 + new_j] = temp;
//}

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
