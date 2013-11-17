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

using namespace std;

int* generate_random_array(int grid_size);
int** convertArraytoMatrix(int* array, int n);
void solveMatrix(int* array, int grid_size);
void solveNextRow(int* array, int grid_size);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattanDistances(int* array, int grid_size);
void swap(int* array, int loc1, int loc2);
bool isSolution(int* array, int grid_size);
void relocate(int* array, int grid_size, int index, int piece);
bool checkLeft(int position, int n);
bool checkUp(int position, int n);
int findSmallest(int* array);
int* findNumLoc(int* array, int grid_size, int target);

void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int grid = pow(3,2);
int counter, last(grid-1);
int loc[2] = {0, 0};
int iterate = 0;
int prev = INT_MAX;

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
  //cout << it->second << " ";

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
    solveNextRow(array, grid_size);
  } else {
    cout <<"\nPuzzle is not solvable :/" << endl;
  }
}

void solveNextRow(int* array, int grid_size) {
  for (int i=0; i<sqrt(grid_size); i++) {
    int* position = findNumLoc(array, grid_size, i);
    int index = position[0]*sqrt(grid_size) + position[1];
    while (index != i) {
      relocate(array, grid_size, index, i);
      index = position[0]*sqrt(grid_size) + position[1];
      break;
    }
  }
}

void relocate(int* array, int grid_size, int index, int piece) {
  int* possible = new int[2];
  int location = loc[0]*sqrt(grid_size) + loc[1];
  if (checkLeft(index-1, sqrt(grid_size))) {
    possible[0] = abs(location - index - 1);
  }
  if (checkUp(index-1, sqrt(grid_size))) {
    possible[1] = abs(index-sqrt(grid_size)*(index / sqrt(grid_size))-sqrt(grid_size)+abs(location-index));
  }

  int shortest = findSmallest(possible);
  
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

int* findNumLoc(int* array, int grid_size, int target) {
  int* pos = new int[2];
  for (int i=0; i<grid_size; i++) {
    if (array[i] = target) {
      pos[0] = i / int(sqrt(grid_size));
      pos[1] = i % int(sqrt(grid_size));
    }
  }
  return pos;
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
    if (array[i] != i) {
      return false;
    }
  }
  return true;
}

void swap(int* array, int loc1, int loc2) {
  int temp = array[loc2];
  array[loc2] = array[loc1];
  array[loc1] = temp;
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
