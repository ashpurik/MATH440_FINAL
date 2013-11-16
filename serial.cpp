# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <stdlib.h>
# include <algorithm>
# include <cmath>
# include <map>
# include <vector>
# include <time.h>

using namespace std;

int* generate_random_array(int grid_size);
int** convertArraytoMatrix(int* array, int n);
void solveMatrix(int* array, int n);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);

void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int counter, last(8);
int loc[2] = {0, 0};

int main ( int argc, char *argv[] ) {
  
  int* numbers = generate_random_array(9);
  int** matrix = convertArraytoMatrix(numbers, 3);
  printMatrix(matrix, 3, 3);
  cout << endl;
  printArray(numbers, 9);
  cout << endl;
  solveMatrix(numbers, 9);  

  delete numbers;
  freeMatrix(matrix, 3);
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
      cout << "( ";
      findCycle(num_loc_pairs, it->first, it->first, sqrt(grid_size));
      cout << ")";
    }
  }

  cout << "\n\nCounter = " << counter;
  int manhattan = 2*(sqrt(grid_size) -1) - loc[0] - loc[1];
  cout << "\n\nManhattan Distance = " << manhattan;
  cout << "\n\nNumber of Interchanges = " << manhattan+counter << endl;

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
  cout << it->second << " ";

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

void solveMatrix(int* array, int n) {
  bool solvable = isSolvable(array, n);
  if (solvable) {
    cout << "\nPuzzle is solvable!" << endl;
  } else {
    cout <<"\nPuzzle is not solvable :/" << endl;
  }
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
