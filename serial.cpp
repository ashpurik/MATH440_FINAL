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
void solve(int* array, int grid_size);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattanDistances(int* array, int grid_size);

void swap(int* array, int loc1, int loc2);
bool checkLeft(int position, int n);
bool checkRight(int position, int n);
bool checkUp(int position, int n);
bool checkDown(int position, int n);
int findSmallest(int* array);

void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int grid = pow(10,2);
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
  int** solved = convertArraytoMatrix(numbers, sqrt(grid));
  printMatrix(solved, sqrt(grid), sqrt(grid));

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
    solve(array, grid_size);
  } else {
    cout <<"\nPuzzle is not solvable :/" << endl;
  }
}

void solve(int* array, int grid_size) {
  cout << "Iteration: " << iterate << endl;
  int location = loc[0]*sqrt(grid_size) + loc[1];
  //cout << "Location : " << location << endl;
  int distances[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};

  if (iterate == 100) {
    return;
  }

  if (manhattanDistances(array, grid_size) == 0) {
    return;
  } else {
    if(checkLeft(location, sqrt(grid_size))) {
      //cout << "in left" << endl;
      swap(array, location, location-1);
      distances[0] = manhattanDistances(array, grid_size);
      swap(array, location, location-1);
    }
    if(checkRight(location, sqrt(grid_size))) {
      //cout << "in right" << endl;
      swap(array, location, location+1);
      distances[1] = manhattanDistances(array, grid_size);
      swap(array, location, location+1);
    }
    if(checkUp(location, sqrt(grid_size))) {
      //cout << "in up" << endl;
      swap(array, location-1, location);
      distances[2] = manhattanDistances(array, grid_size);
      swap(array, location-1, location);
    }
    if(checkDown(location, sqrt(grid_size))) {
      //cout << "in down" << endl;
      swap(array, location+1, location);
      distances[3] = manhattanDistances(array, grid_size);
      swap(array, location+1, location);
    }

    printArray(distances, 4);
    cout << endl;
    //cout << "Previous: " << prev << endl;

    int smallestPos = findSmallest(distances);
    switch (smallestPos) {
      case 0: 
        //cout << "swapped Left\n" << endl;
        swap(array, location, location-1);
        loc[1] = (location-1) % int(sqrt(grid_size));
        prev = 0;
        break;
      case 1:
        //cout << "swapped Right\n" << endl;
        swap(array, location, location+1);
        loc[1] = (location+1) % int(sqrt(grid_size));
        prev = 1;
        break;
      case 2:
        //cout << "swapped Up\n" << endl;
        swap(array, location-sqrt(grid_size), location);
        loc[0] = (location-sqrt(grid_size)) / sqrt(grid_size);
        //cout << "Moving to " << location-1 << "from " << location << endl;
        prev = 2;
        break;
      case 3:
        //cout << "swapped down\n" << endl;
        swap(array, location+sqrt(grid_size), location);
        loc[0] = (location+sqrt(grid_size)) / sqrt(grid_size);
        //cout << "Moving to " << location+1 << " from " << location << endl;
        prev = 3;
        break;
    }

   //int** test = convertArraytoMatrix(array, sqrt(grid_size));
   //printMatrix(test, sqrt(grid_size), sqrt(grid_size));
    //cout << endl;
    //prev = location;
    iterate++;
    solve(array, grid_size);
  }
}

int findSmallest(int* array) {
  int smallest = INT_MAX;
  int index = 0;
  for (int i=0; i<4; i++) {
    if (array[i] < smallest) {
      smallest = array[i];
      index = i;
    }
  }
  return index;
}

void swap(int* array, int loc1, int loc2) {
  int temp = array[loc2];
  array[loc2] = array[loc1];
  array[loc1] = temp;
}

bool checkLeft(int position, int n) {
  //cout << "Left check: " << position << ". Previous: " << prev << endl;
  if (position % n == 0 || prev == 1) {
    return false;
  } else {
    return true;
  }
}

bool checkRight(int position, int n) {
  //cout << "Right check: " << position << ". Previous: " << prev << endl;
  if (position % n == (n-1) || prev == 0) {
    return false;
  } else {
    return true;
  }
}

bool checkUp(int position, int n) {
  //cout << "Up check: " << position << ". Previous: " << prev << endl;
  if (position / n == 0 || prev == 3) {
    return false;
  } else {
    return true;
  }
}

bool checkDown(int position, int n) {
  //cout << "Down check: " << position << ". Previous: " << prev << endl;
  if (position / n == (n-1) || prev == 2) {
    return false;
  } else {
    return true;
  }
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
