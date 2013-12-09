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
# include <stack>

//MPI
# include "mpi.h"

using namespace std;

int* generate_random_array(int grid_size);
int** convertArraytoMatrix(int* array, int n);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattanDistances(int* array, int grid_size);
void swap(int* array, int loc1, int loc2);
bool dfs(int steps, int currManhat, int* array, int grid_size, int nextBound, int bound, int rank, int iter);
int manhattan(int start, int end, int* array, int n);
int convertCoordinatestoIndex(int i, int j, int n);
bool checkLeft(int j);
bool checkUp(int i);
bool checkDown(int i, int n);
bool checkRight(int j, int n);
int findSection(int* array, int grid_size);
bool isSolution(int* array, int grid_size);
void slave(int grid_size, int rank);

void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int counter;
int loc[2] = {0, 0};

int nextBound = INT_MAX;
int newBound;
int bound;
int found;
int master = 0;
int comp;
int cnt;
int sleeptime; 
int cntmax; 

map<char, char> path;
map<unsigned long long, int> visited;
int* numbers;
int iter=0;

int main ( int argc, char *argv[] ) {
  
  int N;
  double start, end;
  int my_rank, num_cores;
  int** matrix;
  int** solved;
  int solve = 0;

  //grabs value of N from command line
  if (argc < 2) {
    cout << "Please specify N" << endl;
  } else {
   N = atoi( argv[1]);
  }

  MPI::Init (argc, argv);

  // Get the number of processes
  num_cores = MPI::COMM_WORLD.Get_size();

  //Get the individual process ID
  my_rank = MPI::COMM_WORLD.Get_rank();

  int GRID = N*N;

  if (my_rank == 0) { 
    start = MPI::Wtime();
  }

  //master generates a random array
  if (my_rank == master) {

    cout << "\n---------------------------------" << endl;

    numbers = new int[GRID];
    numbers = generate_random_array(GRID);
    matrix = convertArraytoMatrix(numbers, N);
    printMatrix(matrix, N, N);
    cout << endl;

    //master checks if the random arrangement produces a solvable puzzle
    bool solvable = isSolvable(numbers, GRID);
    if (!solvable) {
      cout << "Puzzle is not solvable\n" << endl;
      delete numbers;
      freeMatrix(matrix, N);
      solve = 0;
    } else {
      cout << "Puzzle is solvable!\n" << endl;
      solve = 1;

      if (N == 3) { 
	bound = 30;
        sleeptime = 500;
        cntmax = 80;
      } else if (N == 4) {
        bound = 80;
        sleeptime = 5000;
        cntmax = 10000;
      } else if (N == 5) {
        bound = 200;
      }

    }

  }

  //master notifies all cores whether puzzle is solvable or not
  //and specifies their timers
  MPI_Bcast(&solve, 1, MPI_INT, master, MPI_COMM_WORLD);
  MPI_Bcast(&sleeptime, 1, MPI_INT, master, MPI_COMM_WORLD);
  MPI_Bcast(&cntmax, 1, MPI_INT, master, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (solve == 1) {

    if (my_rank == master) {

      int* infoToSend = new int[3];
      bool noSolution = true;

      //find where the blank spot is
      int section = findSection(numbers, GRID);

      //send the appropriate core the grid & starting information to begin solving
      MPI_Send(numbers, GRID, MPI_INT, section+1, 123, MPI_COMM_WORLD);
      infoToSend[0] = 0;
      infoToSend[1] = manhattanDistances(numbers, GRID);
      infoToSend[2] = bound;
      MPI_Send(infoToSend, 3, MPI_INT, section+1, 123, MPI_COMM_WORLD);

      //keep waiting for a core to send the solution (1) 
      while (noSolution) {

        int finished;
        MPI_Status status;

        MPI_Recv(&finished, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (finished==1) {
          noSolution = false;

          end = MPI::Wtime();
          cout << "\nElapsed wall clock time: " << end-start << " seconds. " << endl;

 	  //once a solution has been found, call MPI_Abort on all cores so they stop searching for a solution
	  delete numbers;
          freeMatrix(matrix, N);
          int errorcode = 0;
          MPI_Abort(MPI_COMM_WORLD, errorcode);
        } 
      }

    } else {
      slave(GRID, my_rank);
    }
    
  } else {

    if (my_rank == master) {

      cout << "---------------------------------" << endl;
      end = MPI::Wtime();
      cout << "Elapsed wall clock time: " << end-start << " seconds. " << endl;

    }

    MPI::Finalize();

  }
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

//locates the row of the board where the blank spot is located
int findSection(int* array, int grid_size) {
  int i;
  for (i = 0; i < grid_size; i++) {
    if (array[i] == grid_size-1) {
      break;
    }
  }
  return (i /= sqrt(grid_size));
}

//Each core receives a grid to explore and performs DFS on it
void slave(int grid_size, int rank) {
  int* array = new int[grid_size];
  int* info = new int[3];
  MPI_Status status1, status2;
  MPI_Request first, second;

  bool noSolution = true;

  //keep receiving arrays to explore until a solution has been found
  while(noSolution) {

    MPI_Irecv(array, grid_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &first);
    MPI_Irecv(info, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &second);

    comp = 0;
    cnt = 0;
    //we need to wait some time before canceling MPI_Irecv
    while(!comp && (cnt <= cntmax)) {
      usleep(sleeptime);
      MPI_Test(&first, &comp, &status1);
      MPI_Test(&second, &comp, &status2);
      cnt++;
    }

    if (comp == 0) {
      MPI_Cancel(&first);
      MPI_Cancel(&second);
      MPI_Wait(&first, &status1);
      MPI_Wait(&second, &status2);
      noSolution = false;
      cout << rank << ": NOTHING RECEIVED. " << endl;

    //explore the grid!
    } else {
      if(dfs(info[0], info[1], array, grid_size, nextBound, info[2], rank, 0)) {
        found = 1;
        MPI_Send(&found, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
        int** solved = convertArraytoMatrix(array, sqrt(grid_size));
        printMatrix(solved, sqrt(grid_size), sqrt(grid_size));
        freeMatrix(solved, sqrt(grid_size));
        noSolution = false;
      } else {
        //tell the master that the bound was increased
        if (newBound > bound) {
          found = 0;
          bound = newBound;
          MPI_Send(&found, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
        }
      }
    }

  }
}

//dfs function remains mostly the same as in the serial code
bool dfs(int steps, int currManhat, int* array, int grid_size, int nextBound, int bound, int rank, int iter) {

  if (steps + currManhat > bound) {
    nextBound = min(nextBound, steps + currManhat);
    newBound = nextBound;
    return false;
  }

  if (isSolution(array, grid_size)) {
    cout << "CORE " << rank << " FOUND THE SOLUTION\n";
    cout << "***************************************************\n";
    cout << "***************************************************\n";
    return true;
  }

  unsigned long long state = 0;
  for (int i = 0; i < grid_size; i++) { // transform the grid into an unsigned long long, to store the state
    state <<= int(sqrt(grid_size)); // move to the left n size bits
    state += array[i]; // add it to the state
  }

  if (visited.count(state) && (visited[state] <= steps)) { // we want to prevent cycling on the grid
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

  //keep recursively checking to the right for the given core
  if(checkRight(j+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i, j+1, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1);
    if (dfs(steps + 1, currManhat + dh, array, grid_size, nextBound, bound, rank, iter+1)) { 
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1); // restore
  }
  //send the new grid to the appropriate core for exploring
  if (checkUp(i-1)) {
    int new_index = convertCoordinatestoIndex(i-1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j);
    int section = findSection(array, grid_size);
    MPI_Send(array, grid_size, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    int* infoToSend = new int[3];
    infoToSend[0] = steps+1;
    infoToSend[1] = currManhat+dh;
    infoToSend[2] = bound;
    MPI_Send(infoToSend, 3, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j); // restore
  }
  //keep recursively checking to the left for the given core
  if (checkLeft(j-1)) {
    int new_index = convertCoordinatestoIndex(i, j-1, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); //swap
    if (dfs(steps + 1, currManhat + dh, array, grid_size,  nextBound, bound, rank, iter+1)) { 
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); // restore
  }
  //send the new grid to the appropriate core for exploring
  if (checkDown(i+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i+1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j);
    int section = findSection(array, grid_size);
    MPI_Send(array, grid_size, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    int* infoToSend = new int[3];
    infoToSend[0] = steps+1;
    infoToSend[1] = currManhat+dh;
    infoToSend[2] = bound;
    MPI_Send(infoToSend, 3, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j); // restore
  }

  return false;
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

