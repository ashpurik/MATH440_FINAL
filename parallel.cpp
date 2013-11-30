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

bool dfs(int steps, int currManhat, int* array, int grid_size, int bound, int nextBound, int rank, int iter);
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
void slave(int grid_size, int rank);

void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);

vector<int> searchedElements;
int counter;
int loc[2] = {0, 0};

int nextBound;// = INT_MAX;
int newBound;
int bound;
int found;
int master = 0;
int comp;
int cnt;
int stepSize;

map<char, char> path;
map<unsigned long long, int> visited;
int* numbers;
int iter=0;

int main ( int argc, char *argv[] ) {
  
  int N;
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

  //clock_t start = clock();
  int GRID = N*N;

  if (my_rank == master) {

    cout << "\n---------------------------------" << endl;

    numbers = new int[GRID];
    //numbers = generate_random_array(GRID);
    int stuff[9] = {2, 4, 6, 7, 3, 5, 8, 1, 0};
    //numbers = generate_random_array(GRID);
    for (int i=0; i<GRID; i++) {
      numbers[i] = stuff[i];
    }

    matrix = convertArraytoMatrix(numbers, N);
    printMatrix(matrix, N, N);
    cout << endl;

    bool solvable = isSolvable(numbers, GRID);
    if (!solvable) {
      cout << "Puzzle is not solvable\n" << endl;

      delete numbers;
      freeMatrix(matrix, N);
      solve = 0;
     //MPI::Finalize();
      //return 0;
    } else {
      cout << "Puzzle is solvable!\n" << endl;
      solve = 1;
    }

  }

  MPI_Bcast(&solve, 1, MPI_INT, master, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (solve == 1) {

    if (my_rank == master) {
      int num_steps = ida_star(numbers, GRID);
      cout << "SOLUTION FOUND" << endl;
      cout << "NUM STEPS: " << bound << endl;
    } else {
      slave(GRID, my_rank);
    }

  //if (my_rank == master){
    //printPath(num_steps), printf("\n");
    //solved = convertArraytoMatrix(numbers, N);
    //printMatrix(solved, N, N);
  //}
      MPI::Finalize();

   return 0;
    

  } else {

    if (my_rank == master) {


      cout << "---------------------------------" << endl;

    }

       MPI::Finalize();

   return 0;

  }

  // if (my_rank == master) {
  // delete numbers;
  // freeMatrix(matrix, N);
  // freeMatrix(solved, N);
  // //clock_t end = clock();

  // cout << "---------------------------------" << endl;
  // }

  //cout << setprecision(15) << "EXECUTION TIME: " << double(end-start)/CLOCKS_PER_SEC << " seconds\n" << endl;

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

int findSection(int* array, int grid_size) {
  int i;
  for (i = 0; i < grid_size; i++) {
    if (array[i] == grid_size-1) {
      break;
    }
  }
  return (i /= sqrt(grid_size));
}

int ida_star(int* array, int grid_size) {

  bound = manhattanDistances(array, grid_size);
  cout << "INITIAL BOUND: " << bound << endl;
  int nextBound;
  int* infoToSend = new int[4];
  int oldBound = 0;

  int iter=0;
  bool noSolution = true;

  //while (iter < 5) {
  while (noSolution) {

    //find where the blank spot is
    int section = findSection(array, grid_size);

    nextBound = INT_MAX;
    //path.clear();
    //visited.clear();

    if (bound > oldBound) { 
      oldBound = bound;
      //send starting board & info to appropriate core
      MPI_Send(array, grid_size, MPI_INT, section+1, 123, MPI_COMM_WORLD);
      infoToSend[0] = 0;
      infoToSend[1] = manhattanDistances(array,grid_size);
      infoToSend[2] = bound;
      infoToSend[3] = nextBound;
      cout << "0: I AM SENDING TO " << section+1 << ". ";
      printArray(infoToSend, 4);
      MPI_Send(infoToSend, 4, MPI_INT, section+1, 123, MPI_COMM_WORLD);
    } 
    int finished;
    MPI_Status status;

    MPI_Recv(&finished, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (finished==1) {
      noSolution = false;
    } else {
      bound = finished;
      cout << "0: BOUND: " << bound << endl;
    }

    //iter++;
  }

  return 1;
}

void slave(int grid_size, int rank) {
  int solve;
  int* array = new int[grid_size];
  int* info = new int[4];
  int** arrayInfo = new int*[2];
  stack<int**> unexplored;
  MPI_Status status1, status2;
  MPI_Request first, second;

  bool noSolution = true;

  while(noSolution) {

    MPI_Irecv(array, grid_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &first);
    //cout << rank << ": RECEIVING ARRAY: ";
    //printArray(array, grid_size);
    MPI_Irecv(info, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &second);
    cout << rank << ": RECEIVING: " << info[0] << " " << info[1] << " " << info[2] << " " << info[3] << endl;

    // int tag = status2.MPI_TAG;
    // if (tag == 123) {
    //   bound = info[2];
    // }

    cout << rank << ": SLAVE BOUND: " << bound << endl;

    comp = 0;
    cnt = 0;
    while(!comp && (cnt <= 10)) {
      sleep(1);
      MPI_Test(&first, &comp, &status1);
      MPI_Test(&second, &comp, &status2);
      int tag = status2.MPI_TAG;
      if (tag == 123) {
        bound = info[2];
      }
      cnt++;
    }

    if (comp == 0) {
      MPI_Cancel(&first);
      MPI_Cancel(&second);
      MPI_Wait(&first, &status1);
      MPI_Wait(&second, &status2);
      noSolution = false;
      cout << rank << ": NOTHING RECEIVED" << endl;
    } else {
      cout << rank << ": I RECEIVED AN ARRAY: " << " ";
      printArray(array, grid_size);

      arrayInfo[0] = array;
      arrayInfo[1] = info;
      unexplored.push(arrayInfo);
      //cout << rank << ": STACK SIZE: " << unexplored.size() << endl;

      //cout << "STEPS: " << unexplored.top()[1][0] << endl;
      //cout << "CURRMANHAT: " << unexplored.top()[1][1] << endl;
      //cout << "SLAVE BOUND: " << unexplored.top()[1][2] << endl;
      //cout << "NEXT BOUND: " << unexplored.top()[1][3] << endl;

      ///cout << rank << ": MY STACK SIZE: " << unexplored.size() << " ";
      //printArray(unexplored.top()[0], grid_size);

      while(!unexplored.empty()) {
        int* gridToExplore = unexplored.top()[0];
        if(dfs(unexplored.top()[1][0], unexplored.top()[1][1], gridToExplore, grid_size, unexplored.top()[1][2], unexplored.top()[1][3], rank, 0)) {
          found = 1;
          MPI_Send(&found, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
          noSolution = false;
        } else {
          //cout << rank << ": BOUND: " << bound << ". STEPSIZE: " << stepSize << endl; 
          //cout << rank << ": SLAVE NEWBOUND: " << newBound << endl;
          if (newBound > bound) {
            bound = newBound;
            //MPI_Bcast(&bound, 1, MPI_INT, rank, MPI_COMM_WORLD);
            cout << rank << ": SENDING TO MASTER THE NEW BOUND " << newBound << endl;
            MPI_Send(&newBound, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
            visited.clear();
          }
          //int found = 0;
          //cout << "AFTER DFS NEW BOUND: " << newBound << endl;
          //MPI_Send(&newBound, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
          //visited.clear();
        }
        unexplored.pop();
      }
      //cout << rank << ": AFTER EXPLORING" << endl;
    }

  }

}

bool dfs(int steps, int currManhat, int* array, int grid_size, int bound, int nextBound, int rank, int iter) {

  if (bound == 24) {
    cout << rank << ": ***************************************BOUND IS 24**************************" << endl;
  }

  if (steps + currManhat > bound) {
    nextBound = min(nextBound, steps + currManhat);
    newBound = nextBound;
    //cout << rank << ": DFS NEBOUND: " << newBound << endl;
    //cout << rank << ": IN THE INCREASING BOUND" << endl;
    return false;
  }

  if (isSolution(array, grid_size)) {

    printArray(array, grid_size);
    return true;
  }

  unsigned long long state = 0;
  for (int i = 0; i < grid_size; i++) { // transform the grid into an unsigned long long, to store the state
    state <<= int(sqrt(grid_size)); // move to the left n size bits
    state += array[i]; // add it to the state
  }

  if (visited.count(state) && (visited[state] <= steps)) { // we want to prevent cycling on the grid
    cout << rank << ": I AM CYCLING" << endl;
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
    cout << rank << ": ITERATION: " << iter << ". CURRENT BOUND: " << bound << ". I am exploring to the right: ";
    printArray(array, grid_size);
    path[steps + 1] = 'R';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound, nextBound, rank, iter+1)) { // if ok, no need to restore, just go ahead
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+j+1); // restore
  }
  if (checkUp(i-1)) {
    int new_index = convertCoordinatestoIndex(i-1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j);
    cout << rank << ": ITERATION: " << iter << ". CURRENT BOUND: " << bound << ". I am exploring up ";
    printArray(array, grid_size);
    path[steps + 1] = 'U';
    int section = findSection(array, grid_size);
    MPI_Send(array, grid_size, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    int* infoToSend = new int[4];
    infoToSend[0] = steps+1;
    infoToSend[1] = currManhat+dh;
    infoToSend[2] = bound;
    infoToSend[3] = nextBound;
    cout << rank << ": I AM SENDING TO " << section+1 << ". ";
    printArray(infoToSend, 4);
    MPI_Send(infoToSend, 4, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    swap(array, i*sqrt(grid_size)+j, (i-1)*sqrt(grid_size)+j); // restore
    //return false;
  }
  if (checkLeft(j-1)) {
    int new_index = convertCoordinatestoIndex(i, j-1, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); //swap
    cout << rank << ": ITERATION: " << iter << ". I am exploring to the left ";
    printArray(array, grid_size);
    path[steps + 1] = 'L';
    if (dfs(steps + 1, currManhat + dh, array, grid_size, bound, nextBound, rank, iter+1)) { // if ok, no need to restore, just go ahead
      return true;
    }
    swap(array, i*sqrt(grid_size)+j, i*sqrt(grid_size)+(j-1)); // restore
  }
  if (checkDown(i+1, sqrt(grid_size))) {
    int new_index = convertCoordinatestoIndex(i+1, j, sqrt(grid_size));
    int dh = manhattan(index, new_index, array, sqrt(grid_size));
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j);
    cout << rank << ": ITERATION: " << iter << ". I am exploring down ";
    printArray(array, grid_size);
    path[steps + 1] = 'D';
    int section = findSection(array, grid_size);
    MPI_Send(array, grid_size, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    int* infoToSend = new int[4];
    infoToSend[0] = steps+1;
    infoToSend[1] = currManhat+dh;
    infoToSend[2] = bound;
    infoToSend[3] = nextBound;
    cout << rank << ": I AM SENDING TO " << section+1 << ". ";
    printArray(infoToSend, 4);
    MPI_Send(infoToSend, 4, MPI_INT, section+1, 5, MPI_COMM_WORLD);
    swap(array, i*sqrt(grid_size)+j, (i+1)*sqrt(grid_size)+j); // restore
    //return false;
  }
  //cout << rank << ": FINISHED EXPLORING THIS STATE" << endl;

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

