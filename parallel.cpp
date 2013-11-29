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

# include "mpi.h"
# include <stack>
# include <vector>
# include <queue>

#include <iterator>

using namespace std;

int* generate_random_array(int grid_size);
int** convertArraytoMatrix(int* array, int n);
int convertCoordinatestoIndex(int i, int j, int n);

bool isSolution(int* array, int grid_size);
bool isSolvable(int* array, int grid_size);
void findCycle(std::map<int,int> theMap, int initial_key, int curr_key, int n);
int manhattan(int start, int end, int* array, int n);

int master_tasks(int* array, int grid_size, int num_cores);
void slave_tasks(int my_rank);

bool checkMove(int i, int j, int n);
bool checkLeft(int j);
bool checkUp(int i);
bool checkDown(int i, int n);
bool checkRight(int j, int n);
void swap(int* array, int loc1, int loc2);


void freeMatrix(int** matrix, int row);
void printMatrix(int** matrix, int row, int col);
void printArray(int* array, int size);
//void printPath(int previous);

int N = 3;
int GRID = (N*N);

vector<int> searchedElements;
int counter, last(GRID-1);
int loc[2] = {0, 0};

int nextBound;
int master = 0;
int can_be_solved = 0;
bool solve = false;

int* numbers;
int** matrix;
int** solved;

deque<unsigned long long> already_visited;
queue<int*> master_distributes;
//vector<int*> master_distributes;

//map<int, int> path;
map<unsigned long long, int> visited;

int main ( int argc, char *argv[] ) {

  int my_rank;
  int num_cores;

  MPI::Init(argc, argv);

  my_rank = MPI::COMM_WORLD.Get_rank();
  num_cores = MPI::COMM_WORLD.Get_size();

  if (my_rank == master){

    //master_distributes.push_back(new int[N*N]);
    numbers = new int[GRID];
    numbers = generate_random_array(GRID);

    matrix = convertArraytoMatrix(numbers, N);
    printMatrix(matrix, N, N);
    cout << endl;

    bool solvable = isSolvable(numbers, GRID);
    if (solvable){
  cout << "Puzzle is solvable!" << endl;
  can_be_solved = 1;

    }
    else {
        cout << "NOT A SOLVABLE MATRIX" << endl;

  can_be_solved = 0;

  delete numbers;
        freeMatrix(matrix, N);
    }
  }  

  MPI::COMM_WORLD.Bcast(&can_be_solved, 1, MPI_INT, master);
  MPI::COMM_WORLD.Barrier();

  if (can_be_solved == 1){
    if (my_rank == master){
       unsigned long long state = 0;
       for (int i = 0; i < GRID; i++) { // transform 16 numbers into 6ROW_SIZE bits, exactly into ULL
           state <<= int(sqrt(GRID)); // move left ROW_SIZE bits
           state += numbers[i]; // add this digit (max 15 or 1111)
        }

       already_visited.push_back(state);
       //master_distributes.push_back(numbers); 
       master_distributes.push(numbers);
       int num_steps = master_tasks(numbers, GRID, num_cores);
    }
    else {
       slave_tasks(my_rank);
    }
  }

  else{
    MPI::Finalize();
    return 0;
  }

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
    int initial_i = start / n;   //current i location of space
    int initial_j = start % n;   //current j location of space
    int possible_i = end / n;    //the position of i possible swap
    int possible_j = end % n;    //the position of j possible swap
    int dist_blank_from_goal = (abs(initial_i - goal_i) + abs(initial_j - goal_j));
    int dist_possible_from_goal = (abs(possible_i - goal_i) + abs(possible_j - goal_j));
    return dist_blank_from_goal - dist_possible_from_goal;
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


int master_tasks(int* array, int grid_size, int num_cores) {

  //while (solve == false) {

    int i = 1;
    MPI::Status my_status;
    int *send_state = new int[N*N];
    int iter=0;
    //distribute states to be checked
    //while(iter < 2) {
    while(!master_distributes.empty()){

     // erase first element
     //master_distributes.erase( master_distributes.begin() );
     //cout << "\n\nVector integers after erasing first element: ";
     //std::copy( master_distributes.begin(), master_distributes.end(), output );

      //for(int i=0; i<master_distributes.size(); i++) {
      //  printArray(master_distributes.at(i), N*N);
      //}

      send_state = master_distributes.front();

      master_distributes.pop();

      // master_distributes.erase(master_distributes.begin());
      MPI_Send(send_state, N*N, MPI_INT, i%(num_cores - 1) + 1, 123, MPI_COMM_WORLD);  //i%(num_cores-1) + 1 set to 3
      //recv swapped states
      int *received_state = new int[N*N*4];
      MPI::COMM_WORLD.Recv(received_state, 4*N*N, MPI_INT, i%(num_cores - 1) + 1, 123, my_status);//i%(num_cores-1) + 1 set to 3
      int *individual_state = new int[N*N];

      for (int m = 0; m < 4; m++) {

        for (int i = 0; i < N*N; i++){
          individual_state[i] = received_state[m*N*N +i];
        }

        cout << "MASTER: ";
        printArray(individual_state, N*N);

        //change the received state into a unsigned value to be added to the list of visited 
        unsigned long long state = 0;
        for (int j = 0; j < N*N; j++) {          // transform 16 numbers into 6ROW_SIZE bits, exactly into ULL
          state <<= int(sqrt(N*N));         // move left ROW_SIZE bits
          state += individual_state[j]; // add this digit (max 15 or 1111)
        }

        for (int j = 0; j < N*N; j++){
        //if new, push onto stack of visited and to distribute....... i have a funny feeling this isn't 100% correct
          if (find(already_visited.begin(), already_visited.end(), state) == already_visited.end()){
            already_visited.push_back(state);
            master_distributes.push(individual_state);
            cout << "IN HERE" << endl;
          }
        }
        cout << already_visited.size() << endl;
      }

    i++;
    }
    solve = true;

  //}

  return 1;
}

void slave_tasks(int my_rank){

  MPI_Status my_status;
  
  while (solve == false){
    int *tobechecked = new int[N*N];
    
    MPI_Recv(tobechecked, N*N, MPI_INT, master, 123, MPI_COMM_WORLD, &my_status);

    cout << "MY RANK : " << my_rank << " WHAT I'M SOLVING : ";
    printArray(tobechecked, N*N);

    //check if the received array is the solution 
    bool isSOLUTION = isSolution(tobechecked, N*N);

    //if the solution has been found, do shit
    if (isSOLUTION == true) {
       solve = true; 
       MPI::COMM_WORLD.Bcast(&solve, 1, MPI::BOOL, my_rank);
       break;
    }

    int i, new_i, new_j;
    for (i = 0; i < N*N; i++) {
       if (tobechecked[i] == N*N-1) {
          break;
       }
    }

    int j = i % N;
    i /= N;

    int index = convertCoordinatestoIndex(i, j, N*N);

    int *shiftedgrids = new int[4*N*N];

    if (checkRight(j+1, N)) {
      swap(tobechecked, i*N+j, i*N+j+1);
      for (int p = 0; p < N*N; p++){
        shiftedgrids[p] = tobechecked[p];
      }
      swap(tobechecked, i*N+j, i*N+j+1); // restore
    }

    if (!checkRight(j+1, N)){
      for (int p = 0; p < N*N; p++){
        shiftedgrids[p] = tobechecked[p];
      }
    }

    if (checkUp(i-1)) {
      swap(tobechecked, i*N+j, (i-1)*N+j);
      for (int p = 0; p < N*N; p++){
        shiftedgrids[N*N + p] = tobechecked[p];
      }
      swap(tobechecked, i*N+j, (i-1)*N+j); // restore
    }

    if (!checkUp(i-1)){
      for (int p = 0; p < N*N; p++){
        shiftedgrids[N*N + p] = tobechecked[p];
      }
    }

    if (checkLeft(j-1)) {
      swap(tobechecked, i*N+j, i*N+(j-1)); //swap
      for (int p = 0; p < N*N; p++){
        shiftedgrids[2*N*N + p] = tobechecked[p];
      }
      swap(tobechecked, i*N+j, i*N+(j-1)); // restore
    }

    if (!checkLeft(j-1)){
      for (int p = 0; p < N*N; p++){
        shiftedgrids[2*N*N + p] = tobechecked[p];
      }
    }

    if (checkDown(i+1, N)) {
      swap(tobechecked, i*N+j, (i+1)*N+j);
      for (int p = 0; p < N*N; p++){
        shiftedgrids[3*N*N + p] = tobechecked[p];
      }
      swap(tobechecked, i*N+j, (i+1)*N+j); // restore
    }

    if (!checkDown(i+1, N)){
      for (int p = 0; p < N*N; p++){
        shiftedgrids[3*N*N + p] = tobechecked[p];
      }
    }

    printArray(shiftedgrids, 4*N*N);
    MPI::COMM_WORLD.Send(shiftedgrids, 4*N*N, MPI_INT, master, 123);
  }
}

//void printPath(int previous) {
//  if (previous == 0) {
//    return;
//  }
//  printPath(previous - 1);
//  printf("%c", path[previous]);
//}

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