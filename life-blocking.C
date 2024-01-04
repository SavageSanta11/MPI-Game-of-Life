#include <stdio.h>
#include <mpi.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include<bits/stdc++.h>

using namespace std;

void read_input_file(vector<vector<int> > &life, string const &input_file_name, int rank, int local_rows) {

    // Open the input file for reading.
    ifstream input_file;
    input_file.open(input_file_name);
    if (!input_file.is_open())
        perror("Input file cannot be opened");

    string line, val;
    int x, y;
     // Populate the life matrix

        auto lower_limit = rank*local_rows;
        auto upper_limit = local_rows*(rank+1)-1;


    while (getline(input_file, line)) {
        stringstream ss(line);

        // Read x coordinate.
        getline(ss, val, ',');
        x = stoi(val);

        // Read y coordinate.
        getline(ss, val);
        y = stoi(val);

       
        if (x>=lower_limit && x<= upper_limit){
            
                life[x-(local_rows*rank)][y] = 1;
            
        }

    }
    input_file.close();
}



int main(int argc, char **argv)
   {

      MPI::Init(argc, argv);

      auto mpisize = MPI::COMM_WORLD.Get_size();
      auto mpirank = MPI::COMM_WORLD.Get_rank();

      auto mpiroot = 0;

      int rows;
      int cols;
      int time;
      string input_file_name;
      
      
      if(mpirank == mpiroot){
        //root will get the data
        if (argc != 5){
            cout<<"Incorrect number of command line arguments, 4 expected"<<endl;
            exit(1);
        }
        input_file_name = "data/" + string(argv[1]);;
        time = atoi(argv[2]);
        rows = atoi(argv[3]);
        cols = atoi(argv[4]);
        
        
        if(rows < 1 || cols< 1 || time<1){
            cerr<<"Incvalid command line arguments"<<endl;
        }
    }
    int l = input_file_name.length();
    MPI::COMM_WORLD.Bcast(&rows, 1, MPI::INT, mpiroot);
    MPI::COMM_WORLD.Bcast(&cols, 1, MPI::INT, mpiroot);
    MPI::COMM_WORLD.Bcast(&time, 1, MPI::INT, mpiroot);
    MPI::COMM_WORLD.Bcast(&l, 1, MPI::INT, mpiroot);
    MPI::COMM_WORLD.Bcast(&input_file_name[0], l+1, MPI::CHAR, mpiroot);

    auto local_rows = rows/mpisize;
    if(mpirank == mpisize - 1){
        local_rows += rows %mpisize; //put the remaining rows on the last rank
    }

    auto local_rows_ghost = local_rows + 2;
    auto local_cols_ghost = cols + 2;


    vector<vector<int> > currGrid(local_rows, vector<int>(cols, 0));
    vector<vector<int> > nextGrid(local_rows_ghost, vector<int>(local_cols_ghost, 0));
    read_input_file(currGrid, input_file_name, mpirank, local_rows);

    // padding the currGrid

    for (int g = 0; g < currGrid.size(); ++g){
        currGrid[g].insert(currGrid[g].begin(),1, 0);
        currGrid[g].push_back(0);
    }
    currGrid.insert(currGrid.begin(),1, vector<int>(cols+2, 0));
    currGrid.push_back(vector<int>(cols+2, 0));


    //there shouldn't be a wrap around
    auto upper_neighbour = (mpirank == 0) ? mpisize - 1: mpirank - 1;
    auto lower_neighbour = (mpirank == mpisize - 1) ? 0: mpirank + 1;

    const int ALIVE = 1;
    const int DEAD = 0;

    // time loop
    double t1 = MPI_Wtime();
    for (auto itime = 0; itime < time; ++itime){
        if(mpirank != 0 && mpirank != mpisize -1){
            //send top row above
            //send top row of the current rank to the rank just above
            MPI::COMM_WORLD.Send(&currGrid[1][0], local_cols_ghost, MPI::INT, upper_neighbour, 0); //start from 0th clum, &take the address

            //receive bottom row from below
            MPI::COMM_WORLD.Recv(&currGrid[local_rows+1][0], local_cols_ghost, MPI::INT, lower_neighbour, 0);

            //send bottom row below
            MPI::COMM_WORLD.Send(&currGrid[local_rows][0], local_cols_ghost, MPI::INT, lower_neighbour, 0);

            //receive top row from above
            MPI::COMM_WORLD.Recv(&currGrid[0][0], local_cols_ghost, MPI::INT, upper_neighbour, 0);
        }
        else if (mpirank == 0) {
            //send bottom row below
            MPI::COMM_WORLD.Send(&currGrid[local_rows][0], local_cols_ghost, MPI::INT, lower_neighbour, 0);
            //receive bottom row from below
            MPI::COMM_WORLD.Recv(&currGrid[local_rows+1][0], local_cols_ghost, MPI::INT, lower_neighbour, 0);
        }
        else if (mpirank == mpisize-1) {
            //send top row above
            //send top row of the current rank to the rank just above
            MPI::COMM_WORLD.Send(&currGrid[1][0], local_cols_ghost, MPI::INT, upper_neighbour, 0); //start from 0th clum, &take the address
            // receive top row from rank-1
            //receive top row from above
            MPI::COMM_WORLD.Recv(&currGrid[0][0], local_cols_ghost, MPI::INT, upper_neighbour, 0);
        }

        //update the grid
        for (auto irow = 1; irow <= local_rows; ++irow){
            for (auto icol = 1; icol <= cols; ++icol){
                auto alive_neighbours = 0;

                for (auto jrow = irow - 1; jrow <= irow+1; ++jrow){
                    for (auto jcol = icol-1; jcol<= icol+1; ++jcol){
                        if ((jrow != irow || jcol != icol) &&(currGrid[jrow][jcol] == ALIVE)){
                            alive_neighbours += 1;
                        }
                    }
                }
                if (alive_neighbours < 2){
                    nextGrid[irow][icol] = DEAD;
                }
                if (currGrid[irow][icol] == ALIVE && (alive_neighbours == 2 || alive_neighbours == 3)){
                    nextGrid[irow][icol] = ALIVE;
                }
                if (alive_neighbours > 3){
                    nextGrid[irow][icol] = DEAD;
                }
                if (currGrid[irow][icol] == DEAD && (alive_neighbours == 3)){
                    nextGrid[irow][icol] = ALIVE;
                }
            }
        }
        for (auto irow = 1; irow<=local_rows; ++irow){
            for (auto icol = 1; icol <= cols; ++icol){
                currGrid[irow][icol] = nextGrid[irow][icol];
            }
        }

        if (itime == time - 1){
            //display current grid on screen
            vector< pair <int,int> > storage;
        if (mpirank != mpiroot){
            for (int irow = 1; irow <= local_rows; ++irow){
                MPI::COMM_WORLD.Send(&currGrid[irow][1], cols, MPI::INT, mpiroot, 0);
            }
        }

        if (mpirank == mpiroot){
            for (auto irow = 1; irow<=local_rows; ++irow){
                for (auto icol = 1; icol <= cols; ++icol){
               
                    if (currGrid[irow][icol] == 1){
                       
                        storage.push_back(make_pair((irow-1), (icol-1)));
                    }
            }

        }

        for (auto sourcerank = 1; sourcerank < mpisize; ++sourcerank){
            auto nrecv = rows/mpisize;
            if (sourcerank == mpisize - 1){
                nrecv += rows % mpisize;
            }
            vector<int> buff(cols,0);
            for (auto irecv = 0; irecv < nrecv; ++irecv){
                MPI::COMM_WORLD.Recv(&buff[0], cols, MPI::INT, sourcerank, 0);
                for (int jcol = 0; jcol<cols; ++jcol){
                    
                    if (buff[jcol] == 1){
                       
                        storage.push_back(make_pair((irecv + (local_rows*sourcerank)), jcol));
                    }
                }
                
            }

        }
        }
            if (mpirank == mpiroot){
                sort(storage.begin(),storage.end());
        

            ofstream output_file;
            string output_file_name = input_file_name.substr(0, input_file_name.length() - 5);
            output_file.open(output_file_name + "." + to_string(time) +
                    ".csv");
            if (!output_file.is_open())
                perror("Output file cannot be opened");

            // Output each live cell on a new line.
            for (int i=0; i<storage.size(); i++) {
                    output_file << storage[i].first << "," << storage[i].second << "\n";
                    
            }
            output_file.close();
            }
        }
        
        
    }
      double t2 = MPI_Wtime();
      
      float global_time_taken;
        float max_time_taken;
        float local_time_taken = t2-t1;
        float min_time_taken;
        MPI_Reduce(&local_time_taken, &min_time_taken, 1, MPI_FLOAT, MPI_MIN, 0,MPI::COMM_WORLD);
        MPI_Reduce(&local_time_taken, &global_time_taken, 1, MPI_FLOAT, MPI_SUM, 0, MPI::COMM_WORLD);
        MPI_Reduce(&local_time_taken, &max_time_taken, 1, MPI_FLOAT, MPI_MAX, 0,MPI::COMM_WORLD);
        
        // Write out the final state to the output file.
        if (mpirank ==mpiroot){
            cout << "TIME: Min: " << min_time_taken << " s Avg: " <<  (global_time_taken)/mpisize << " s Max: " << max_time_taken << " s"<<'\n';
        }
      MPI_Finalize();
     
   
      return 0;

   }
