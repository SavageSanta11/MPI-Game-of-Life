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
        // MPI_Request send1, send2, send3, send4, recv1, recv2, recv3, recv4;
        // MPI_Status statuses[8];
        if(mpirank != 0 && mpirank != mpisize -1){
            MPI_Request send1, recv1, send2, recv2;
            //send top row above
        //send top row of the current rank to the rank just above
            MPI_Isend(&currGrid[1][0], local_cols_ghost, MPI::INT, upper_neighbour, 0,MPI_COMM_WORLD, &send1); //start from 0th clum, &take the address

            //receive bottom row from below
            MPI_Irecv(&currGrid[local_rows+1][0], local_cols_ghost, MPI::INT, lower_neighbour, 0,MPI_COMM_WORLD, &recv1);

            //send bottom row below
            MPI_Isend(&currGrid[local_rows][0], local_cols_ghost, MPI::INT, lower_neighbour, 0,MPI_COMM_WORLD, &send2);

            //receive top row from above
            MPI_Irecv(&currGrid[0][0], local_cols_ghost, MPI::INT, upper_neighbour, 0,MPI_COMM_WORLD, &recv2);

            MPI_Request requests[4] = { send1, recv1, send2, recv2 };
            MPI_Status statuses[4];
            MPI_Waitall(4, requests, statuses);
        }
        else if (mpirank == 0) {
            
            MPI_Request send3, recv3;

            //send bottom row below
            MPI_Isend(&currGrid[local_rows][0], local_cols_ghost, MPI::INT, lower_neighbour, 0, MPI_COMM_WORLD,&send3);
            //receive bottom row from below
            MPI_Irecv(&currGrid[local_rows+1][0], local_cols_ghost, MPI::INT, lower_neighbour, 0,MPI_COMM_WORLD, &recv3);

            MPI_Request requests[2] = { send3, recv3};
            MPI_Status statuses[2];
            MPI_Waitall(2, requests, statuses);
        }
        else if (mpirank == mpisize-1) {
            
            MPI_Request send4, recv4;

            //send top row above
            //send top row of the current rank to the rank just above
            MPI_Isend(&currGrid[1][0], local_cols_ghost, MPI::INT, upper_neighbour, 0, MPI_COMM_WORLD, &send4); //start from 0th clum, &take the address
            // receive top row from rank-1
            //receive top row from above
            MPI_Irecv(&currGrid[0][0], local_cols_ghost, MPI::INT, upper_neighbour, 0, MPI_COMM_WORLD, &recv4);

            MPI_Request requests[2] = { send4, recv4};
            MPI_Status statuses[2];
            MPI_Waitall(2, requests, statuses);
        }


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
            MPI_Request requests[mpisize];
            MPI_Status statuses[mpisize];
        if (mpirank != mpiroot){
            for (int irow = 1; irow <= local_rows; ++irow){
                MPI_Isend(&currGrid[irow][1], cols, MPI::INT, mpiroot, 0, MPI_COMM_WORLD, &requests[mpirank-1]);
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
                MPI_Irecv(&buff[0], cols, MPI::INT, sourcerank, 0, MPI_COMM_WORLD, &requests[sourcerank-1]);
                MPI_Wait(&requests[sourcerank - 1], &statuses[sourcerank - 1]);
                for (int jcol = 0; jcol<cols; ++jcol){
                    
                    if (buff[jcol] == 1){
                        //cout << (irecv + (local_rows*sourcerank)) << "," << jcol <<'\n';
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