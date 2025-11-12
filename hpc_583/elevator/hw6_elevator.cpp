#include "hw6_elevator.hpp"

"""
The exercise is to write a threaded elevator scheduler for a 50-story business building with 6 elevators. 
Up to 7500 people may visit the building in a day. 
Your functions should use the C++ threads library to simulate the behavior of the elevators and visitors of the building. 
The driver code takes as input the number of people to service (hint: makes it a bounded buffer). 
After servicing all person requests, have each elevator report the number of passengers it transported, and be certain 
that your code exits cleanly. Please print events to cout in a manner that is thread safe. 

For each person, report their ID and current floor and destination floor (i.e. Person 0 wants to go from floor 7 to floor 9), 
report the elevator number they get on (i.e. Person 1 entered elevator 0), 
report when they arrive at their target floor (i.e. Person 98 arrived at floor 18). 

For each elevator, report the direction of movement from current floor to next floor 
(i.e. Elevator 2 moving from floor 12 to floor 19), join all the elevator threads once all visitors have arrived at their floors, 
and report success (i.e. Job completed!).

1) randomly select the current floor and destination floor for each person
2) I let each elevator start from the ground floor if it matters
3) issue a person thread in the detach state (i.e. thread(person, i).detach(); )
4) each elevator will be its own joinable thread
5) keep track of the current floor and direction of each elevator
6) keep track of the destinations of the visitors who are waiting for elevators on each floor
7) determine which elevator should respond to a given request, based on its current location, 
direction, and the number of passengers it is already carrying
8) the maximum occupancy of the elevators is fixed, say 10
9) move the elevators up or down to their next destination, picking up and dropping off passengers as necessary
10) print events to screen as they occur
11) add some sleep functions to mimic the action of the elevators moving between floors

"""

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " npeople" << endl;
        return 1;
    }

    npeople = atoi(argv[1]);
    if (npeople > 7500) {
        cerr << "Number of people exceeds maximum limit of 7500" << endl;
        return 1;
    }

    thread elevators[NUM_ELEVATORS];
    for (int i = 0; i < NUM_ELEVATORS; ++i) {
        elevators[i] = thread(elevator, i);
    }

    default_random_engine gen;
    uniform_int_distribution<int> dist(0, MAX_WAIT_TIME);

    for (int i = 0; i < npeople; ++i) {
        int wait_time = dist(gen);
        this_thread::sleep_for(chrono::milliseconds(wait_time));
        thread(person, i).detach();
    }
    for (auto &e : elevators) {
        e.join();
    }

    cout << "Job completed!" << endl;

    int total_passengers_serviced = 0;
    for (int i = 0; i < NUM_ELEVATORS; ++i) {
        total_passengers_serviced += global_passengers_serviced[i];
    }
    cout << "Total passengers serviced by all elevators: " << total_passengers_serviced << endl << flush;
    
    return 0;
}
