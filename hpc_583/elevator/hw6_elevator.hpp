// multiple rider case attempt
/*
    Passenger List: elevator maintains list of passengers (pairs of person_id and
    destination floor)
    Batch Processing: elevator picks up as many passengers as possible from the queue
    if they are on the same floor, respecting the maximum occupancy
    Route Management: elevator manages passengers' routes, moving to each passenger's
    destination floor and logging the necessary information
    Complete Logging: detailed logs for each step including entering and exiting the
    elevator are added to provide comprehensive traceability
    also
    Occupancy Declaration: occupancy variable initialized to 0 at the start of the
    elevator function
    Occupancy Management: increment occupancy each time a passenger enters the elevator
    and decrement each time a passenger exits.
*/

#ifndef ELEVATOR_HPP
#define ELEVATOR_HPP
#include <iostream>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <random>
#include <atomic>
#include <vector>
#include <condition_variable>
#include <tuple>
#include <unordered_map>

using namespace std;

const int NUM_FLOORS = 50;
const int NUM_ELEVATORS = 6;
const int MAX_OCCUPANCY = 5;
const int MAX_WAIT_TIME = 5000; // milliseconds

mutex cout_mtx;
mutex queue_mtx;
condition_variable cv;
queue<tuple<int, int, int>> global_queue; // person_id, start_floor, dest_floor
vector<int> elevator_positions(NUM_ELEVATORS, 0);
atomic<int> num_people_serviced(0);
vector<int> global_passengers_serviced(NUM_ELEVATORS, 0);
int npeople;

void elevator(int id) {
    int occupancy = 0; // initialize occupancy for each elevator
    int curr_floor = 0;
    vector<pair<int, int>> passengers; // pairs of (person_id, dest_floor)
    const int time_per_floor = 10;
    // please complete the code segment

    while (true){
        tuple<int, int, int> request;
        bool get_req = false;
        {
            unique_lock<mutex> lock(queue_mtx);
            // wait until elevators get requests or finish all
            cv.wait(lock, [&]() {return !global_queue.empty() || !passengers.empty() || num_people_serviced.load() >= npeople;});
            if (global_queue.empty() && passengers.empty() && num_people_serviced.load() >= npeople) {
                cv.notify_all();
                break; // all done
            }

            // elevator can pick passengers
            if (!global_queue.empty() && occupancy < MAX_OCCUPANCY){
                request = global_queue.front();
                global_queue.pop();
                get_req = true;
            }
        } 
        // elevator receive requests
        if (get_req){
            int person_id, start_floor, dest_floor;
            tie(person_id, start_floor, dest_floor) = request;
            // pickup
            if (curr_floor != start_floor){
                {
                    lock_guard<mutex> lock(cout_mtx);
                    cout << "Elevator(PICKUP) " << id << " moving from FLOOR " << curr_floor << " to FLOOR " << start_floor << endl;
                }
                // moving
                this_thread::sleep_for(chrono::milliseconds(time_per_floor * abs(start_floor - curr_floor)));
                curr_floor = start_floor;
            }
            // passenger enters
            {
                lock_guard<mutex> lock(cout_mtx);
                cout << "Person " << person_id << " entered elevator " << id << endl;
            }
            // occupantship:
            occupancy++;
            passengers.emplace_back(person_id, dest_floor);
            global_passengers_serviced[id]++;
        }

        // elevator service passengers
        if (!passengers.empty()) {
            auto next = passengers.front();
            passengers.erase(passengers.begin());
            int person_id = next.first;
            int dest_floor = next.second;
            // sending passengers
            if (curr_floor != dest_floor) {
                {
                    lock_guard<mutex> lock(cout_mtx);
                    cout << "Elevator(SENDING) " << id << " moving from FLOOR " << curr_floor << " to FLOOR " << dest_floor << endl;
                }
                this_thread::sleep_for(chrono::milliseconds(time_per_floor * abs(dest_floor - curr_floor)));
                curr_floor = dest_floor;
            }
            // person exits
            {
                lock_guard<mutex> lock(cout_mtx);
                cout << "Person " << person_id << " arrived at FLOOR " << dest_floor << endl;
            }
            occupancy--;
            num_people_serviced.fetch_add(1);
        }

    }
    {
        lock_guard<mutex> lock(cout_mtx);
        cout << "Elevator " << id << " has finished servicing all people." << endl;
        cout << "Elevator " << id << " serviced " << global_passengers_serviced[id] << " passengers." << endl;
    }
}

void person(int id) {
    int curr_floor = rand() % NUM_FLOORS;
    int dest_floor = rand() % NUM_FLOORS;
    while (dest_floor == curr_floor) {
        dest_floor = rand() % NUM_FLOORS;
    }
    {
        lock_guard<mutex> lock(cout_mtx);
        cout << "Person " << id << " wants to go from floor " << curr_floor << " to floor " << dest_floor << endl;
    }
    // please complete the code segment
    // request for elevators
    {
        lock_guard<mutex> lock(queue_mtx);
        global_queue.emplace(id, curr_floor, dest_floor);
    }
    cv.notify_one();
}
#endif // ELEVATOR_HPP