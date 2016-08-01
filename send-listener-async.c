// file: send-listener-async.c
//
// This program demonstrates how to do send and async receive on
//  multiple threads

// equation for updated location: x(t) = v*t + (1/2)a*t^2 // https://www.grc.nasa.gov/www/k-12/airplane/motion.html


#include <stdio.h>
#include <inttypes.h>
#include <sys/select.h>
#include <chrono>
#include <thread>
#include <lcm/lcm.h>
#include "exlcm_location_t.h"
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>
#include <cmath>

// for openMP/multithreading support
#include <omp.h>

#include <iostream>

using namespace std;

bool stopRecv = false; // indicate to stop the message receive event loop
bool ignoreSelfMessages = true; // ignore messages from the same pid
const int DIM = 3;
const int RADIUS = 1; // the radius around a position is used to determine potential collisions

// -----

struct PLoc {
	int pid; // process id
	
	double myLocation[3];
	double myVelocity[3];
	int logicalTimeStamp;

	// store each process' location data:
	//   key   .. process id
	//   value .. location data
	unordered_map<int, const _exlcm_location_t*> locationMap;
	PLoc(int p) : pid(p), logicalTimeStamp(0) {};
};

// -----


long getSystemTimeMS() {
	auto now = chrono::system_clock::now();
    auto now_ms = chrono::time_point_cast<chrono::milliseconds>(now);
    auto value = now_ms.time_since_epoch();
    long duration = value.count();
    chrono::milliseconds dur(duration);
    return duration;
}

long getSystemTimeNS() {
	auto now = chrono::system_clock::now();
    auto now_ns = chrono::time_point_cast<chrono::nanoseconds>(now);
    auto value = now_ns.time_since_epoch();
    long duration = value.count();
    chrono::nanoseconds dur(duration);
    return duration;
}

void print3DArray(
	const double array[3]
) {
	printf("[%f, %f, %f]", array[0], array[1], array[2]);
}

double distance3D(
	const double* A,
	const double* B
) {
	return sqrt( pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2) );
}

void copy3DArray(
	double* src,
	double* dest
) {
	for(int i = 0; i < DIM; ++i)
		dest[i] = src[i];
}

int randomInt(int neg, int max) {
	return ( rand() % (-1*neg+max) + neg + 1 );
}

int PlusOrMinusOne() {
    // return ((getSystemTimeMS() % 2) * 2 - 1 );
    return ((rand() % 2) * 2 - 1 );
}


void printMyLocation(
	PLoc* pLoc
) {
	cout << endl; cout << " ****** ";
    cout << " My location at log.time= " << pLoc->logicalTimeStamp << " (pid=" << pLoc->pid << ") = "; print3DArray( pLoc->myLocation );
    cout << " ****** " << endl; cout << endl;
}

void printLocation(
	const exlcm_location_t * msg,
	const char * channel,
	PLoc* pLoc
) {
    if(!ignoreSelfMessages || pLoc->pid != (int)msg->pid) {
    	printf("Received message on channel \"%s\":\n", channel);
    	printf("  pid              = %"PRId64"\n", msg->pid);
    	printf("  timestampLogical = %"PRId64"\n", msg->timestampLogical);
    	printf("  timestampMS      = %"PRId64"\n", msg->timestampMS);
    	printf("  position         = "); print3DArray( msg->position ); printf("\n");
    	printf("  velocity         = "); print3DArray( msg->velocity ); printf("\n");
    	// printf("  destination      = "); print3DArray( msg->destination ); printf("\n");
    }
}


// returns true if the distance between two center points is less than twice the RADIUS
//  (immediate since a collision is immediately to expect if no stop)
bool immediateCollisionCheck(
	const double* A,
	const double* B,
	int p1,
	int p2,
	int p2LogTime
) {
	double distance = distance3D( A, B );
	if( distance < 2*RADIUS ) {
		cout << " WARN2. An immediate collision between process: " << p1 << " and process: " << p2 << " is about to happen! distance = " << distance << " / log.time of p2 = " << p2LogTime << endl;
		cout << "   Position A: "; print3DArray( A );
		cout << "   Position B: "; print3DArray( B );
		cout << "   distance: " << distance; cout << endl;
		return true;
	}
	return false;
}


void futureCollisionCheck(
	const double* A,
	double* vA, // might be changed if future collision
	const double* B,
	const double* vB,
	int p1,
	int p2,
	int p2LogTime
) {

	// collision check. example:
	// for process A: A0 = (12,6), vA = (-1,0).
	// for process B: B0 = (5,2), vB = (1,1).
	// ==> then a collision would happen in time step = 3. since A at (9,6) and B at (8,5)
	double diffA0B0[3] = { A[0]-B[0], A[1]-B[1], 0 };
	double diffvBvA[3] = { vB[0]-vA[0], vB[1]-vA[1], 0 };

	int factorX = (int)(diffA0B0[0] / diffvBvA[0]);
	if( factorX >= 1 && factorX < 100 && abs( factorX * diffvBvA[1] - diffA0B0[1] ) <= 1 ) {
		cout << " WARN3. A future collision between process: " << p1 << " and process: " << p2 << " could happen at log.timeStamp = " << factorX << " / current log.time of p2 = " << p2LogTime << endl;
		cout << "   Position A: "; print3DArray( A );
		cout << "   Position B: "; print3DArray( B );
		cout << "   diffA0B0: "; print3DArray( diffA0B0 );
		cout << "   diffvBvA: "; print3DArray( diffvBvA );
		
		// only the lower process id is allowed to change its velocity vector
		//  (avoid both processes change the vel. vector and create another collision scenario)
		if( p1 < p2 ) {
			// change the velocity to course correct
			// TODO: here just a simple (non smart) change of the vel. vector. in the future:
			//  - take into account the destination
			//  - take into account other processes as well. (since we do not want to crash into those as well)
			vA[0] = ((int)vA[0] + 1) % 2; // velocity vector passed by reference. so changes are seen by send thread
			vA[1] = ((int)vA[1] + 1) % 2;
		}
		
	}

}


static void my_handler(const lcm_recv_buf_t *rbuf,
                       const char * channel, 
                       const exlcm_location_t * msg,
                       void * user) {
    
    PLoc* pLoc = (PLoc*)user;
    
    printMyLocation(pLoc);
    printLocation(msg, channel, pLoc);
    
    int pId = (int)msg->pid;

    pLoc->locationMap[pId] = msg; // store the remote process' location info in the map
    
    // measure the delay in which messages get processed by the receiving process
    // (note: for local tests ok. but tests with remote lcm processes the clock might not be in sync!)
    long delay = getSystemTimeMS() - msg->timestampMS;
    if (delay > 10)
    	cout << " WARN1. message delay > 10 sec. " << endl;

	// check for immediate collision (between my current position and the position of the partner process)
	bool collisionCheck = false;
	if( pLoc->pid != (int)msg->pid ) {
		collisionCheck = immediateCollisionCheck(
			&pLoc->myLocation[0],
			&pLoc->locationMap[pId]->position[0],
			pLoc->pid,
			pId,
			msg->timestampLogical
		);
	}
	
	if( collisionCheck ) {
		// TODO: do something. stop, etc.
	}
	
	// check for future likely collisions and react if likely by changing the velocity vector
	futureCollisionCheck(
		&pLoc->myLocation[0],
		&pLoc->myVelocity[0],
		&pLoc->locationMap[pId]->position[0],
		&pLoc->locationMap[pId]->velocity[0],
		pLoc->pid,
		pId,
		msg->timestampLogical
	);
	
	// TODO: only one should correct. idea: only process with lower id changes the course
	
	// TODO: when course correction, use the information of all other processes in the map
	
	// TODO: check for potential collision using the velocity (+ acc.) vectors of other processes

/*    if(!ignoreSelfMessages || pLoc->pid != (int)msg->pid) {
    	printf("Received message on channel \"%s\":\n", channel);
    	printf("  pid                = %"PRId64"\n", msg->pid);
    	printf("  timestampLogical   = %"PRId64"\n", msg->timestampLogical);
    	printf("  timestampMS        = %"PRId64"\n", msg->timestampMS);
    	printf("  message delay [MS] = %"PRId64"\n", delay);
    	printf("  position           = "); print3DArray( pLoc->locationMap[pId]->position ); printf("\n");
    } */

}

void asyncReceive(
	lcm_t* lcm,
	const int pid,
	const int timeoutS,
	const int timeoutMS
) {

    // setup the LCM file descriptor for waiting.
    int lcm_fd = lcm_get_fileno(lcm);
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(lcm_fd, &fds);

    // wait a limited amount of time for an incoming message
    struct timeval timeout = { 
        timeoutS,  // seconds
        timeoutMS  // microseconds
    };
    int status = select(lcm_fd + 1, &fds, 0, 0, &timeout); // http://www.tutorialspoint.com/unix_system_calls/_newselect.htm

    if(0 == status) {
        // no messages
        printf("waiting for message\n");
    } else if(FD_ISSET(lcm_fd, &fds)) {
        // LCM has events ready to be processed.
        lcm_handle(lcm);
    }

}

void asyncReceiveLoop(
	lcm_t* lcm,
	const int pid,
	const int timeoutS,
	const int timeoutMS
) {

	// async event loop: every timeoutS/timeoutMS time period check if there is any
	// message to be processed.
    while(!stopRecv) {
    	asyncReceive(lcm, pid, timeoutS, timeoutMS);
    }

}

// determine our current position (neglecting acceleration, other forces, etc.) by
//  applying our velocity vector to our current position
void updateSendDataLocation(
	exlcm_location_t* send_data,
	PLoc* pLoc,
	int logicalTimeStamp
) {
	// my new position is determined by my latest location + the velocity (which could have been recently updated by the receiving thread to avoid a collision)
	for(int j = 0; j < DIM; ++j)
		send_data->position[j] = pLoc->myLocation[j] + pLoc->myVelocity[j];

	// make sure to also update the shared location object (representing the local position)
	// (this is needed for computing if a collision is happening or not)
	copy3DArray(&send_data->position[0], &pLoc->myLocation[0]);
	
	// make sure to communicate the change in velocity (changed by receiving thread to the outside world)
	copy3DArray(&pLoc->myVelocity[0], &send_data->velocity[0]);
	
	pLoc->logicalTimeStamp = logicalTimeStamp;
}

void updateSendDataTime(
	exlcm_location_t* send_data,
	int timestampLogical
) {
	send_data->timestampLogical = timestampLogical;
	send_data->timestampMS = getSystemTimeMS();
}

exlcm_location_t initSendData(
	int pid,
	PLoc* pLoc
) {

	double initPosition[3] = { (rand()+pid) % 17 + 1, (rand()+pid) % 13 + 1, 0 };
	// double initVelocity[3] = { PlusOrMinusOne(), PlusOrMinusOne(), 0 };
	
	int r1 = rand() + pid;
	int r2 = rand() - pid;
	// double initVelocity[3] = { (r1 % 5) - 2, r2 % 5 - 2, 0 };
	double initVelocity[3] = { (r1 % 3) - 1, r2 % 3 - 1, 0 };

	// copy the init position and velocity into the shared data-structure
	// (we need it there since the receive thread will read from it and might change things if course correction needed)
	copy3DArray(&initPosition[0], &pLoc->myLocation[0]);
	copy3DArray(&initVelocity[0], &pLoc->myVelocity[0]);

 	exlcm_location_t send_data = {
        .timestampLogical = 0,
        .timestampMS = 0,
        .pid = pid,
        .position = {0,0,0},
        .velocity = {0,0,0},
        .destination = { 2, 2, 2 }
    };
    
    copy3DArray(&initPosition[0], &send_data.position[0]);
    copy3DArray(&initVelocity[0], &send_data.velocity[0]);

    return send_data;
}


void sendLoop(
	int tid,
	int pid,
	int N,
	const int sendFSec,
	lcm_t * lcm,
	PLoc* pLoc
) {

	// initialize publish data
	exlcm_location_t send_data = initSendData( pid, pLoc );
    
	cout << " Master with TID = " << tid << " starting publishing messages. PID = " << pid << endl;
			
	// repeat the send/recv procedure for N many times
	for(int i = 0; i < N; ++i) {
	
	    // position data published every "sendFSec" many seconds
    	this_thread::sleep_for(chrono::seconds(sendFSec));

    	updateSendDataLocation( &send_data, pLoc, i );
		updateSendDataTime( &send_data, i );

    	exlcm_location_t_publish(lcm, "LOCATION", &send_data);
		
	}
    this_thread::sleep_for(chrono::seconds(sendFSec));
			
	cout << " Master with TID = " << tid << " stopping receive threads." << endl;
	stopRecv = true;

}

void subscribeToChannels(
	lcm_t * lcm,
	PLoc* pLoc
) {
	// subscribe to the "LOCATION" channel
    exlcm_location_t_subscription_t * sub =
        exlcm_location_t_subscribe(lcm, "LOCATION", &my_handler, pLoc);

}

int main(
	int argc,
	char ** argv
) {

    int pid = getpid(); // process id
    int tid = 0; // openMP thread id
	int N = 10; // how often to publish/receive the message

	if (argc > 1 && sscanf (argv[1], "%i", &N)!=1) {
		N = 10; printf ("error - N not an integer. N will be set to 10.");
	}
	if (argc > 2 && sscanf (argv[2], "%i", &pid)!=1) {
		pid = getpid(); printf ("warn - no process id given. taken system process id.");
	}

    lcm_t * lcm = lcm_create(NULL);
    if(!lcm)
        return 1;
    
    PLoc* pLoc = new PLoc(pid);
    
    subscribeToChannels(lcm, pLoc);
	
	const int timeoutS  = 2;
	const int timeoutMS = 0;//500;
	const int sendFSec = 1;

	// split up in a send / receive thread:
	//  receive thread:         will constantly collect incoming messages (and
	//                          update shared data structures)
	//  send thread (master=0): will publish its position data according to a frequency
	#pragma omp parallel num_threads(2) private(tid) shared(stopRecv, pLoc)
	{
		tid = omp_get_thread_num();
		// master thread does the send
		if( tid == 0 ) {
			sendLoop(tid, pid, N, sendFSec, lcm, pLoc);
    	} else {
			// other threads can receive
			asyncReceiveLoop(lcm, pid, timeoutS, timeoutMS);
    	}
    }

    lcm_destroy(lcm); // clean up any resources used by LCM
    
    delete pLoc;
    
    return 0;

}

