//
// Created by sylwester on 7/16/19.
//

#ifndef EXACTVERTEXCOVER_MYPQ_H
#define EXACTVERTEXCOVER_MYPQ_H

#include <utils/TimeMeasurer.h>
#include "Makros.h"

using namespace std;

/**
 * My version of priority queue (or rather a heap). It enables changing the value of an element and simultaneously
 * rearranging it in the priority queue.
 *
 * To set value of element with given id to val, use setPriority( id,val ) function.
 * To extract id of an element with maximum value, use top() function.
 * To get the value of an element with given id, use getPriority(id)
 * THAT'S IT!
 *
 *  This class may be seen as an array, in which modification of an element is done pessimistically in lgN time, but we
 *  have acces to the largest element in O(1) time.
 */
 template<class _T>
class MyPQ {
public:

    /**
     * Helper class. Each item in the queue is stored in an array (to enable quick access via index) and in queue via pointer to that item.
     */
    class Item{
        public:
            Item( int id, _T prior, int indInQueue ) : id(id), priority(prior), index(indInQueue) {}
            int id;
            _T priority;
            int index;
    };
    // END OF ITEM CLASS

    /**
     * Creates priority queue for N elements. Their initial value is set to initVal
     * @param N
     */
    MyPQ( int N, _T initVal ) : N(N), validN(N) {
        for(int i=0; i<N; i++){
            items.push_back( Item(i,initVal,i) );
        }

        itemQueue = vector<Item*>(N);
        for(int i=0; i<N; i++){
            itemQueue[i] = &items[i];
        }
    }


    /**
     * Function used whenever a value of element is changed. It sifts up or down the element in the queue.
     * Worst case complexity: O( lgN )
     * @param id
     * @param prior
     * @return true if the value was changed, false otherwise
     */
    bool setPriority(int id, _T prior) {

        if( items[id].priority == prior ) return false;
        else{
            items[id].priority = prior;
            if( !shiftUp(id) ) shiftDown(id);

            return true;
        }
    }


    /**
     * Function does not remove anything, just returns the id of the top of the queue.
     * @return id of the element with graetest priority
     */
    int top() {
        return itemQueue[0]->id;
    }

    /*_T remove() {
        int a = 0;
        int b = validN-1;
        int id = itemQueue[b]->id;
        _T prior = itemQueue[b]->priority;

        swapElements(a,b);
        validN--;
        setPriority(id, prior);
    }*/

    /**
     *
     * @param id
     * @return priority of element with gven id
     */
    _T getPriority( int id ){ return items[id].priority; }

    int size(){ return N; }



    static void test(){

        int l_zest = 20;
        for(int lz = 0; lz<l_zest;lz++) {


            int N = 1000;
            MyPQ<PII> pq(N,{0,0});

            int ITER = N * N;
            VPII values(N, {0,0});

            priority_queue<PII> pqstl;
            for (int i = 0; i < N; i++) pqstl.push({0,0});

            for (int i = 0; i < ITER; i++) {
//                cerr << "\r" << ((double) i / ITER) << flush;

                int ind = rand() % N;
                PII prior = {rand() % (20), rand() % (20) };

//        cerr << "setting priority of " << ind << " to " << prior << endl;

                TimeMeasurer::startMeasurement("MyPQ");
                pq.setPriority(ind, prior);
                TimeMeasurer::stopMeasurement("MyPQ");

                TimeMeasurer::startMeasurement("pqstl");
                pqstl.pop();
                pqstl.push(prior);
                TimeMeasurer::stopMeasurement("pqstl");


                values[ind] = prior;

//        for (int k = 0; k < N; k++) {
//            cerr << k << " -> (prior, index) = (" << pq.items[k].priority << "," << pq.items[k].index << ")" << endl;
//        }
//        cerr << endl;
//
//        for( int j=0; j<pq.validN; j++ ){
//            cerr << "itemQueue[" << j << "].id: " << pq.itemQueue[j]->id << "   prior: " << pq.itemQueue[j]->priority << endl;
//        }

//                pq.isHeap();

                if (pq.getPriority(pq.top()) != *max_element(ALL(values))) {
                    cerr << "pq max: " << pq.getPriority(pq.top()) << endl;
                    cerr << "values max: " << *max_element(ALL(values)) << endl;

//                    cerr << "pq values: " << endl;
//                    for (int k = 0; k < N; k++) {
//                        cerr << k << " -> (prior, index) = (" << pq.items[k].priority << "," << pq.items[k].index << ")"
//                             << endl;
//                    }
                    exit(1);
                }


//        cerr << "\ttop id: " << pq.top() << "    priority: " << pq.getPriority(pq.top()) << endl;
//        ENDL(2);

            }

            cerr << "\r" << ((double) lz / l_zest) << flush;
        }


        cerr << "TEST PASSED!" << endl;

        TimeMeasurer::writeAllMeasurements();
    }


private:

    /**
     * Debug function.
     * @return true if itemQueue forms a heap, false otherwise
     */
    bool isHeap(){
        for( int i=0; i<validN; i++ ){
            if( (i << 1) >= validN ) break;
            else{
                if( itemQueue[i<<1]->priority > itemQueue[i]->priority ) return false;

                if( (i<<1)+1 >= validN ) continue;
                else{
                    if( itemQueue[(i<<1) + 1]->priority > itemQueue[i]->priority ) return false;
                }
            }
        }
        return true;
    }

    /**
     *
     * @param id
     * @return true if the element was shifted up, false otherwise
     */
    bool shiftUp(int id) {
        Item & it = items[id];
        int cnt = 0;
        while( it.index != 0 && itemQueue[ it.index >> 1 ]->priority < it.priority ){
            int a = it.index;
            int b = it.index >> 1;

            itemQueue[ b ]->index = a;
            it.index = b;

            swap( itemQueue[ a ], itemQueue[ b ]  ); // swapping POINTERS! Items only change index value, do not move anywhere :)
            cnt++;
        }
        return cnt != 0;
    }

    /**
     * Tries to sift down the element with given id
     * @param id
     * @return true if the element was shifted down, false otherwise
     */
    bool shiftDown(int id) {
        Item & it = items[id];
        int cnt = 0;
        int indToSwap = it.index;

        while( it.index << 1 < validN ){

            int L = ( it.index << 1 );
            int R = ( it.index << 1 ) + 1;

            if( L < validN && itemQueue[L]->priority > it.priority ) indToSwap = L;

            if( R < validN && itemQueue[R]->priority > itemQueue[indToSwap]->priority ) indToSwap = R;

            if( indToSwap == it.index ) return false;
            else{
                swapElements( it.index, indToSwap );
            }

        }
    }

    /**
     * Function swaps the elements on given two position in the priority queue
     * @param a
     * @param b
     */
    void swapElements( int posInQueue1, int posInQueue2 ){
        itemQueue[posInQueue1]->index = posInQueue2;
        itemQueue[posInQueue2]->index = posInQueue1;

        swap( itemQueue[posInQueue1], itemQueue[posInQueue2] );
    }

    int N; // capacity of the queue
    int validN; // elements that are in the queue, that were not removed.
    vector<Item> items;
    vector<Item*> itemQueue;



};


#endif //EXACTVERTEXCOVER_MYPQ_H
