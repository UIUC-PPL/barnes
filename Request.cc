/*
 * CharmBH: Request.cc
 * Helps DataManager (DM) manage requests for remote data from TreePieces 
 * on its PE. Requests are maintained in hashtables mapping requested Keys
 * to Request records. The methods below embody actions to be taken 
 * upon receipt of remote data.
 */

#include "Request.h"

#include "State.h"
#include "Worker.h"
#include "Traversal_decls.h"
#include "Node.h"

#include "TreePiece.h"

/*
 * When remote particles are received, we must supply them to all workers
 * on this PE that requested them. These workers will then use them for
 * computations that couldn't previously perform due to missing data.
 */
void Request::deliverParticles(int num){
  // For each requestor of this set of particles
  for(int i = 0; i < requestors.length(); i++){
    Requestor &req = requestors[i];
    // Which worker requested these particles?
    CutoffWorker<ForceData> *worker = req.worker;

    /*
     * Workers maintain the current local bucket for which they
     * are performing the traversal (either remote or local). 
     * Currently, the worker could be in the middle of a traversal
     * for some other local bucket. We save the pointer to this 
     * local bucket here, and replace it with the one which caused
     * the worker to request these particles in the first place.
     */
    void *saveContext = worker->getContext();
    worker->setContext(req.context);

    ExternalParticle *externalParticles = (ExternalParticle *)data;
    // so that the worker may do some book-keeping, if required 
    for(int j = 0; j < num; j++){ 
      // supply each particle in the bucket to the worker
      worker->work(externalParticles+j);
    }

    // Restore the current bucket of the worker, in case
    // it is in the middle of some other traversal.
    worker->setContext(saveContext);
    // Update the corresponding state to reflect the fact that a previously
    // outstanding remote data request is now complete.
    State *state = req.state;
    state->decrPending();
    // If there are no more outstanding data requests for the 
    // traversal, we are done with it.
    if(state->complete()) worker->done();
  }

  // No more requestors of this bucket to keep track of 
  requestors.clear();
}

/*
 * When the subtree under a requested node is obtained on the PE, find the 
 * workers that requested it, and perform a traversal on this subtree
 * for each one. Note that the worker may be in the midst of a traversal for
 * a different bucket from the one that generated this request. Therefore,
 * we must save and restore its current bucket through get/setContext.
 */
void Request::deliverNode(){
  for(int i = 0; i < requestors.length(); i++){
    Requestor &req = requestors[i];
    CutoffWorker<ForceData> *worker = req.worker;

    // Save a pointer to the current bucket of the worker.
    void *saveContext = worker->getContext();
    // Replace it with the one that generated the request
    // for this subtree
    worker->setContext(req.context);

    // Recall the kind of traversal the worker was performing
    // when it generated the request... 
    Traversal<ForceData> *traversal = req.traversal; 
    // And a pointer to its current book-keeping state
    State *state = req.state;
    // For each child of the parent whose children were requested...
    Node<ForceData> *firstChild = (Node<ForceData>*) data;
    for(int j = 0; j < BRANCH_FACTOR; j++){
      // perform a top-down traversal of this subtree.
      traversal->topDownTraversal(firstChild+j,worker,state);
    }

    // Restore the current bucket of the worker
    worker->setContext(saveContext);
    // Check whether this traversal has completed.
    state->decrPending();
    if(state->complete()) worker->done();
  }

  requestors.clear();

}

#include "Traversal_defs.h"
