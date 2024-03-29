#ifndef RecoMuon_GlobalTrackingTools_GlobalMuonRefitter_H
#define RecoMuon_GlobalTrackingTools_GlobalMuonRefitter_H

/** \class GlobalMuonRefitter
 *  class to build muon trajectory
 *
 *  $Date: 2010/09/15 09:18:59 $
 *  $Revision: 1.5 $
 *
 *  \author N. Neumeister 	 Purdue University
 *  \author C. Liu 		 Purdue University
 *  \author A. Everett 		 Purdue University
 */

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMuon/TrackingTools/interface/MuonTrajectoryBuilder.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

//for seed regeneration
#include "TrackingTools/DetLayers/interface/NavigationDirection.h"
#include "RecoTracker/TrackProducer/interface/TrackProducerBase.h"
#include "RecoTracker/SpecialSeedGenerators/interface/SeedFromGenericPairOrTriplet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"


namespace edm {class Event;}
namespace reco {class TransientTrack;}

class TrajectoryStateOnSurface;

class MuonDetLayerMeasurements;
class MuonServiceProxy;
class Trajectory;

class TrajectoryFitter;

class GlobalMuonRefitter {

  public:

    typedef TransientTrackingRecHit::RecHitContainer RecHitContainer;
    typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
    typedef TransientTrackingRecHit::RecHitPointer RecHitPointer;
    typedef TransientTrackingRecHit::ConstRecHitPointer ConstRecHitPointer;

    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    typedef std::vector<Trajectory> TC;
    typedef TC::const_iterator TI;

    enum subDetector { PXB = 1, PXF = 2, TIB = 3, TID = 4, TOB = 5, TEC = 6 };

  public:

    /// constructor with Parameter Set and MuonServiceProxy
    GlobalMuonRefitter(const edm::ParameterSet&, const MuonServiceProxy*);
          
    /// destructor
    virtual ~GlobalMuonRefitter();

    /// pass the Event to the algo at each event
    virtual void setEvent(const edm::Event&);

    /// set the services needed by the TrackTransformer
    void setServices(const edm::EventSetup&);

    /// build combined trajectory from sta Track and tracker RecHits
    std::vector<Trajectory> refit(const reco::Track& globalTrack, const int theMuonHitsOption) const;

    /// build combined trajectory from subset of sta Track and tracker RecHits
    std::vector<Trajectory> refit(const reco::Track& globalTrack,
				  const reco::TransientTrack track,
				  TransientTrackingRecHit::ConstRecHitContainer allRecHitsTemp,
				  const int theMuonHitsOption) const;

    /// refit the track with a new set of RecHits
    std::vector<Trajectory> transform(const reco::Track& newTrack,
                                      const reco::TransientTrack track,
                                      TransientTrackingRecHit::ConstRecHitContainer recHitsForReFit) const;
    
    // get rid of selected station RecHits
    ConstRecHitContainer getRidOfSelectStationHits(ConstRecHitContainer hits) const;

    //split hits in even or odd
    ConstRecHitContainer getEvenOddHits(ConstRecHitContainer hits, bool const getEven) const;

    //additional to create new seeds
    TrajectoryStateOnSurface TSOSFromNewSeed(
        TrajectorySeed* newseed,ConstRecHitContainer newseedtrackerhits,
        const MagneticField* myMF) const;
    TrajectorySeed* NewSeedFromPairOrTriplet(
        TransientTrackingRecHit::ConstRecHitContainer& recHitsForReFit, 
        PropagationDirection& propDir,ConstRecHitContainer& newseedtrackerhits 
        ) const;

    //jitter TSOS
    TrajectoryStateOnSurface scaleTSOS(TrajectoryStateOnSurface tsosIn, 
                                                           double const scale) const;
    




  protected:

    enum RefitDirection{RinsideOut,RoutsideIn,Rundetermined};
    
    /// check muon RecHits, calculate chamber occupancy and select hits to be used in the final fit
    void checkMuonHits(const reco::Track&, ConstRecHitContainer&, 
                       std::vector<int>&) const;

    /// get the RecHits in the tracker and the first muon chamber with hits 
    void getFirstHits(const reco::Track&, ConstRecHitContainer&, 
                       ConstRecHitContainer&) const;
 
    /// select muon hits compatible with trajectory; check hits in chambers with showers
    ConstRecHitContainer selectMuonHits(const Trajectory&, 
                                        const std::vector<int>&) const;
 
    /// print all RecHits of a trajectory
    void printHits(const ConstRecHitContainer&) const;

    RefitDirection checkRecHitsOrdering(const ConstRecHitContainer&) const;

    const MuonServiceProxy* service() const { return theService; }

  protected:
    std::string theCategory;
    bool theTkTrajsAvailableFlag;
    float thePtCut;

  private:
  
    int   theMuonHitsOption;
    float theProbCut;
    int   theHitThreshold;
    float theDTChi2Cut;
    float theCSCChi2Cut;
    float theRPCChi2Cut;
    bool  theCosmicFlag;
    std::string theHitsToKeep;
    std::string theNewSeed;
    double theJitterScale;
    

    edm::InputTag theDTRecHitLabel;
    edm::InputTag theCSCRecHitLabel;
    edm::Handle<DTRecHitCollection>    theDTRecHits;
    edm::Handle<CSCRecHit2DCollection> theCSCRecHits;

    int	  theSkipStation;
    int   theTrackerSkipSystem;
    int   theTrackerSkipSection;

    unsigned long long theCacheId_TRH;        

    std::string thePropagatorName;
  
    bool theRPCInTheFit;

    RefitDirection theRefitDirection;

    std::string theFitterName;
    edm::ESHandle<TrajectoryFitter> theFitter;
  
    std::string theTrackerRecHitBuilderName;
    edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  
    std::string theMuonRecHitBuilderName;
    edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;

    const MuonServiceProxy* theService;
    const edm::Event* theEvent;
};
#endif
