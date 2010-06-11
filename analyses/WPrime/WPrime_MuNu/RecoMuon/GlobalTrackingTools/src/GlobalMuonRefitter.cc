/**
 *  Class: GlobalMuonRefitter
 *
 *  Description:
 *
 *
 *  $Date: 2009/09/12 20:33:33 $
 *  $Revision: 1.11 $
 *
 *  Authors :
 *  P. Traczyk, SINS Warsaw
 *
 **/

#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"

//---------------
// C++ Headers --
//---------------

#include <iostream>
#include <iomanip>
#include <algorithm>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TrackFitters/interface/RecHitLessByDet.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TRandom3.h"




using namespace std;
using namespace edm;

//----------------
// Constructors --
//----------------

GlobalMuonRefitter::GlobalMuonRefitter(const edm::ParameterSet& par,
				       const MuonServiceProxy* service) : 
  theCosmicFlag(par.getParameter<bool>("PropDirForCosmics")),
  theDTRecHitLabel(par.getParameter<InputTag>("DTRecSegmentLabel")),
  theCSCRecHitLabel(par.getParameter<InputTag>("CSCRecSegmentLabel")),
  theService(service) {

  theCategory = par.getUntrackedParameter<string>("Category", "Muon|RecoMuon|GlobalMuon|GlobalMuonRefitter");

  theHitThreshold = par.getParameter<int>("HitThreshold");
  theDTChi2Cut  = par.getParameter<double>("Chi2CutDT");
  theCSCChi2Cut = par.getParameter<double>("Chi2CutCSC");
  theRPCChi2Cut = par.getParameter<double>("Chi2CutRPC");

  // Refit direction
  string refitDirectionName = par.getParameter<string>("RefitDirection");
  //cout<<"refitDirectionName = "<<refitDirectionName<<endl;
  
  LogDebug(theCategory)<<"GloblaMuonRefitter constructor called"<<endl;
  

  if (refitDirectionName == "RinsideOut" ) theRefitDirection = RinsideOut;
  else if (refitDirectionName == "RoutsideIn" ) theRefitDirection = RoutsideIn;
  else 
    throw cms::Exception("TrackTransformer constructor") 
      <<"Wrong refit direction chosen in TrackTransformer ParameterSet"
      << "\n"
      << "Possible choices are:"
      << "\n"
      << "RefitDirection = RinsideOut or RefitDirection = RoutsideIn";
  
  theFitterName = par.getParameter<string>("Fitter");  
  thePropagatorName = par.getParameter<string>("Propagator");

  theSkipStation        = par.getParameter<int>("SkipStation");
  theTrackerSkipSystem	= par.getParameter<int>("TrackerSkipSystem");
  theTrackerSkipSection	= par.getParameter<int>("TrackerSkipSection");//layer, wheel, or disk depending on the system

  theTrackerRecHitBuilderName = par.getParameter<string>("TrackerRecHitBuilder");
  theMuonRecHitBuilderName = par.getParameter<string>("MuonRecHitBuilder");
  theStab 		= par.getUntrackedParameter<int>("StabStab", 1);
  theStabScale 	= par.getUntrackedParameter<double>("StabScale", 1);
  theNewSeed 	= par.getUntrackedParameter<bool>("NewSeed", false);

  theRPCInTheFit = par.getParameter<bool>("RefitRPCHits");

  theCacheId_TRH = 0;

}

//--------------
// Destructor --
//--------------

GlobalMuonRefitter::~GlobalMuonRefitter() {
}


//
// set Event
//
void GlobalMuonRefitter::setEvent(const edm::Event& event) {

  event.getByLabel(theDTRecHitLabel, theDTRecHits);
  event.getByLabel(theCSCRecHitLabel, theCSCRecHits);
}


void GlobalMuonRefitter::setServices(const EventSetup& setup) {

  theService->eventSetup().get<TrajectoryFitter::Record>().get(theFitterName,theFitter);

  // Transient Rechit Builders
  unsigned long long newCacheId_TRH = setup.get<TransientRecHitRecord>().cacheIdentifier();
  if ( newCacheId_TRH != theCacheId_TRH ) {
    LogDebug(theCategory) << "TransientRecHitRecord changed!";
    setup.get<TransientRecHitRecord>().get(theTrackerRecHitBuilderName,theTrackerRecHitBuilder);
    setup.get<TransientRecHitRecord>().get(theMuonRecHitBuilderName,theMuonRecHitBuilder);
  }
}


//
// build a combined tracker-muon trajectory
//
vector<Trajectory> GlobalMuonRefitter::refit(const reco::Track& globalTrack, 
                                             const int theMuonHitsOption) const {
  LogTrace(theCategory) << " *** GlobalMuonRefitter *** option " << theMuonHitsOption << endl;
    
  ConstRecHitContainer allRecHitsTemp; // all muon rechits temp

  reco::TransientTrack track(globalTrack,&*(theService->magneticField()),theService->trackingGeometry());
  
  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit)
      if((*hit)->isValid())
          if ( (*hit)->geographicalId().det() == DetId::Tracker )
              allRecHitsTemp.push_back(theTrackerRecHitBuilder->build(&**hit));
          else if ( (*hit)->geographicalId().det() == DetId::Muon ){
              if( (*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit){
                  LogTrace(theCategory) << "RPC Rec Hit discarged"; 
                  continue;
              }
              allRecHitsTemp.push_back(theMuonRecHitBuilder->build(&**hit));
          }
          
  

  vector<Trajectory> refitted = refit(globalTrack,track,allRecHitsTemp,theMuonHitsOption);
  return refitted;
}

//
// build a combined tracker-muon trajectory
//
vector<Trajectory> GlobalMuonRefitter::refit(const reco::Track& globalTrack,
					     const reco::TransientTrack track,
					     TransientTrackingRecHit::ConstRecHitContainer allRecHitsTemp,
                                             const int theMuonHitsOption) const {

  // MuonHitsOption: 0 - tracker only
  //                 1 - include all muon hits
  //                 2 - include only first muon hit(s)
  //                 3 - include only selected muon hits

  vector<int> stationHits(4,0);

  ConstRecHitContainer allRecHits; // all muon rechits
  ConstRecHitContainer fmsRecHits; // only first muon rechits
  ConstRecHitContainer selectedRecHits; // selected muon rechits

   LogTrace(theCategory) << " *** GlobalMuonRefitter *** option " << theMuonHitsOption << endl;

   LogTrace(theCategory) << " Track momentum before refit: " << globalTrack.pt() << endl;

  LogTrace(theCategory) << " Hits size before : " << allRecHitsTemp.size() << endl;
  printHits(allRecHitsTemp);
  allRecHits = getRidOfSelectStationHits(allRecHitsTemp);  

  //std::cout << " Hits size before : " << allRecHitsTemp.size() << endl;
  //std::cout<<" theStab = "<<theStab<<std::endl;
  //ending in 2 is event, ending in 3 is odd.
	if (theStab ==12) allRecHits = getEvenOddHits(allRecHitsTemp, true);  
	if (theStab ==13) allRecHits = getEvenOddHits(allRecHitsTemp, false);  
	if (theStab ==121) allRecHits = getEvenOddLayers(allRecHitsTemp, true);  
	if (theStab ==131) allRecHits = getEvenOddLayers(allRecHitsTemp, false);  
    //std::cout << " Hits size after : " << allRecHits.size() << endl;

  //    printHits(allRecHits);
  LogTrace(theCategory) << " Hits size after: " << allRecHits.size() << endl;

  vector <Trajectory> outputTraj;

  if ((theMuonHitsOption == 1) || (theMuonHitsOption == 3)) {
    // refit the full track with all muon hits
      vector <Trajectory> globalTraj = transform(globalTrack, track, allRecHits);

    if (!globalTraj.size()) {
      LogTrace(theCategory) << "No trajectory from the TrackTransformer!" << endl;
      return vector<Trajectory>();
    }

    LogTrace(theCategory) << " Initial trajectory state: " 
                          << globalTraj.front().lastMeasurement().updatedState().freeState()->parameters() << endl;
  
    if (theMuonHitsOption == 1 )
      outputTraj.push_back(globalTraj.front());
    
    if (theMuonHitsOption == 3 ) {
      checkMuonHits(globalTrack, allRecHits, stationHits);
      selectedRecHits = selectMuonHits(globalTraj.front(),stationHits);
      LogTrace(theCategory) << " Selected hits size: " << selectedRecHits.size() << endl;  
      outputTraj = transform(globalTrack, track, selectedRecHits);
    }     
  } else if (theMuonHitsOption == 2 )  {
      getFirstHits(globalTrack, allRecHits, fmsRecHits);
      outputTraj = transform(globalTrack, track, fmsRecHits);
    } 

  if (outputTraj.size()) {
    LogTrace(theCategory) << "Refitted pt: " << outputTraj.front().firstMeasurement().updatedState().globalParameters().momentum().perp() << endl;
    return outputTraj;
  } else {
    LogTrace(theCategory) << "No refitted Tracks... " << endl;
    return vector<Trajectory>();
  }
  
}


//
//
//
void GlobalMuonRefitter::checkMuonHits(const reco::Track& muon, 
				       ConstRecHitContainer& all,
				       std::vector<int>& hits) const {

  LogTrace(theCategory) << " GlobalMuonRefitter::checkMuonHits " << endl;

  float coneSize = 20.0;
  int dethits[4];
  for ( int i=0; i<4; i++ ) hits[i]=dethits[i]=0;

  // loop through all muon hits and calculate the maximum # of hits in each chamber
  for (ConstRecHitContainer::const_iterator imrh = all.begin(); imrh != all.end(); imrh++ ) {
        
    if ( (*imrh != 0 ) && !(*imrh)->isValid() ) continue;
  
    int station = 0;
    int detRecHits = 0;
    MuonRecHitContainer dRecHits;
      
    DetId id = (*imrh)->geographicalId();

    // Skip tracker hits
    if (id.det()!=DetId::Muon) continue;

    if ( id.subdetId() == MuonSubdetId::DT ) {
      DTChamberId did(id.rawId());
      DTLayerId lid(id.rawId());
      station = did.station();

      // Get the 1d DT RechHits from this layer
      DTRecHitCollection::range dRecHits = theDTRecHits->get(lid);

      for (DTRecHitCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
	double rhitDistance = fabs(ir->localPosition().x()-(**imrh).localPosition().x());
	if ( rhitDistance < coneSize ) detRecHits++;
        LogTrace(theCategory)	<< "       " << (ir)->localPosition() << "  " << (**imrh).localPosition()
               << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
      }
    }// end of if DT
    else if ( id.subdetId() == MuonSubdetId::CSC ) {
    
      CSCDetId did(id.rawId());
      station = did.station();

      // Get the CSC Rechits from this layer
      CSCRecHit2DCollection::range dRecHits = theCSCRecHits->get(did);      

      for (CSCRecHit2DCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
	double rhitDistance = (ir->localPosition()-(**imrh).localPosition()).mag();
	if ( rhitDistance < coneSize ) detRecHits++;
        LogTrace(theCategory)	<< ir->localPosition() << "  " << (**imrh).localPosition()
	       << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
      }
    }
    else {
      if ( id.subdetId() != MuonSubdetId::RPC ) LogError(theCategory)<<" Wrong Hit Type ";
      continue;      
    }
      
    if ( (station > 0) && (station < 5) ) {
      if ( detRecHits > hits[station-1] ) hits[station-1] = detRecHits;
    }

  } // end of loop over muon rechits

  for ( int i = 0; i < 4; i++ ) 
    LogTrace(theCategory) <<" Station "<<i+1<<": "<<hits[i]<<" "<<dethits[i] <<endl; 

  LogTrace(theCategory) << "CheckMuonHits: "<<all.size();

  // check order of muon measurements
  if ( (all.size() > 1) &&
       ( all.front()->globalPosition().mag() >
	 all.back()->globalPosition().mag() ) ) {
    LogTrace(theCategory)<< "reverse order: ";
    stable_sort(all.begin(),all.end(),RecHitLessByDet(alongMomentum));
  }
}


//
// Get the hits from the first muon station (containing hits)
//
void GlobalMuonRefitter::getFirstHits(const reco::Track& muon, 
				       ConstRecHitContainer& all,
				       ConstRecHitContainer& first) const {
    
    LogTrace(theCategory) << " GlobalMuonRefitter::getFirstHits\nall rechits length:" << all.size() << endl;
    first.clear();
    
    int station_to_keep = 999;
    vector<int> stations;
    for (ConstRecHitContainer::const_iterator ihit = all.begin(); ihit != all.end(); ++ihit) {
        
        int station = 0;
        bool use_it = true;
        DetId id = (*ihit)->geographicalId();
        unsigned raw_id = id.rawId();
        if (!(*ihit)->isValid()) station = -1;
        else {
            if (id.det() == DetId::Muon) {
                switch (id.subdetId()) {
                    case MuonSubdetId::DT:  station = DTChamberId(raw_id).station(); break;
                    case MuonSubdetId::CSC: station = CSCDetId(raw_id).station(); break;
                    case MuonSubdetId::RPC: station = RPCDetId(raw_id).station(); use_it = false; break;
                }
            }
        }
        
        if (use_it && station > 0 && station < station_to_keep) station_to_keep = station;
        stations.push_back(station);
        
        LogTrace(theCategory) << "rawId: " << raw_id << " station = " << station << " station_to_keep is now " << station_to_keep;
    }
    
    if (station_to_keep <= 0 || station_to_keep > 4 || stations.size() != all.size())
        LogInfo(theCategory) << "failed to getFirstHits (all muon hits are outliers/bad ?)! station_to_keep = " 
                             << station_to_keep << " stations.size " << stations.size() << " all.size " << all.size();
    
    for (unsigned i = 0; i < stations.size(); ++i)
        if (stations[i] >= 0 && stations[i] <= station_to_keep) first.push_back(all[i]);
    
    return;
}


//
// select muon hits compatible with trajectory; 
// check hits in chambers with showers
//
GlobalMuonRefitter::ConstRecHitContainer 
GlobalMuonRefitter::selectMuonHits(const Trajectory& traj, 
                                   const std::vector<int>& hits) const {

  ConstRecHitContainer muonRecHits;
  const double globalChi2Cut = 200.0;

  vector<TrajectoryMeasurement> muonMeasurements = traj.measurements(); 

  // loop through all muon hits and skip hits with bad chi2 in chambers with high occupancy      
  for (std::vector<TrajectoryMeasurement>::const_iterator im = muonMeasurements.begin(); im != muonMeasurements.end(); im++ ) {

    if ( !(*im).recHit()->isValid() ) continue;
    if ( (*im).recHit()->det()->geographicalId().det() != DetId::Muon ) {
      //      if ( ( chi2ndf < globalChi2Cut ) )
      muonRecHits.push_back((*im).recHit());
      continue;
    }  
    ConstMuonRecHitPointer immrh = dynamic_cast<const MuonTransientTrackingRecHit*>((*im).recHit().get());

    DetId id = immrh->geographicalId();
    int station = 0;
    int threshold = 0;
    double chi2Cut = 0.0;

    // get station of hit if it is in DT
    if ( (*immrh).isDT() ) {
      DTChamberId did(id.rawId());
      station = did.station();
      threshold = theHitThreshold;
      chi2Cut = theDTChi2Cut;
    }
    // get station of hit if it is in CSC
    else if ( (*immrh).isCSC() ) {
      CSCDetId did(id.rawId());
      station = did.station();
      threshold = theHitThreshold;
      chi2Cut = theCSCChi2Cut;
    }
    // get station of hit if it is in RPC
    else if ( (*immrh).isRPC() ) {
      RPCDetId rpcid(id.rawId());
      station = rpcid.station();
      threshold = theHitThreshold;
      chi2Cut = theRPCChi2Cut;
    }
    else
      continue;

    double chi2ndf = (*im).estimate()/(*im).recHit()->dimension();  

    bool keep = true;
    if ( (station>0) && (station<5) ) {
      if (hits[station-1]>threshold) keep = false;
    }   
    
    if ( (keep || (chi2ndf<chi2Cut)) && (chi2ndf<globalChi2Cut) ) {
      muonRecHits.push_back((*im).recHit());
    } else {
      LogTrace(theCategory)
	<< "Skip hit: " << id.det() << " " << station << ", " 
	<< chi2ndf << " (" << chi2Cut << " chi2 threshold) " 
	<< hits[station-1] << endl;
    }
  }
  
  // check order of rechits
  reverse(muonRecHits.begin(),muonRecHits.end());
  return muonRecHits;
}


//
// print RecHits
//
void GlobalMuonRefitter::printHits(const ConstRecHitContainer& hits) const {

  LogTrace(theCategory) << "Used RecHits: " << hits.size();
  for (ConstRecHitContainer::const_iterator ir = hits.begin(); ir != hits.end(); ir++ ) {
    if ( !(*ir)->isValid() ) {
      LogTrace(theCategory) << "invalid RecHit";
      continue; 
    }
    
    const GlobalPoint& pos = (*ir)->globalPosition();
    
    LogTrace(theCategory) 
      << "r = " << sqrt(pos.x() * pos.x() + pos.y() * pos.y())
      << "  z = " << pos.z()
      << "  dimension = " << (*ir)->dimension()
      << "  det = " << (*ir)->det()->geographicalId().det()
      << "  subdet = " << (*ir)->det()->subDetector()
      << "  raw id = " << (*ir)->det()->geographicalId().rawId();
  }

}


//
// add Trajectory* to TrackCand if not already present
//
GlobalMuonRefitter::RefitDirection
GlobalMuonRefitter::checkRecHitsOrdering(const TransientTrackingRecHit::ConstRecHitContainer& recHits) const {

  if (!recHits.empty()){
    ConstRecHitContainer::const_iterator frontHit = recHits.begin();
    ConstRecHitContainer::const_iterator backHit  = recHits.end() - 1;
    while( !(*frontHit)->isValid() && frontHit != backHit) {frontHit++;}
    while( !(*backHit)->isValid() && backHit != frontHit)  {backHit--;}

    double rFirst = (*frontHit)->globalPosition().mag();
    double rLast  = (*backHit) ->globalPosition().mag();

    if(rFirst < rLast) return RinsideOut;
    else if(rFirst > rLast) return RoutsideIn;
    else {
      LogError(theCategory) << "Impossible determine the rechits order" <<endl;
      return Rundetermined;
    }
  } else {
    LogError(theCategory) << "Impossible determine the rechits order" <<endl;
    return Rundetermined;
  }
}


//-----------------------------------------------------------------------
TrajectoryStateOnSurface GlobalMuonRefitter::TSOSFromNewSeed(
    TrajectorySeed* newseed,ConstRecHitContainer newseedtrackerhits,
    const MagneticField* myMF
    ) const
{
//-----------------------------------------------------------------------
    LogDebug(theCategory)<< "==============Satarting the new TSOS sequence"<<endl;
    TrajectoryStateOnSurface badgerTSOS = 0;
      
    PTrajectoryStateOnDet ptbadgerTSOS = newseed->startingState();
    LogDebug(theCategory) <<"got the newseed starting state successfully"<<endl;
    
    TrajectoryStateTransform mytrajsformer;
    LogDebug(theCategory) <<"got my trajectory transformer"<<endl;
    badgerTSOS = mytrajsformer.transientState(ptbadgerTSOS,newseedtrackerhits.front()->surface(),&*myMF);
    
    LogDebug(theCategory) <<"TrajectoryStateOnSurface was created"<<endl;

    return badgerTSOS;
    
}//------------------------------- TSOSFromNewSeed()




//-----------------------------------------------------------------------
TrajectorySeed* GlobalMuonRefitter::NewSeedFromPairOrTriplet(
    TransientTrackingRecHit::ConstRecHitContainer& recHitsForReFit, 
    PropagationDirection& propDir,ConstRecHitContainer& newseedtrackerhits 
    ) const
{
//----------------------------------------------------------------------
    LogDebug(theCategory)<< "==============Starting the New Seed sequence"<<endl; 

   //my brand new seed
    TrajectorySeed* newseed = 0;

    //Get stuff needed from the event setup
    edm::ESHandle<MagneticField> theMF = theService->magneticField();
    edm::ESHandle<TrackerGeometry> theG; 
    theService->eventSetup().get<TrackerDigiGeometryRecord>().get(theG);
    edm::ESHandle<Propagator> thePropagator = theService->propagator(thePropagatorName);
    edm::ESHandle<Propagator> theRevPropagator = theService->propagator(thePropagatorName);
    const TrackerGeometry* myG = theG.product();
    const MagneticField* myMF = theMF.product();
    const Propagator* myPropagator = thePropagator.product();
    const Propagator* myRevPropagator = theRevPropagator.product();
    //dummy vector of charges
    vector<int> dummycharges;
    
    //Try to use the code to make a new seed from hits
    SeedFromGenericPairOrTriplet seedObjt(&*myMF,&*myG,&*theTrackerRecHitBuilder,&*myPropagator,&*myRevPropagator,dummycharges,false, 1000.);
    
    //Get the first 2 or 3 hits from the tracker to be used in
    //the new seed
    newseedtrackerhits.clear();
    LogTrace(theCategory) <<"recHitsForReFit.size() = "
                          <<recHitsForReFit.size()<<endl;
    LogDebug(theCategory) << "=========About to get hits for my new seed"
                          <<endl;
    for (ConstRecHitContainer::const_iterator it = recHitsForReFit.begin(),
             itend = recHitsForReFit.end();it!=itend;++it){
        if ((*it)->geographicalId().det() == DetId::Muon){ 
            LogDebug(theCategory) <<"Found a Muon hit for seed hits :("<<endl;
        }
        else if ((*it)->geographicalId().det() == DetId::Tracker){
            LogDebug(theCategory) << "Found a Tracker hit for seed hits"<<endl;
            newseedtrackerhits.push_back(*it);
          }
        else { 
            LogDebug(theCategory) <<"WARNING: Neither tracker nor muon hit"
                                  <<endl;
        }
        //break the lopp if max 3 hits from tracker were found
        if (newseedtrackerhits.size() == 3) break;
    }//------- loop over recHitsForReFit
    
    //if there are not enough tracker hits, return zero
    if (newseedtrackerhits.size()<2){
        LogTrace(theCategory)<<"Size of newseedtrackerhits is = "
                             <<newseedtrackerhits.size()
                             <<"; too few tracker hits to make a new seed"
                             <<endl;
        return newseed;
    }
    
    LogTrace(theCategory)<<"Size of newseedtrackerhits is = "
                         <<newseedtrackerhits.size()<<endl;
    printHits(newseedtrackerhits);
    
    //make a new set of seeds for seed generation
    SeedingHitSet newseedset(newseedtrackerhits);
    //Get the Navigation direction (this is very annoying, I had to 
    //rename the direction of propagation enum to avoid conflict)
    NavigationDirection seedDir;
    if (theRefitDirection == RinsideOut ) seedDir = insideOut;
    else if (theRefitDirection == RoutsideIn ) seedDir = outsideIn;
    else 
        throw cms::Exception("RefitDirection problem") 
            <<"Something went wrong with the NavigationDirection";
    LogDebug(theCategory) <<"================ About to create a new seed"
                          <<endl;
    //try to create the new seed
    std::vector<TrajectorySeed*> newseedvect = seedObjt.seed(newseedset,propDir,
                            seedDir,theService->eventSetup());
    if (newseedvect.empty()) {
        LogDebug(theCategory) <<"No new seeds could be created"
                              <<endl;
        return newseed;
    }
    else{
        LogDebug(theCategory) <<"newseedvect.size() = "
                              <<newseedvect.size()<<endl;
        newseed = newseedvect.front();
    }
    
    
    
    //check if the trajectory seed was created successfully
    if (newseed != 0){
        LogDebug(theCategory) <<"New trajectory seed created successfully"
                              <<endl;
        //check the new seed information
        LogTrace(theCategory) << "                      Starting newseed: "
                              << " nHits= " << newseed->nHits()
                              << " tsos: "
                              << newseed->startingState().parameters().position() 
                              << "  p="
                              << newseed->startingState().parameters().momentum() 
                              << endl;
    
    }
    else{
        LogDebug(theCategory) <<"Failed to create a good new trajectory seed"
                              <<endl;
        return newseed;
    }
   
    return newseed;
    

}//--------------------NewSeedFromPairOrTriplet





//
// Convert Tracks into Trajectories with a given set of hits
//
vector<Trajectory> GlobalMuonRefitter::transform(const reco::Track& newTrack,
                                                 const reco::TransientTrack track,
                                                 TransientTrackingRecHit::ConstRecHitContainer recHitsForReFit) const {

    LogDebug(theCategory) <<"Going into the transform fucntion"<<endl;
    
  LogTrace(theCategory) << "GlobalMuonRefitter::transform: " << recHitsForReFit.size() << " hits:";
  printHits(recHitsForReFit);

  if(recHitsForReFit.size() < 2) return vector<Trajectory>();

  // Check the order of the rechits
  RefitDirection recHitsOrder = checkRecHitsOrdering(recHitsForReFit);

  LogTrace(theCategory) << "checkRecHitsOrdering() returned " << recHitsOrder
			<< ", theRefitDirection is " << theRefitDirection
			<< " (RinsideOut == " << RinsideOut << ", RoutsideIn == " << RoutsideIn << ")";

  // Reverse the order in the case of inconsistency between the fit direction and the rechit order
  if(theRefitDirection != recHitsOrder) reverse(recHitsForReFit.begin(),recHitsForReFit.end());

  // Even though we checked the rechits' ordering above, we may have
  // already flipped them elsewhere (getFirstHits() is such a
  // culprit). Use the global positions of the states and the desired
  // refit direction to find the starting TSOS.
  TrajectoryStateOnSurface firstTSOS, lastTSOS;
  unsigned int innerId, outerId;
  bool order_swapped = track.outermostMeasurementState().globalPosition().mag() < track.innermostMeasurementState().globalPosition().mag();
  bool inner_is_first;
  LogTrace(theCategory) << "order swapped? " << order_swapped;

  // Fill the starting state, depending on the ordering above.
  if ((theRefitDirection == RinsideOut && !order_swapped) || (theRefitDirection == RoutsideIn && order_swapped)) {
    innerId   = newTrack.innerDetId();
    outerId   = newTrack.outerDetId();
    firstTSOS = track.innermostMeasurementState();
    lastTSOS  = track.outermostMeasurementState();
    inner_is_first = true;
  }
  else {
    innerId   = newTrack.outerDetId();
    outerId   = newTrack.innerDetId();
    firstTSOS = track.outermostMeasurementState();
    lastTSOS  = track.innermostMeasurementState();
    inner_is_first = false;
  } 

  LogTrace(theCategory) << "firstTSOS: inner_is_first? " << inner_is_first
			<< " globalPosition is " << firstTSOS.globalPosition()
			<< " innerId is " << innerId;

  if(!firstTSOS.isValid()){
    LogWarning(theCategory) << "Error wrong initial state!" << endl;
    return vector<Trajectory>();
  }

  firstTSOS.rescaleError(1000.);

  

  // This is the only way to get a TrajectorySeed with settable propagation direction
  PTrajectoryStateOnDet garbage1;
  edm::OwnVector<TrackingRecHit> garbage2;
  PropagationDirection propDir = 
    (firstTSOS.globalPosition().basicVector().dot(firstTSOS.globalMomentum().basicVector())>0) ? alongMomentum : oppositeToMomentum;

  // These lines cause the code to ignore completely what was set
  // above, and force propDir for tracks from collisions!
//  if(propDir == alongMomentum && theRefitDirection == RoutsideIn)  propDir=oppositeToMomentum;
//  if(propDir == oppositeToMomentum && theRefitDirection == RinsideOut) propDir=alongMomentum;

  const TrajectoryStateOnSurface& tsosForDir = inner_is_first ? lastTSOS : firstTSOS;
  propDir = (tsosForDir.globalPosition().basicVector().dot(tsosForDir.globalMomentum().basicVector())>0) ? alongMomentum : oppositeToMomentum;
  LogTrace(theCategory) << "propDir based on firstTSOS x dot p is " << propDir
			<< " (alongMomentum == " << alongMomentum << ", oppositeToMomentum == " << oppositeToMomentum << ")";

  // Additional propagation diretcion determination logic for cosmic muons
  if (theCosmicFlag) {
    PropagationDirection propDir_first = (firstTSOS.globalPosition().basicVector().dot(firstTSOS.globalMomentum().basicVector()) > 0) ? alongMomentum : oppositeToMomentum;
    PropagationDirection propDir_last  = (lastTSOS .globalPosition().basicVector().dot(lastTSOS .globalMomentum().basicVector()) > 0) ? alongMomentum : oppositeToMomentum;
    LogTrace(theCategory) << "propDir_first " << propDir_first << ", propdir_last " << propDir_last
			  << " : they " << (propDir_first == propDir_last ? "agree" : "disagree");

    int y_count = 0;
    for (TransientTrackingRecHit::ConstRecHitContainer::const_iterator it = recHitsForReFit.begin(); it != recHitsForReFit.end(); ++it) {
      if ((*it)->globalPosition().y() > 0) ++y_count;
      else --y_count;
    }
    
    PropagationDirection propDir_ycount = alongMomentum;
    if (y_count > 0) {
      if      (theRefitDirection == RinsideOut) propDir_ycount = oppositeToMomentum;
      else if (theRefitDirection == RoutsideIn) propDir_ycount = alongMomentum;
    }
    else {
      if      (theRefitDirection == RinsideOut) propDir_ycount = alongMomentum;
      else if (theRefitDirection == RoutsideIn) propDir_ycount = oppositeToMomentum;
    }
    
    LogTrace(theCategory) << "y_count = " << y_count
			  << "; based on geometrically-outermost TSOS, propDir is " << propDir << ": "
			  << (propDir == propDir_ycount ? "agrees" : "disagrees")
			  << " with ycount determination";
    
    if (propDir_first != propDir_last) {
      LogTrace(theCategory) << "since first/last disagreed, using y_count propDir";
      propDir = propDir_ycount;
    }
  }


  TrajectorySeed seed(garbage1,garbage2,propDir);
  ConstRecHitContainer newseedtrackerhits; 
  TrajectorySeed* newseed = 0 ;
  vector<Trajectory> trajectories;
  TrajectoryStateOnSurface badgerTSOS = 0; 
  //regenerate the seed only if requested, otherwise do the usual thing
  if(theNewSeed){
      newseed= NewSeedFromPairOrTriplet(recHitsForReFit,propDir,newseedtrackerhits);
      //If the new seed was created, try to extract the TSOS
      if (newseed != 0) {
          badgerTSOS = 
              TSOSFromNewSeed(newseed,newseedtrackerhits,&*theService->magneticField().product());
      }  
      else{
          LogDebug(theCategory) <<"There is no new seed or it is invalid"
                                <<endl;
          return vector<Trajectory>();
      }
  }//-----if theNewSeed
  else{
      newseed = &seed;
      if(recHitsForReFit.front()->geographicalId() != DetId(innerId)){
          LogDebug(theCategory)<<"Propagation occured"<<endl;
          LogTrace(theCategory) << "propagating firstTSOS at " << firstTSOS.globalPosition()
                                << " to first rechit with surface pos " << recHitsForReFit.front()->det()->surface().toGlobal(LocalPoint(0,0,0));
          firstTSOS = theService->propagator(thePropagatorName)->propagate(firstTSOS, recHitsForReFit.front()->det()->surface());
          if(!firstTSOS.isValid()){
      LogDebug(theCategory)<<"Propagation error!"<<endl;
      return vector<Trajectory>();
          }
      }

      LogDebug(theCategory) <<"Passed the old firstTSOS"<<endl;
      
      LogDebug(theCategory) << " Old GlobalMuonRefitter : theFitter " 
                            << propDir << endl;
      LogDebug(theCategory) << "                      Old First TSOS: " 
       << firstTSOS.globalPosition() << "  p="
                            << firstTSOS.globalMomentum() << " = "
                            << firstTSOS.globalMomentum().mag() << endl;
      
      badgerTSOS = scaleTSOS(firstTSOS, theStabScale);
      //badgerTSOS = firstTSOS;

  }//----------else if theNewSeed

  //build the trajectory 
  trajectories = theFitter->fit(*newseed,recHitsForReFit,badgerTSOS);
  LogDebug(theCategory) <<"Trajectory fitter called successfully"<<endl;
  
  
  if(trajectories.empty()){
    LogDebug(theCategory) << "No Track refitted!" << endl;
    return vector<Trajectory>();
  }
  
  return trajectories;
}


//
// Remove Selected Station Rec Hits
//
GlobalMuonRefitter::ConstRecHitContainer GlobalMuonRefitter::getRidOfSelectStationHits(ConstRecHitContainer hits) const
{
  ConstRecHitContainer results;
  ConstRecHitContainer::const_iterator it = hits.begin();
  for (; it!=hits.end(); it++) {

    DetId id = (*it)->geographicalId();
//	if (id.det() != DetId::Tracker) continue;
    //Check that this is a Muon hit that we're toying with -- else pass on this because the hacker is a moron / not careful

    if (id.det() == DetId::Tracker && theTrackerSkipSystem > 0) {
      int layer = -999;
      int disk  = -999;
      int wheel = -999;
      if ( id.subdetId() == theTrackerSkipSystem){
	//                              continue;  //caveat that just removes the whole system from refitting

	if (theTrackerSkipSystem == PXB) {
	  PXBDetId did(id.rawId());
	  layer = did.layer();
	}
	if (theTrackerSkipSystem == TIB) {
	  TIBDetId did(id.rawId());
	  layer = did.layer();
	}

	if (theTrackerSkipSystem == TOB) {
	  TOBDetId did(id.rawId());
	  layer = did.layer();
	}
	if (theTrackerSkipSystem == PXF) {
	  PXFDetId did(id.rawId());
	  disk = did.disk();
	}
	if (theTrackerSkipSystem == TID) {
	  TIDDetId did(id.rawId());
	  wheel = did.wheel();
	}
	if (theTrackerSkipSystem == TEC) {
	  TECDetId did(id.rawId());
	  wheel = did.wheel();
	}
	if (theTrackerSkipSection >= 0 && layer == theTrackerSkipSection) continue;
	if (theTrackerSkipSection >= 0 && disk == theTrackerSkipSection) continue;
	if (theTrackerSkipSection >= 0 && wheel == theTrackerSkipSection) continue;
      }
    }

    if (id.det() == DetId::Muon && theSkipStation) {
      int station = -999;
      int wheel = -999;
      if ( id.subdetId() == MuonSubdetId::DT ) {
	DTChamberId did(id.rawId());
	station = did.station();
	wheel = did.wheel();
      } else if ( id.subdetId() == MuonSubdetId::CSC ) {
	CSCDetId did(id.rawId());
	station = did.station();
      } else if ( id.subdetId() == MuonSubdetId::RPC ) {
	RPCDetId rpcid(id.rawId());
	station = rpcid.station();
      }
      if(station == theSkipStation) continue;
    }
    results.push_back(*it);
  }
  return results;
}

GlobalMuonRefitter::ConstRecHitContainer GlobalMuonRefitter::getEvenOddLayers(ConstRecHitContainer hits, bool const getEven) const{

  ConstRecHitContainer results;
  ConstRecHitContainer::const_iterator it = hits.begin();
  int ihit =0;
  for (; it!=hits.end(); it++) {
      
      int iLayer = -1;
      DetId id = (*it)->geographicalId();
      bool keepHit = true;
      if (id.det() == DetId::Tracker) {
          int subdetid = id.subdetId();
          if (subdetid == PXB) iLayer = PXBDetId(id.rawId()).layer();
          if (subdetid == TIB) iLayer = TIBDetId(id.rawId()).layer();
          if (subdetid == TOB) iLayer = TOBDetId(id.rawId()).layer();
          if (subdetid == PXF) iLayer = PXFDetId(id.rawId()).disk();
          if (subdetid == TID) iLayer = TIDDetId(id.rawId()).wheel();
          if (subdetid == TEC) iLayer = TECDetId(id.rawId()).wheel();
      }

      LogTrace(theCategory)<<"GetEvenOddLayers: Tracker iLayer = "<<iLayer<<endl;
      
      if (iLayer >= 0){
          if (getEven) {
              keepHit &= (iLayer%2==0);
          }else {
              keepHit &= (iLayer%2!=0);
          }
      }//tracker layer found
      else if (iLayer < 0){ //keep muon even/odd hits as well
          if (!(id.det() == DetId::Muon)) continue;
          if (getEven) {
              keepHit &= ((ihit%2) == 0);
                  } 
          else{ 
              keepHit &= ((ihit%2) != 0);
          }
          ++ihit;
      }//check if its a muon hit
    
      
      if (keepHit) results.push_back(*it);
  }//loop over hits

  return results;
}



GlobalMuonRefitter::ConstRecHitContainer GlobalMuonRefitter::getEvenOddHits(ConstRecHitContainer hits, bool const getEven) const{

	int ihit = 0;
  ConstRecHitContainer results;
  ConstRecHitContainer::const_iterator it = hits.begin();

//	++it;
//	ihit = 1;

  for (; it!=hits.end(); it++) {
    	
//    DetId id = (*it)->geographicalId();
//	if (id.det() != DetId::Tracker) continue;

	if (getEven && (ihit%2) ==0) {
		results.push_back(*it);
	} else{ 
   		if (!getEven && (ihit%2) !=0) results.push_back(*it);
	}
	++ihit;
 }
	return results;
}
TrajectoryStateOnSurface GlobalMuonRefitter::scaleTSOS(TrajectoryStateOnSurface tsosIn, double const scale) const{


	GlobalTrajectoryParameters params = tsosIn.globalParameters(); 

	TRandom3 rand3jit(0);
/*
	double jit = TRandom3(0).Gaus(0,0.1);	
//	jit+= 2.0;
	jit+=scale;
//	jit+= 9999.0;
//	jit=0.;
*/
/*
	GlobalPoint glbPoint(
	params.position().x()*jit, 	
	params.position().y()*jit, 	
	params.position().z()*jit);	
*/

//	double jitval = rand3jit.Gaus(scale,0.1);

	GlobalVector glbVec(
		params.momentum().x()*scale*rand3jit.Gaus(1.,0.1),	//*scale*rand3jit.Gaus(0,0.1), 	
		params.momentum().y()*scale*rand3jit.Gaus(1.,0.1),	//*scale*rand3jit.Gaus(0,0.1), 	
		params.momentum().z()*scale*rand3jit.Gaus(1.,0.1));	//*scale*rand3jit.Gaus(0,0.1));	

	GlobalTrajectoryParameters badger(
		params.position(),
		glbVec,
		params.charge(),
		&(*(theService->magneticField()))
//		*mf
//		params.magneticField()		
	); 

	TrajectoryStateOnSurface badgerTSOS(badger,tsosIn.curvilinearError(), tsosIn.surface(), tsosIn.surfaceSide());

	return badgerTSOS;
}




