/*
 * mqPMTSD.cc
 *
 *  Created on: 22-4-2019
 *      Author: schmitz
 */
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "mqPMTSD.hh"
#include "mqPMTHit.hh"
#include "mqDetectorConstruction.hh"

//#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include <math.h>
#include <string>
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
mqPMTSD::mqPMTSD(G4String name)
  :G4VSensitiveDetector(name),pmtHitCollection(0),
//	pmtPositionsX(0),
//	pmtPositionsY(0),
//	pmtPositionsZ(0),
	verbose(0)
{

	collectionName.insert("pmtHitCollection");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
mqPMTSD::~mqPMTSD()

{

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void mqPMTSD::Initialize(G4HCofThisEvent* HCE){
  pmtHitCollection = new mqPMTHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //Store collection with event and keep ID
  static G4int HCID = -1;
  if(HCID<0){
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection( HCID, pmtHitCollection );
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4bool mqPMTSD::ProcessHits(G4Step* ,G4TouchableHistory* ){
  return false;
}

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4bool mqPMTSD::ProcessHits_constStep(const G4Step* aStep,
				       G4TouchableHistory* ){

  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

  //this is the copy number within a layer
  G4int pmtCopyNumber =
    aStep->GetPostStepPoint()->GetTouchable()->GetCopyNumber();

  //extract layer number of hit
  char layerNumberChar = aStep->GetPostStepPoint()->GetTouchable()->GetVolume()->GetName().back();
  //convert layer number to char (ASCII characters work this way)
  G4int layerNumber = layerNumberChar-48;


  //output the results (you can check that this works if you're curious)
  //G4cout << "Layer number, char: " << layerNumberChar << " int: " << layerNumber << G4endl;
/* 
  //get the detector to retrieve layer information
  mqDetectorConstruction* detector = (mqDetectorConstruction*)G4RunManager::GetRunManager()
					->GetUserDetectorConstruction();

  //get number of bars per layer
  G4int nBarPerLayer = detector->GetNBarPerLayer();
*/
  G4int nBarPerLayer = 6;
  //this is the copy number, independent of layer
  //this ID number is using my numbering convention, 0 is bottom left, 1 is middle left, 2 is upper left, 3 is bottom right, 4 is middle right, 5 is upper right. Add 0 for front layer, 6 for mid layer, 12 for back layer
  G4int pmtNumber=0;
  if(pmtCopyNumber<6) pmtNumber = nBarPerLayer*layerNumber+pmtCopyNumber;
  else pmtNumber = pmtCopyNumber;
  
  //R878 PMTs: 0,2,5,6,7,8,9,11,12,14,16,17
  //R7725 PMTs: 13,15
  //ET9814B PMTs: 1,3,4,10

  G4double particleEnergy=aStep->GetTrack()->GetTotalEnergy();
  G4double particleWavelength=(1240/particleEnergy)*eV; //wavelength in nm
  //G4cout << "particle wavelength is: " << particleWavelength << G4endl;
  G4double detProb=0;
  G4bool det=false;
  
  if(pmtNumber==1|pmtNumber==3|pmtNumber==4|pmtNumber==10){
  	detProb = GetET9814B_QE().Value(particleWavelength);
//	G4cout << "Hit ET9814B, detProb is: " << detProb << G4endl;
  } else if(pmtNumber==13|pmtNumber==15){
  	detProb = GetR7725_QE().Value(particleWavelength);
//	G4cout << "Hit R7725, detProb is: " << detProb << G4endl;
  } else {
  	detProb = GetR878_QE().Value(particleWavelength);
//	G4cout << "Hit R878, detProb is: " << detProb << G4endl;
  }	

	det=(G4UniformRand()<detProb);
  if(det){
//	G4cout << "successfully detected in PMT: " << pmtNumber << G4endl;
//Get information from the optical photon which hit the pmt

  G4ThreeVector myEndDirection = aStep->GetPostStepPoint()->GetMomentumDirection();
  G4double theta = acos(std::abs(myEndDirection.dot(G4ThreeVector(0,0,1))));
  G4double energyDeposit = aStep->GetTotalEnergyDeposit();
  G4double globalTime = aStep->GetPostStepPoint()->GetGlobalTime(); //time since the event was created

  G4ThreeVector hitPos = aStep->GetPostStepPoint()->GetPosition();
  G4double xy_mag = sqrt(pow(hitPos.x(),2)+pow(hitPos.y(),2));
  G4int pid = aStep->GetTrack()->GetParentID();

  //Find the correct hit collection

  G4int n = pmtHitCollection->entries();
  //mqPMTHit* hit=NULL;
  //for(G4int i=0; i < n; i++){
  //  if((*pmtHitCollection)[i]->GetPMTNumber() == pmtNumber){
  //    hit=(*pmtHitCollection)[i];
  //    break;
  //  }
  //}

  //if(hit==NULL){//this pmt wasnt previously hit in this event
  //
	mqPMTHit*
	hit = new mqPMTHit(); //so create new hit

    hit->SetPMTNumber(pmtNumber);
    hit->SetHitAngle(theta);
    //hit->SetPMTPhysVol(physVol);
    hit->SetParentID(pid);
    hit->SetTrackID(aStep->GetTrack()->GetTrackID());
    hit->SetInitialEDep(energyDeposit); // TODO this is EDep of first hit
    //hit->IncEDep(energyDeposit);

    hit->SetScintToPMT(true);
    hit->SetHitTime(globalTime);// this is hit time of first hit
    hit->SetHitPosition(hitPos);//this is only hit position of first hit
    //hit->SetTrackLength(aStep->GetTrack()->GetTrackLength());
    pmtHitCollection->insert(hit);
    //hit->SetPMTPos((*pmtPositionsX)[pmtNumber],(*pmtPositionsY)[pmtNumber],
	//	   (*pmtPositionsZ)[pmtNumber]);

  // }
  //increment hit for the selected pmt
  //
    //hit->IncEDep(energyDeposit);
  	//hit->SetLastHitTime(globalTime);
  	//hit->IncPhotonCount();
//////////////////////////
//G4cout << "X-Position: " << hitPos.getX() << G4endl;
//G4cout << "Y-Position: " << hitPos.getY() << G4endl;
//G4cout << "Z-Position: " << hitPos.getZ() << G4endl;
/////////////////////////
 } //else {G4cout << "not detected" << G4endl;}
  return true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void mqPMTSD::EndOfEvent(G4HCofThisEvent* ){
	G4String collName = pmtHitCollection->GetName() ;
	G4int NbHits = pmtHitCollection->entries();
	if (verbose > 0){
			G4cout << "\n-------->Hits Collection "<< collName << ": in this event they are " << NbHits
					<< " hits in the pmts: " << G4endl;

			for (G4int i=0;i<NbHits;i++){
				(*pmtHitCollection)[i]->Print();
				//G4cout <<"Hit time : " << (*pmtHitCollection)[i]->GetHitTime() << G4endl;;

			}

	}



}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void mqPMTSD::clear(){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void mqPMTSD::DrawAll(){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void mqPMTSD::PrintAll(){
}



