#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

// Definition of a Fractionally Charged Particle (FCP)
// to be used for cosmic FCP studies

void DefineFCP()
{
    const G4String name = "fcp";

    // If already in table: do nothing
    if (G4ParticleTable::GetParticleTable()->FindParticle(name)) {
        return;
    }

    const G4double mass   = 0.1 * GeV;
    const G4double charge = 0.1 * eplus;
    const G4double width  = 0.;
    const G4double tau    = -1.0;
    const G4int    pdg    = 100003;

    new G4ParticleDefinition(
        name,
        mass,
        width,
        charge,
        1, 0, 0,
        0, 0, 0,
        "lepton", 0, 0,
        pdg,
        true,
        tau,
        nullptr
    );
}
