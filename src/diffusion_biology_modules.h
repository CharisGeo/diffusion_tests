#ifndef DIFFUSION_BIOLOGY_MODULES_H_
#define DIFFUSION_BIOLOGY_MODULES_H_

#include "biodynamo.h"

namespace bdm {

// List the extracellular substances
enum Substances { kSubstance, kYsubstance, kZsubstance };

// Define displacement behavior:
// Cells move along the diffusion gradient (from low concentration to high)
struct Chemotaxis : public BaseBiologyModule {
  Chemotaxis() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    auto* sim = TSimulation::GetActive();
    auto* rm = sim->GetResourceManager();
    static auto* dg = rm->GetDiffusionGrid(kSubstance);
    dg->SetConcentrationThreshold(1e15);
    static auto* yDg = rm->GetDiffusionGrid(kYsubstance);
    yDg->SetConcentrationThreshold(1e15);
    static auto* zDg = rm->GetDiffusionGrid(kZsubstance);
    zDg->SetConcentrationThreshold(1e15);

    auto& position = cell->GetPosition();
    std::array<double, 3> gradient;
    dg->GetGradient(position, &gradient);
    gradient[0] *= 0.5;
    gradient[1] *= 0.5;
    gradient[2] *= 0.5;

    auto& place = cell->GetPosition();
    std::array<double, 3> gradient2;
    yDg->GetGradient(place, &gradient2);
    gradient2[0] *= 1.5;
    gradient2[1] *= 1.5;
    gradient2[2] *= 1.5;

    auto& wall = cell->GetPosition();
    std::array<double, 3> gradient3;
    zDg->GetGradient(wall, &gradient3);
    gradient3[0] *= 1.5;
    gradient3[1] *= 1.5;
    gradient3[2] *= 1.5;


    cell->UpdatePosition(gradient);
    cell->UpdatePosition(gradient2);
    cell->UpdatePosition(gradient3);
  }

  ClassDefNV(Chemotaxis, 1);
};


}  // namespace bdm

#endif  // DIFFUSION_BIOLOGY_MODULES_H_
