#ifndef DEMO_DIFFUSION_MODULE_H_
#define DEMO_DIFFUSION_MODULE_H_
#ifndef INTEGRATION_SUBSTANCE_INITIALIZATION_H_
#define INTEGRATION_SUBSTANCE_INITIALIZATION_H_

#include <vector>
#include <string>

#include "biodynamo.h"
#include "my_cell.h"
#include "diffusion_biology_modules.h"
#include "substance_initializers.h"

namespace bdm {

using std::array;
using std::vector;
using std::string;

// -----------------------------------------------------------------------------
// This model creates 8 cells at each corner of a cube, and one in the middle.
// The cell in the middle secretes a substance. The cells are modeled to
// move according to the extracellular gradient; in this case to the middle.
// -----------------------------------------------------------------------------

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<Chemotaxis, RegulateGenes, GrowDivide>;
  using AtomicTypes = VariadicTypedef<MyCell>;
};

inline int Simulate(int argc, const char** argv) {
  // Initialize BioDynaMo
  Simulation<> simulation(argc, argv);
  auto* param = simulation.GetParam();

  // 3. Define initial model
  // Create an artificial bounds for the simulation space
  param->bound_space_ = true;
  param->min_bound_ = 0;
  param->max_bound_ = 1000;
  param->run_mechanical_interactions_ = false;

  RegulateGenes regulate_example;
  regulate_example.AddGene(
      [](double curr_time, double last_concentration) {
        return curr_time * last_concentration + 0.2f;
      },
      1);
  regulate_example.AddGene(
      [](double curr_time, double last_concentration) {
        return last_concentration * last_concentration * curr_time;
      },
      5);
  regulate_example.AddGene(
      [](double curr_time, double last_concentration) {
        return last_concentration + curr_time + 3;
      },
      7);

  auto create = [&](const std::array<double, 3>& place){
   MyCell cell(place);
   cell.SetDiameter(30);
   cell.SetAdherence(0.4);
   cell.SetMass(1.0);
   cell.SetCellType(1);
   cell.AddBiologyModule(regulate_example);
   cell.AddBiologyModule(Chemotaxis());
   cell.AddBiologyModule(GrowDivide(35, 3000, {gAllBmEvents}));
   return cell;
   };
 ModelInitializer::CreateCellsRandom(param->min_bound_, param->max_bound_, 100, create);

  auto construct = [&](const std::array<double, 3>& position) {
    MyCell cell(position);
    cell.SetDiameter(50);
    cell.SetMass(1.0);
    cell.SetAdherence(0.4);
    cell.SetCellType(-1);
    cell.AddBiologyModule(regulate_example);
    cell.AddBiologyModule(Chemotaxis());
    cell.AddBiologyModule(GrowDivide(52, 3000, {gAllBmEvents}));
    return cell;
  };
 ModelInitializer::CreateCellsRandom(param->min_bound_, param->max_bound_, 100, construct);
 

  // 3. Define the substances in our simulation
  // Order: substance id, substance_name, diffusion_coefficient, decay_constant,
  // resolution
  ModelInitializer::DefineSubstance(kSubstance, "Substance", 0.5, 0.005, 10);
  ModelInitializer::DefineSubstance(kYsubstance, "Ysubstance", 0.5, 0.005, 10);
  ModelInitializer::DefineSubstance(kZsubstance, "Zsubstance", 0.5, 0.005, 10);

  // Order: substance id, substance name, initialization model, along which axis
  // (0 = x, 1 = y, 2 = z). See the documentation of `GaussianBand` for
  // information about its arguments
  ModelInitializer::InitializeSubstance(kSubstance, "Substance",
                                        GaussianBand(50, 250, Axis::kXAxis));
  ModelInitializer::InitializeSubstance(kSubstance, "Substance",
                                        GaussianBand(50, 250, Axis::kYAxis));
  ModelInitializer::InitializeSubstance(kSubstance, "Substance",
                                        GaussianBand(50, 250, Axis::kZAxis));
  ModelInitializer::InitializeSubstance(kYsubstance, "Ysubstance",
                                        GaussianBand(200, 250, Axis::kXAxis));
  ModelInitializer::InitializeSubstance(kYsubstance, "Ysubstance",
                                        GaussianBand(200, 250, Axis::kYAxis));
  ModelInitializer::InitializeSubstance(kYsubstance, "Ysubstance",
                                        GaussianBand(200, 250, Axis::kZAxis));
  ModelInitializer::InitializeSubstance(kZsubstance, "Zsubstance",
                                        GaussianBand(400, 250, Axis::kXAxis));
  ModelInitializer::InitializeSubstance(kZsubstance, "Zsubstance",
                                        GaussianBand(400, 250, Axis::kYAxis));
  ModelInitializer::InitializeSubstance(kZsubstance, "Zsubstance",
                                        GaussianBand(400, 250, Axis::kZAxis));

  // Run simulation for N timesteps
  auto* scheduler = simulation.GetScheduler();
  scheduler->Simulate(2000);
  
  auto* rm = simulation.GetResourceManager();
  auto&& cell = (*(rm->Get<MyCell>()))[0];
  const auto* regulate_genes = cell.GetBiologyModules<RegulateGenes>()[0];
  const auto& concentrations = regulate_genes->GetConcentrations();
  std::cout << "Gene concentrations after " << scheduler->GetSimulatedSteps()
            << " time steps" << std::endl;
  for (double concentration : concentrations) {
    std::cout << concentration << std::endl;
   }

  std::cout << "Simulation completed successfully!\n";
  return 0;
  }
}  // namespace bdm

#endif  // DEMO_DIFFUSION_MODULE_H_
#endif 
