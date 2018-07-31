#ifndef MY_CELL_H_
#define MY_CELL_H_

#include <biodynamo.h>

namespace bdm {

// Define my custom cell, which extends Cell by adding an extra
// data member cell_type.
BDM_SIM_OBJECT(MyCell, Cell) {
  BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cell_type_);

 public:
  MyCellExt() {}
  explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {}

  void SetCellType(int t) { cell_type_[kIdx] = t; }
  int GetCellType() const { return cell_type_[kIdx]; }

 private:
  vec<int> cell_type_;
};

}  // namespace bdm

#endif // MY_CELL_H_
