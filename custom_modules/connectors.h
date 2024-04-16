
#ifndef __CONNECTORS_H__
#define __CONNECTORS_H__

#include <memory>
#include <stdexcept>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include<unordered_map>


#include "./follicle_utilities.h"
#include "./springs.h"
#include <vector>
using namespace BioFVM;
using namespace PhysiCell;
namespace Connectors{
// eventually Connector should use template args to handle any kind of connector
class Connector
{
  private:
  //
  public:
  enum connector_type{simple_spring, crosslink};
  // crosslink 0 length attachment point
  double m_rest_length;
  Cell* m_pNeighbor;
  connector_type m_type;
  std::vector<double> m_attachment_point;
  std::vector<double> m_termination_point; 
  Connector* m_my_ptr;
  int pCell_index;
  Connector(Cell* pNeighbor, double length, connector_type type); 
  ~Connector();
  const connector_type& get_connector_type();
  const Connector& get_connector();
  const Cell* get_neighbor();
  const double& get_rest_length();
  void set_attachment_point(std::vector<double> point);
  void set_termination_point(std::vector<double> point);
};
void delete_connector(Connector& connector);

struct Connector_List{
  Cell* member_pCell;
  std::vector<Connector*> connection_list;
  const Connector_List& get_connector_list();//return *this;
  Connector_List(Cell* m_pCell);
  ~Connector_List();
  void add_connection(Connector* pCon);
  void remove_connection(Connector* pCon);
  bool has_connector_list(Cell* pCell);
};
bool is_valid_pCell();
void remove_all_connections();
void invalid_pCell_cleanup();
// function to create connection list, attached to object and store it's own pointer

extern std::vector <Connector_List*> connector_lists_by_pCell_index;
};


#endif // __CONNECTORS_H__
