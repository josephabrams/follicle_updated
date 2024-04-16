#include "./connectors.h"
#include <iterator>
constexpr bool DEBUG_CONNECTOR{true};
std::string connector_debug_file = " ../debug_logs/connector.out";
namespace Connectors{

std::vector <Connector_List*> connector_lists_by_pCell_index;

Connector::Connector(Cell* pNeighbor, double length, connector_type type)
: m_rest_length{length}, m_pNeighbor{pNeighbor}, m_type{type}, m_attachment_point{0.0, 0.0, 0.0},
  m_termination_point{0.0, 0.0, 0.0}, m_my_ptr{this}
{ 
  if(DEBUG_CONNECTOR==true) {

    std::ofstream ofs;
    ofs.open (connector_debug_file, std::ofstream::out | std::ofstream::app);
    std::cout<< "Connector created and attached to "<< pNeighbor<<"\n"; 
  }
}

Connector::~Connector()
{}

const Connector::connector_type& Connector::get_connector_type()
{ return this->m_type; }

const Cell* Connector::get_neighbor()
{ return this->m_pNeighbor; }

const double& Connector::get_rest_length()
{return this->m_rest_length;}
const Connector& Connector::get_connector()
{ return *this;}

void Connector::set_attachment_point(std::vector <double> point)
{this->m_attachment_point= point;}
void Connector::set_termination_point(std::vector <double> point)
{this->m_termination_point= point;}

void create_simple_spring(Cell* pCell, Cell* pNeighbor)
{
  if(true)
  {
    double membrane_distance=pCell->phenotype.geometry.radius+pNeighbor->phenotype.geometry.radius;
    double initial_length= norm(pCell->position-pNeighbor->position)-membrane_distance;
    Connector(pNeighbor,initial_length,Connector::connector_type::simple_spring);
  }
}
void delete_connector(Connector& connector)
{
  Connector* connector_ptr{connector.m_my_ptr};
  delete[] connector_ptr;
}

Connector_List::Connector_List(Cell* m_pCell)
: member_pCell{m_pCell}
{
  //TODO: resize connector list to all cells and then just specify index
  // this means resizing the list when new agents are made/removed
  //
  connection_list.clear();
  std::vector<Connector_List*>::iterator it;
  it = connector_lists_by_pCell_index.begin();
  connector_lists_by_pCell_index.insert ( it+m_pCell->index , this );
}

const Connector_List& Connector_List::get_connector_list()
{return *this;}//return *this;
Connector_List::~Connector_List()
{
      
    // Using iterator to handle empty vector case
  std::vector<Connector*>::iterator it;
  int index=0;
  for (it = connection_list.begin(); it < connection_list.end(); it++) 
    { 
       Connector* connector= connection_list[index];
       delete_connector(*connector);
       index++;
    }
}

bool Connector_List::has_connector_list(Cell* pCell){
  if(connector_lists_by_pCell_index[pCell->index]->member_pCell==pCell)
  {
    return true;
  }
  else {
    return false;
  }
}
void Connector_List::add_connection(Connector* pCon){}
void Connector_List::remove_connection(Connector* pCon){}
bool is_valid_pCell(){return true;}
void remove_all_connections(){}
void invalid_pCell_cleanup(){}



};
